#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include "data_structure.h"
#include "pugixml.hpp"
#include "screen.h"
#include "utility_functions.h"
//#define TIXML_USE_STL

using namespace std;
using namespace pugi;

//typedef TiXmlElement tag;
int str2int(string ss){
   stringstream strValue;
   strValue << ss;
   int intValue;
   strValue >> intValue;
   return intValue;
}


bool update_3d_vector(vector<threevector>* objectToUpdate, float x, float y, float z) {
   threevector r;
   r.push_back(x);r.push_back(y);r.push_back(z);
   (*objectToUpdate).push_back(r); //[0] = x;  (*objectToUpdate)[1] = y; (*objectToUpdate)[2] = z;
   return true;
}

int readXML(VasprunXML *vasprun, Configuration *config) {
   float x, y, z;
   vector<float> r;
   r.push_back(0);
   r.push_back(0);
   r.push_back(0);

   xml_document doc;
   screen.status << "Loading file into memory";
   xml_parse_result result = doc.load_file(vasprun->input_filename.c_str() );
   if (result) {
      screen.status << "File loaded";
   } else {
      screen.error << "WARNING: Parsed with errors.";
      screen.error << result.description();
   }
   
  
   screen.status << "Beginning to parse"   ;

   xml_node n_incar = doc.child("modeling").child("incar");
   key thisKey;
   for (xml_node i = n_incar.child("i"); i; i=i.next_sibling("i")) {
      thisKey.name = i.attribute("name").value();
      thisKey.value = i.child_value(); 
      vasprun->incar.keys.push_back(thisKey);

      /*copy some parameters into special objects*/
      if (0==strcmp("POTIM",thisKey.name.c_str())) {
         vasprun->dt = stod(thisKey.value);
      } else if (0==strcmp("SYSTEM",thisKey.name.c_str())) {
         vasprun->system_name = thisKey.value;
      } else if (0==strcmp("TEBEG",thisKey.name.c_str())) {
         vasprun->starting_temperature = stod(thisKey.value);
      }
   }
   


   xml_node n_atominfo = doc.child("modeling").child("atominfo");
   vasprun->numatoms = str2int(n_atominfo.child("atoms").child_value());
   vasprun->numtypes = str2int(n_atominfo.child("types").child_value());

   //<array name="atomtypes">  
   //get the first <array> tag
   xml_node n_temp = n_atominfo.select_nodes("//atominfo/array[@name='atomtypes']")[0].node() ;
   //from the <array> tag, select .//set/rc
   xpath_node_set ns_atomtypes = n_temp.select_nodes(".//set/rc");
   

   //iterate over each <rc> tag and and put the contents (<c></c>) into objects
   for (xpath_node_set::const_iterator it = ns_atomtypes.begin(); it != ns_atomtypes.end(); ++it) {
      xpath_node_set ns_c = it->node().select_nodes(".//c");
      atomType atomTypeData;
      atomTypeData.atomspertype    = str2int(it->node().child("c").child_value() );
      atomTypeData.element         = trim(it->node().child("c").next_sibling().child_value());
      atomTypeData.mass            = stod(it->node().child("c")\
                                   .next_sibling().next_sibling().child_value() );
      vasprun->atomic_mass[atomTypeData.element] = atomTypeData.mass;
      atomTypeData.valence         = stod(it->node().child("c").next_sibling()\
                                   .next_sibling().next_sibling().child_value() ); 
      atomTypeData.pseudopotential = it->node().child("c").next_sibling().next_sibling()\
                                   .next_sibling().next_sibling().child_value() ;
//      cout << "TEST~" << endl;
      vasprun->atoms.push_back(atomTypeData);   
   }

   //select the 'atoms' array (containing ion symbols)
   n_temp = n_atominfo.select_nodes("//atominfo/array[@name='atoms']")[0].node();
   //from the <array> tag, select .//set/rc
   xpath_node_set ns_ionlist = n_temp.select_nodes(".//set/rc");
   int counter=0;
   string symbol;
   for (xpath_node_set::const_iterator it = ns_ionlist.begin(); it != ns_ionlist.end(); ++it) {
      symbol=it->node().child("c").child_value();
      vasprun->ion_symbols.push_back( ltrim(rtrim(symbol)) );
//      cout << counter  << "    " <<  symbol <<  endl;
      counter++;
   }




   //Get the lattice vectors
   n_temp = doc.select_nodes("//modeling/structure[@name='initialpos']/crystal/varray[@name='basis']")[0].node();
   sscanf( n_temp.child("v").child_value() , "  %f  %f  %f " , &x, &y, &z );
   update_3d_vector(&vasprun->latt, x,y,z);
   sscanf( n_temp.child("v").next_sibling().child_value() , "  %f  %f  %f " , &x, &y, &z );
   update_3d_vector(&vasprun->latt, x,y,z);
   sscanf( n_temp.child("v").next_sibling().next_sibling().child_value() , "  %f  %f  %f " , &x, &y, &z );
   update_3d_vector(&vasprun->latt, x,y,z);


   xpath_node_set ns_temp;
   //Iterate over each timestep
   xpath_node_set ns_timesteps = doc.select_nodes("//modeling/calculation");

   int timecounter = 0;
   //for each timestep
   for (xpath_node_set::const_iterator it = ns_timesteps.begin(); it != ns_timesteps.end(); ++it) {
      TimeStep thisStep;
      thisStep.realTime = timecounter;
      timecounter++;
      //get stress tensor:
      ns_temp = it->node().select_nodes(".//varray[@name='stress']/v");
      vector<string> this_row;
      for (int i =0; i<3; i++) {
         this_row = str2vec(ns_temp[i].node().child_value()," ");
         thisStep.stress_tensor[0][i] = stod(this_row[0]);
         thisStep.stress_tensor[1][i] = stod(this_row[1]);
         thisStep.stress_tensor[2][i] = stod(this_row[2]);         
      }

      //get energies:
      ns_temp = it->node().select_nodes("./energy/i");
      for (int i=0; i<ns_temp.size(); i++) {
         thisStep.energy[ns_temp[i].node().attribute("name").value()] = stod(ns_temp[i].node().child_value());
      }
      


      ns_temp = it->node().select_nodes(".//structure/varray[@name='positions']");
      if (ns_temp.size() > 0) {
         n_temp = it->node().select_nodes(".//structure/varray[@name='positions']")[0].node();
         for (xml_node thisNode = n_temp.child("v"); thisNode; thisNode = thisNode.next_sibling("v")) {
            sscanf(thisNode.child_value(), "  %f  %f  %f ", &x, &y, &z );
            r[0] = 0; r[1] = 0; r[2] = 0;
            r[0] = x*vasprun->latt[0][0];
            r[1] = y*vasprun->latt[1][1];
            r[2] = z*vasprun->latt[2][2];
            thisStep.ppp.push_back(r);
            }
         n_temp = it->node().select_nodes(".//varray[@name='forces']")[0].node();
         for (xml_node thisNode = n_temp.child("v"); thisNode; thisNode = thisNode.next_sibling("v")) {
            sscanf(thisNode.child_value(), "  %f  %f  %f ", &x, &y, &z );
            r[0] = 0; r[1] = 0; r[2] = 0;
            r[0] = x*vasprun->latt[0][0];
            r[1] = y*vasprun->latt[1][1];
            r[2] = z*vasprun->latt[2][2];
            thisStep.fff.push_back(r);
            }
         vasprun->timesteps.push_back(thisStep);
      }
   }



   int ts_count = 0;
   int time_start = config->time_start;
   int time_end = config->time_end;
   vector<TimeStep>  new_timesteps;
   if (time_end > vasprun->timesteps.size()) {
      screen.error << "Specified end time greater than total timesteps.";
      time_end = vasprun->timesteps.size() - 1; 
   }
   if (time_start > vasprun->timesteps.size()) {
      screen.error << "Specified start time (time_start) greater than total timesteps";
      time_start=0;
   }
   if (time_end <= time_start) {
      screen.error << "Specified start time (time_start) is greater than end time (time_end). Disabling this feature";
      time_start=0;
      time_end=vasprun->timesteps.size()-1;
   }


   for (int t=time_start; t<=time_end; t++) {
      new_timesteps.push_back(vasprun->timesteps[t]);
   }
   vasprun->timesteps = new_timesteps;







   cout << "DONE PARSING" << endl;
   
   
   int index_counter=0;
   for (unsigned  i=0; i < vasprun->atoms.size(); i++ ) {
      vasprun->atoms[i].sindex = index_counter;
      vasprun->atoms[i].eindex = index_counter + vasprun->atoms[i].atomspertype-1;
      index_counter = vasprun->atoms[i].eindex+1;
   }

   vasprun->ntimesteps = vasprun->timesteps.size();

   screen.data("SYSTEM", vasprun->system_name);
   screen.data("NATOMS", vasprun->numatoms);
   screen.data("NTYPES", vasprun->numtypes);
   screen.data("NSTEPS", vasprun->ntimesteps);
   screen.data("DT", to_string(vasprun->dt) + " femtosecond");
   screen.data("REALTIME", to_string(vasprun->dt*vasprun->ntimesteps/1000.0) + " picoseconds") ;
   screen.data("LATTICE DIMENSIONS","");
   screen.data("   x",vasprun->latt[0][0]);
   screen.data("   y",vasprun->latt[1][1]);
   screen.data("   z",vasprun->latt[2][2]);


   screen.data("ATOMS","");
   for (int i =0; i<vasprun->numtypes; i++) {
      screen.data("   " + vasprun->atoms[i].element, "");
      screen.data("      COUNT", vasprun->atoms[i].atomspertype);
      screen.data("      MASS", vasprun->atoms[i].mass);
      screen.data("      VALENCE",vasprun->atoms[i].valence );
      screen.data("      PP", vasprun->atoms[i].pseudopotential);
      screen.data("      sINDEX", vasprun->atoms[i].sindex);
      screen.data("      eINDEX", vasprun->atoms[i].eindex);
   }

   
   screen.status << "Splitting up data by atom type";
   vasprun->dataIntoAtoms();

   screen.status << "Unwrapping coordinates";
   vasprun->unwrap();
   vasprun->unwrap_atomType(&vasprun->allatoms);
   

   screen.status << "Weighing the system";
   vasprun->mass_system(); //"weigh" the system.

// Fill the atomindex variable with that atom types' index.
   for (int i=0; i<vasprun->atoms.size(); i++ ) {
      vasprun->atoms[i].atomindex=i;
   }

/* This section writes the unwrapped coordinates to a series of files for animated plotting to ensure proper wrapping */ 
/*
   for (int t=0; t<vasprun->ntimesteps-1; t++ ) {
   ofstream of;
   of.open("diagnostic_scripts/animated_wrapped" + static_cast<ostringstream*>( &(ostringstream() << t) )->str() + ".data");
      of << "\n";
      for (int i=0; i<vasprun->atoms[2].atomspertype-1; i++ ) {
         of << "\t";
         for (int x=0;x<3;x++){
            of << vasprun->atoms[2].timesteps[t].ppp_uw[i][x] << "\t\t";
         }
         of << "\n";
      }
   of.close();
   }
*/


   return 0;
}




/* int main() {
   cout << "\n Starting XML read"<<endl;
   VasprunXML v;
   v.input_filename="vasprun.xml";
   readXML(&v);
   cout << "\n done"<<endl;
}*/





   






















































