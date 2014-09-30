#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include "data_structure.h"
#include "tinyxml/tinyxml.h"
#include "tinyxml/tinystr.h"
//#define TIXML_USE_STL
using namespace std;

typedef TiXmlElement tag;
int str2int(string ss){
   stringstream strValue;
   strValue << ss;
   int intValue;
   strValue >> intValue;
   return intValue;
}

void parse_atomtypes_tag(tag* atomtypesTag, FileInfo *vasprun ){
   for (tag* level1 = atomtypesTag->FirstChildElement(); level1 != NULL; level1=level1->NextSiblingElement()) {
      if (0==strcmp(level1->Value(),"set")) {
         for (tag* rc = level1->FirstChildElement(); rc != NULL; rc = rc->NextSiblingElement()) {
            atomType atomTypeData; //struct to hold the data from this line
            int fieldcounter=0;
            for (tag* c = rc->FirstChildElement(); c!=NULL; c = c->NextSiblingElement()) {
               if (fieldcounter==0) {
                  atomTypeData.atomspertype = str2int(c->FirstChild()->ToText()->Value());
               } else if (fieldcounter==1) {
                  atomTypeData.element = c->FirstChild()->ToText()->Value();
               } else if (fieldcounter==2) {
                  atomTypeData.mass = stod(c->FirstChild()->ToText()->Value());
               } else if (fieldcounter==3) {
                  atomTypeData.valence = stod(c->FirstChild()->ToText()->Value());
               } else if (fieldcounter==4) {
                  atomTypeData.pseudopotential = c->FirstChild()->ToText()->Value();
               }
               fieldcounter++; 
               }
               vasprun->atoms.push_back(atomTypeData);
         }
      }
   }
}

bool update_3d_vector(vector<threevector>* objectToUpdate, float x, float y, float z) {
   threevector r;
   r.push_back(x);r.push_back(y);r.push_back(z);
   (*objectToUpdate).push_back(r); //[0] = x;  (*objectToUpdate)[1] = y; (*objectToUpdate)[2] = z;
   return true;
}

//int readXML(FileInfo& info) {
int readXML(FileInfo *vasprun) {
   int timestep_counter = 0;
   TiXmlDocument doc;
//   if (doc.LoadFile(info.input_filename.c_str())) {
   cout << "Loading file into memory" << endl;
   if (!doc.LoadFile(vasprun->input_filename.c_str())) {
      cerr << doc.ErrorDesc() << endl;
      return 0;
   }
   cout << "File Loaded. Beginning to parse" <<endl;
  

   for (tag* level1 = doc.FirstChildElement(); level1 != NULL; level1 = level1->NextSiblingElement()) {
      if (0==strcmp(level1->Value(),"modeling")) {
         for (tag* level2 = level1->FirstChildElement(); level2 != NULL; level2 = level2->NextSiblingElement()) {
            if (0==strcmp(level2->Value(),"incar")) {
               for (tag* level3 = level2->FirstChildElement(); level3!=NULL; level3 = level3->NextSiblingElement()) {
                  if (0==strcmp(level3->Value(), "i") && 0==strcmp(level3->Attribute("name"),"POTIM")) {
                     vasprun->dt = stod(level3->FirstChild()->ToText()->Value());
                  } else if (0==strcmp(level3->Value(), "i") && 0==strcmp(level3->Attribute("name"),"SYSTEM")) {
                     vasprun->system_name = level3->FirstChild()->ToText()->Value();
                  } else if (0==strcmp(level3->Value(), "i") && 0==strcmp(level3->Attribute("name"),"TEBEG")) {
                     vasprun->starting_temperature = stod(level3->FirstChild()->ToText()->Value());
                  }
               }
            } else if (0==strcmp(level2->Value(),"atominfo")) {
               for (tag* level3 = level2->FirstChildElement(); level3 != NULL; level3 = level3->NextSiblingElement()) {
                  if (0==strcmp(level3->Value(),"atoms")) {
                     vasprun->numatoms = str2int(level3->FirstChild()->ToText()->Value());
                  } else if (0==strcmp(level3->Value(),"types")) {
                     vasprun->numtypes = str2int(level3->FirstChild()->ToText()->Value());
                  } else if (0==strcmp(level3->Value(),"array")) {
                     if (0==strcmp(level3->Attribute("name"), "atomtypes")) {
                        parse_atomtypes_tag(level3, vasprun);
                     }
                  }
               }
            } else if (0==strcmp(level2->Value(),"structure")){
               if (0==strcmp(level2->Attribute("name"),"initialpos")) {
                  for (tag* level3 = level2->FirstChildElement(); level3 != NULL; level3 = level3->NextSiblingElement()) {
                     if (0==strcmp(level3->Value(),"crystal")) {
                        for (tag* level4 = level3->FirstChildElement(); level4 != NULL; level4 = level4->NextSiblingElement()) {
                           if (0==strcmp(level4->Value(),"varray")) {
                              if (0==strcmp(level4->Attribute("name"),"basis")) {
                                 int dimension=1; //are we looking at x, y, or z vector?
                                 for (tag* v = level4->FirstChildElement(); v!=NULL; v = v->NextSiblingElement()) {
                                    float x,y,z;
                                    sscanf(v->FirstChild()->ToText()->Value(), "  %f  %f  %f ", &x, &y, &z );
                                    if (dimension==1){
                                       update_3d_vector(&vasprun->latt, x,y,z);
                                    } else if (dimension == 2) {
                                       update_3d_vector(&vasprun->latt, x,y,z);
                                    } else if (dimension ==3) {
                                       update_3d_vector(&vasprun->latt, x,y,z);
                                    }
                                    dimension++;
                                 }
                              }
                           }
                        }
                     }
                  }              
               }
            } else if (0==strcmp(level2->Value(),"calculation")) {
               TimeStep thisStep;
               for (tag* level3 = level2->FirstChildElement(); level3 != NULL; level3 = level3->NextSiblingElement()) {
                  if (0==strcmp(level3->Value(),"structure")) {
                     for (tag* level4 = level3->FirstChildElement("varray"); level4 != NULL; level4 = level4->NextSiblingElement()) {
                        if (0==strcmp(level4->Attribute("name"),"positions")) {
                           for (tag* vec = level4->FirstChildElement("v"); vec !=NULL; vec = vec->NextSiblingElement()) {
                              float x,y,z;
                              vector<float> r;
                              sscanf(vec->FirstChild()->ToText()->Value(), "  %f  %f  %f ", &x, &y, &z );
                              r.push_back(x*vasprun->latt[0][0]);
                              r.push_back(y*vasprun->latt[1][1]);
                              r.push_back(z*vasprun->latt[2][2]);
                              thisStep.ppp.push_back(r);
                           }
                        } 
                     }  
                  } else if (0==strcmp(level3->Value(),"varray") && 0==strcmp(level3->Attribute("name"),"forces")) {
                     for (tag* vec = level3->FirstChildElement("v"); vec !=NULL; vec = vec->NextSiblingElement()) {
                        float x,y,z;
                        vector<float> r;
                        sscanf(vec->FirstChild()->ToText()->Value(), "  %f  %f  %f ", &x, &y, &z );
                        r.push_back(x);
                        r.push_back(y);
                        r.push_back(z);
                        thisStep.fff.push_back(r);
                     }
                  }
               }
               vasprun->timesteps.push_back(thisStep);
               timestep_counter++;          
            }
         }
      } 
   }

   vasprun->ntimesteps = vasprun->timesteps.size();
   cout << "SYSTEM: " << vasprun->system_name << endl;
   cout << "Number of atoms: " << vasprun->numatoms << endl;
   cout << "Number of atom types: " << vasprun->numtypes << endl;
   for (int i =0; i<vasprun->numtypes; i++) {
      cout << "\tThere are " << vasprun->atoms[i].atomspertype << " " << vasprun->atoms[i].element << " atoms each with mass " 
           << vasprun->atoms[i].mass << " and " << vasprun->atoms[i].valence 
           << " valence electrons.\n\t\t--> The " << vasprun->atoms[i].pseudopotential << " pseudopotential was used."  << endl;
   }

   cout << "The x lattice vector is <" << vasprun->latt[0][0] << ", " << vasprun->latt[0][1] << ", " << vasprun->latt[0][2] << ">." <<endl;
   cout << "The y lattice vector is <" << vasprun->latt[1][0] << ", " << vasprun->latt[1][1] << ", " << vasprun->latt[1][2] << ">." <<endl;
   cout << "The z lattice vector is <" << vasprun->latt[2][0] << ", " << vasprun->latt[2][1] << ", " << vasprun->latt[2][2] << ">." <<endl;

   //atomType* theatomIwant;
   int index_counter=0;
   for (unsigned  i=0; i < vasprun->atoms.size(); i++ ) {
      vasprun->atoms[i].sindex = index_counter;
      vasprun->atoms[i].eindex = index_counter + vasprun->atoms[i].atomspertype-1;
      cout << "   Atom "<<i<<" will go from index "<<vasprun->atoms[i].sindex<<" to " << vasprun->atoms[i].eindex << "." <<endl;
      index_counter = vasprun->atoms[i].eindex+1;
   }

   int junk=vasprun->dataIntoAtoms();
   atomType* theatomIwant = vasprun->GetAtom("Al");
   cout << "I chose the " << theatomIwant->element << " atoms which starts at index " << theatomIwant->sindex << endl ;
   cout << "  For the " << theatomIwant->element << " atoms, there are " << theatomIwant->timesteps.size()  << " timesteps with timestep 0 having " << theatomIwant->timesteps[0].ppp.size() << " ppp elements." << endl;



   cout << "There are " << vasprun->ntimesteps << " timesteps in the file with dt=" << vasprun->dt << " for a total time of " << vasprun->dt*vasprun->ntimesteps/1000.0 << " picoseconds." << endl;
   cout << "For timestep 0 (the first timestep):" <<endl;
   cout << "\tThere are "<<vasprun->timesteps[0].ppp.size() << " position vectors." << endl;
   cout << "\tThere are "<<vasprun->timesteps[0].fff.size() << " force vectors." << endl;
   cout << "\tThe first x coordinate of the position vector of the first atom in this timestep is " << vasprun->timesteps[0].ppp[0][0] << "\n";


vasprun->unwrap();
vasprun->mass_system(); //"weigh" the system.

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
   FileInfo v;
   v.input_filename="vasprun.xml";
   readXML(&v);
   cout << "\n done"<<endl;
}*/





   






















































