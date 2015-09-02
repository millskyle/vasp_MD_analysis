#ifndef STORAGE_H_
#define STORAGE_H_

#include <string>
#include <map>
#include <algorithm>
#include <vector>
#include <string.h>
#include <cmath>
#include <regex>
#include "utility_functions.h"
#include "screen.h"

using namespace std;




/* Code I copied/modified to display some information about maps */
typedef int DEBUG_LEVEL;
int d_level = 0;

template<typename Key, typename Value>
std::ostream& operator<<(std::ostream& os, const std::pair<const Key, Value>& p)
{
       os << p.first << " => " << &p.second;
           return os;
}

template<typename Container>
void printmap(DEBUG_LEVEL level, const Container& c) {
    if (level >= d_level) {
        for(typename Container::const_iterator it = c.begin();
            it != c.end(); ++it)
        std::cout << *it << '\n';
    }
} 
/* Code I copied/modified to display some information about maps */



typedef vector<float> threevector;

struct key {
   string name;
   string value;
};



struct INCAR {
   vector<key> keys;
};

struct TimeStep {
   vector<threevector> ppp;
   vector<threevector> ppp_uw; //unwrapped
   vector<threevector> fff;
   double MSD=-999;
   vector<double> COM;
   double density=0;
   TimeStep() {
      COM.push_back(0);
      COM.push_back(0);
      COM.push_back(0);
   }
};


struct atomType {
   int atomindex; //the index of the atom referenced by the parent.
   int atomspertype, sindex, eindex;
   string element,pseudopotential;
   double mass, valence;
   vector<string> symbols;
   vector<TimeStep> timesteps; //holds all time dependent data
   bool COM_already=false;
   atomType () {
      atomspertype=0; element="X  ";mass = 0.00; valence = 0.00; pseudopotential="garbage";
   }
}; 


struct atomSet {
   vector<string> symbols;
   atomType atoms;
};


struct VasprunXML {
   //File input/output parameters
   string input_filename,output_data_location,system_name;
   INCAR incar;
   //File data structures
   map <string, double> atomic_mass;
   int numatoms,numtypes,ntimesteps;
   bool unwrapped_already=false;
   vector<int> atom_count;
//   vector<double> latt_x,latt_y,latt_z;
   vector<threevector> latt;
   vector<threevector> COM; //centre of mass
   vector<atomType> atoms;
   vector<TimeStep> timesteps;
   vector<string> ion_symbols;
   float totalMass=0;
   double dt,starting_temperature; //timestep length, delta t

   atomType allatoms;
       
  
   int dataIntoAtoms(){
      int counter=0;
      for (unsigned i=0; i<atoms.size(); i++) {
         vector<TimeStep> alltimes;
         for (unsigned t=0;t<ntimesteps-1; t++) {
            TimeStep ts;
            for (unsigned a=atoms[i].sindex;a<=atoms[i].eindex; a++) {
               ts.ppp.push_back(timesteps[t].ppp[a] ) ;
               ts.fff.push_back(timesteps[t].fff[a] ) ;
            }
            alltimes.push_back(ts);
         }
         atoms[i].timesteps = alltimes;
      }
      allatoms.timesteps = timesteps;
      return 1;
   }



  
   atomType GetAtomByIndex(int sID, int eID) {
      atomType filtered;
      vector<TimeStep> alltimes;
      for (unsigned t=0;t<ntimesteps-1; t++) {
         TimeStep ts;
         for (unsigned i=sID; i<=eID; i++)  {
            ts.ppp.push_back(  allatoms.timesteps[t].ppp[i]  );
            ts.fff.push_back(  allatoms.timesteps[t].fff[i]  );
         }
         alltimes.push_back(ts);
      }
      filtered.timesteps = alltimes;
      return filtered;
   }
      


   

   atomType* GetAtom(string element) { 
   //returns a pointer to a specific atoms object, based on the passed element symbol
   //ie:   vasprun->GetAtom("Al")   returns a pointer the object holding the Aluminum atom info
      for (unsigned i=0; i<atoms.size(); i++) {
         if (atoms[i].element==element) {
            return &atoms[i];
         }
      }
      //If the function hasn't returned yet, then the atom type was not found. Kill the program since
      //we'll later get a Segmentation Fault for pointing to a non-existent object.
      cout << endl << endl << "FATAL ERROR!!!" <<endl<< "ATOM TYPE "<<element<<" NOT FOUND. EXITING." << endl;
      exit(0);
   }


   int unwrap_atomType(atomType *atoms) {
      //Copy the position vectors
      int sign;
      for (unsigned t=0; t < atoms->timesteps.size(); t++) {
         atoms->timesteps[t].ppp_uw = atoms->timesteps[t].ppp;
      }

      for (unsigned t=1; t < atoms->timesteps.size(); t++) {//for each timestep
         vector<threevector> &x0 = atoms->timesteps[t-1].ppp_uw;
         vector<threevector> &x1 = atoms->timesteps[t].ppp_uw;
         for (unsigned a=0; a<x0.size(); a++) { //for atom in ppp vector
            for (int x=0; x<3; x++) {  //for each dimension (0=x,1=y,2=z) {
               if ((abs(x0[a][x] -x1[a][x])) > latt[x][x]/2 ) { //if the difference is greater than half lv
                  if (x0[a][x] < x1[a][x]) {sign=-1;} //if 
                  else {sign=1;}
                  x1[a][x] = x1[a][x] + sign*latt[x][x];
               }
            }
         }
      }
      return 0;
   }


   int unwrap() {
      int sign;
      if (unwrapped_already) {return 0;}

      for (unsigned i=0; i<atoms.size(); i++) {  // for each atom type
         unwrap_atomType(&atoms[i]);
      }

      unwrapped_already=true;
      return 0;
   }


   int mass_system(){
      //for each atom type, add the product of the mass and the number of atoms to the totalMass variable.
      if (totalMass==0) {  //Only do this if totalMass is 0 (meaning it hasn't already been calculated).
         for (unsigned i=0; i<atoms.size(); i++) {
            totalMass+=atoms[i].mass * atoms[i].atomspertype;
         }
         return 0;
      } else {
         return 1;}
   }




   VasprunXML() {
      numatoms=0;
      input_filename = "/tmp/garbage";
      output_data_location = "/tmp/";
 
      /*latt_x.push_back(0);
      latt_x.push_back(0);
      latt_x.push_back(0);
      latt_y.push_back(0);
      latt_y.push_back(0);
      latt_y.push_back(0);
      latt_z.push_back(0);
eparation distance [Angstrom]
6
      latt_z.push_back(0);*/

   }
};


int calculate_COM(atomType *atoms, VasprunXML *vasprun) {
   if (atoms->COM_already) {return 0;}
   

   double atomTotalMass = 0;
   double COMx,COMy,COMz = 0;


   for (unsigned a=0; a<atoms->timesteps[0].ppp.size(); a++) {
      atomTotalMass+=vasprun->atomic_mass[atoms->symbols[a]];
   }



   for (unsigned t=0; t < atoms->timesteps.size(); t++) { //for each timestep
      COMx = 0;
      COMy = 0;
      COMz = 0;
         vector<threevector> &r = atoms->timesteps[t].ppp_uw;
         for (unsigned a=0; a<r.size(); a++) { //for atom in ppp vector
            COMx +=  (vasprun->atomic_mass[atoms->symbols[a]] * r[a][0]);
            COMy +=  (vasprun->atomic_mass[atoms->symbols[a]] * r[a][1]);
            COMz +=  (vasprun->atomic_mass[atoms->symbols[a]] * r[a][2]);
            
         }
      atoms->timesteps[t].COM[0] = COMx / atomTotalMass;
      atoms->timesteps[t].COM[1] = COMy / atomTotalMass;
      atoms->timesteps[t].COM[2] = COMz / atomTotalMass;

   }
   atoms->COM_already=true;

//print it out to check:


return 0;


}

struct atomfilter {
   string name;
   string filter_type;
   string criteria;
   atomSet atoms;
   vector<int> filter_indices;
   

   //////////////
   // Funtion to execute the filter, filling the atoms object specified a couple of lines above with only
   // the atoms that are included in the filter.  This function called from main.cpp after file parsing done.
   //////////////
   int execute_filter(VasprunXML* vasprun) {
      if (filter_type == "index_range") {
         //split the string on commas to get the ranges in separate elements
         vector<string> ranges = str2vec(criteria,",");
         //for each range
         for (int i=0; i < ranges.size(); i++) {
            //split each range now on a colon, denoting the
            vector<string> thisRange = str2vec(ranges[i], ":");
            //get the starting/ending indices
            int sind = stoi(thisRange[0]);
            int eind = stoi(thisRange[1]);
            //make sure that the starting index is smaller than (or equal to) the ending index
            if (sind > eind) {
               screen.error << "INVALID FILTER INDICES. Start index must be less than or equal to end index";
            } else {
               for (int jj=sind; jj<=eind; jj++) {
                  filter_indices.push_back(jj);
               }
            }
         }
      
      } else if ( filter_type == "symbol" ) {
         //split the criteria string up on spaces or commas
         vector<string> desiredSymbols = str2vec(criteria, ", ");
         //for each symbol in the specified filter symbols
         for (int s=0; s<desiredSymbols.size(); s++) {
            //for each ion in the datafile
            for (int i = 0; i < vasprun->ion_symbols.size(); i++ ) {
               //if the ion is of the desired type, add its index to the list
               if (vasprun->ion_symbols[i]==desiredSymbols[s]) {
                  filter_indices.push_back(i);
               }
            }
         }
      } else if ( filter_type == "index" || filter_type == "indices" ) {
         //split the criteria string up on spaces or commas
         vector<string> desiredIndices = str2vec(criteria, ", ");
         for (int i=0; i < desiredIndices.size(); i++) {
            filter_indices.push_back( stoi(rtrim(ltrim(desiredIndices[i]))));
         }
      }
      //sort the vector by index.
      sort(filter_indices.begin(), filter_indices.end() );
   
   
      //Now we need to take the vector of indices to include and actually filter the timestep data
    
      ///////
      //Fill a vector with the symbols of the atoms contained in the filter
      double mass = 0;
      int j=0;
      vector<string> filtered_symbols;
      for (int i=0; i<filter_indices.size(); i++ ) {
         j = filter_indices[i];
         filtered_symbols.push_back(vasprun->ion_symbols[j]);
         mass += vasprun->atomic_mass[vasprun->ion_symbols[j]];

      }
      ///////


      ///////
      //Fill the atomType object with time-dependent data
      atomType filtered_atoms;
      vector<TimeStep> alltimes;

      for (unsigned t=0; t<vasprun->allatoms.timesteps.size(); t++) {
         TimeStep ts;
         for (int i=0; i<filter_indices.size(); i++ ) {
//            cout << t << ":  " <<  filter_indices[i] << endl;
            ts.ppp.push_back( vasprun->allatoms.timesteps[t].ppp[filter_indices[i]]) ;
            ts.fff.push_back( vasprun->allatoms.timesteps[t].fff[filter_indices[i]]) ;
            ts.ppp_uw.push_back( vasprun->allatoms.timesteps[t].ppp_uw[filter_indices[i]]) ;
         }
         alltimes.push_back(ts);
      }
      filtered_atoms.timesteps = alltimes;
      filtered_atoms.symbols = filtered_symbols;




      //update the filter.
      atoms.symbols = filtered_symbols;
      atoms.atoms = filtered_atoms;
      
      //weigh the atoms in the filter
      atoms.atoms.mass=mass;
//      for (int i=0; i<filtered_atoms.timesteps[0].ppp.size(); i++) {
//         atoms.atoms.mass += vasprun->atomic_mass[filtered_symbols[i]];
//      }

      cout << "MASS: " << atoms.atoms.mass << endl;

     
      return 0; 
   }

};



struct Configuration {
   string tempstr;
   
   bool msd;
   string msd_data_prefix;
   string msd_filter;
   double msd_reference_D;

   bool rho;
   string rho_data_prefix;
   string rho_filter;
   double rho_experimental;


   bool spatial_distribution_projection;
   bool spatial_distribution_lattice;
   int collapse_dimension;
   vector<string> lattice_atoms;
   vector<string> liquid_atoms;
   int nbins_x;
   int nbins_y;

   bool rdf;
   string rdf_filter;
   string rdf_data_prefix;
   int rdf_bins;
   bool rdf_cut_half_lv;
   string rdf_plot_type;   
   

   bool forces;
   int forces_bins;
   string forces_set_1;
   string forces_set_2;

   vector<string> filter_name_list;

   bool index_show;

   map <string, atomfilter> atomfilters;   

   ofstream script_wrapper;
   string script_wrapper_location = "output/make_all_plots.sh";

   ofstream log;
   string log_file_location = "log";
  
      


   
} config;
   

























struct GnuPlotScript {
   string name;
   string title;
   string xlabel;
   string ylabel;
   string cmap;
   string output;
   int linestyle;
   ofstream script; 
   
   string style() {
      linestyle++;
      if (linestyle >=8) {
         linestyle=1;
      }
      return to_string(linestyle);
   }
   
   void command(string str,bool newline=true) {
      script << str;
      if (newline) {
         script << "\n";
      } 
   }
   void initialise(string iname,string ititle,string ixlabel,string iylabel,string ioutput, string icmap="Set2") {
      name = iname;
      title = ititle;
      xlabel = ixlabel;
      ylabel = iylabel;
      cmap = icmap;
      output = ioutput;
      linestyle = 1;
 
      script.open("output/plot_" + name + ".gnu");
      script << "#!/usr/bin/gnuplot\n\nreset\n"
              << "set title \"" << title << "\"\n"
              << "set term pdf font \"Times,8\"\n"
              << "set output \"" << output << "\"\n"
              << "set xlabel \"" << xlabel << "\"\n"
              << "set ylabel \"" << ylabel << "\"\n"
              << "load '../supporting/plotting/gnuplot_colormaps/" + cmap + ".plt' \n";
   }
   void close() {
      script.close();
   }
};

#endif
