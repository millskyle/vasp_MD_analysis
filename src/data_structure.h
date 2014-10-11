#ifndef STORAGE_H_
#define STORAGE_H_

#include <string>
#include <vector>
#include <string.h>
#include <cmath>

using namespace std;

typedef vector<float> threevector;

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
   vector<TimeStep> timesteps; //holds all time dependent data
   bool COM_already=false;
   atomType () {
      atomspertype=0; element="X  ";mass = 0.00; valence = 0.00; pseudopotential="garbage";
   }
}; 


struct FileInfo {
   //File input/output parameters
   string input_filename,output_data_location,system_name;
   //File data structures
   int numatoms,numtypes,ntimesteps;
   bool unwrapped_already=false;
   vector<int> atom_count;
//   vector<double> latt_x,latt_y,latt_z;
   vector<threevector> latt;
   vector<threevector> COM; //centre of mass
   vector<atomType> atoms;
   vector<TimeStep> timesteps;
   float totalMass=0;
   double dt,starting_temperature; //timestep length, delta t

   int dataIntoAtoms(){
      for (unsigned i=0; i<atoms.size(); i++) {
         vector<TimeStep> alltimes;
         for (unsigned t=0;t<ntimesteps-1; t++) {
            TimeStep ts;
            for (unsigned a=atoms[i].sindex;a<atoms[i].eindex; a++) {
               ts.ppp.push_back(timesteps[t].ppp[a] ) ;
               ts.fff.push_back(timesteps[t].fff[a] ) ;
            }
            alltimes.push_back(ts);
         }
         atoms[i].timesteps = alltimes;
      }
      return 1;
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

   int unwrap() {
      int sign;
      if (unwrapped_already) {return 0;}

      //copy the position vectors
      for (unsigned i=0; i<atoms.size(); i++) {  // for each atom type
         for (unsigned t=0; t < ntimesteps-1; t++) { //for each timestep
            atoms[i].timesteps[t].ppp_uw = atoms[i].timesteps[t].ppp;
         }
      }
      for (unsigned i=0; i<atoms.size(); i++) {  // for each atom type
         for (unsigned t=1; t < ntimesteps-1; t++) { //for each timestep
//            cout << "t=" <<  t << "\n";
//            atoms[i].timesteps[t].ppp_uw = atoms[i].timesteps[t].ppp;
            vector<threevector> &x0 = atoms[i].timesteps[t-1].ppp_uw;
            vector<threevector> &x1 = atoms[i].timesteps[t].ppp_uw;
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


   int calculate_COM(atomType *atoms) {
      if (atoms->COM_already) {return 0;}
      
      double atomTotalMass = atoms->mass * atoms->atomspertype;
      double COMx,COMy,COMz = 0;
      for (unsigned t=0; t < ntimesteps-1; t++) { //for each timestep
         COMx = 0;
         COMy = 0;
         COMz = 0;
            vector<threevector> &r = atoms->timesteps[t].ppp_uw;
            for (unsigned a=0; a<r.size(); a++) { //for atom in ppp vector
               COMx +=  (atoms->mass * r[a][0]);
               COMy +=  (atoms->mass * r[a][1]);
               COMz +=  (atoms->mass * r[a][2]);
               
            }
         atoms->timesteps[t].COM[0] = COMx / atomTotalMass;
         atoms->timesteps[t].COM[1] = COMy / atomTotalMass;
         atoms->timesteps[t].COM[2] = COMz / atomTotalMass;

      }
      atoms->COM_already=true;

//print it out to check:


   return 0;


   }
   




   FileInfo() {
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
      latt_z.push_back(0);
      latt_z.push_back(0);*/

   }
};

struct Configuration {
   string tempstr;
   
   bool msd;
   string msd_data_prefix;
   vector<string> msd_atoms;

   bool rho;
   string rho_data_prefix;
   vector<string> rho_atoms;

   bool spatial_distribution_projection;
   bool spatial_distribution_lattice;
   int collapse_dimension;
   vector<string> lattice_atoms;
   vector<string> liquid_atoms;
   int nbins_x;
   int nbins_y;

   bool rdf;
   string rdf_data_prefix;
   vector<string> rdf_atoms;
   int rdf_bins;
   bool rdf_cut_half_lv;
   string rdf_plot_type;   
   
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
