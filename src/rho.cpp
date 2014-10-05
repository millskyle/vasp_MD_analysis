#ifndef RHO_H
#define RHO_H

#include "data_structure.h"
#include "utility_functions.h"

Conversions converter;

int global_density(FileInfo *vasprun, Configuration *config) {
   if (!config->rho) {cout << "\nDensity function called but not requested in configuration. Exiting"; return 1;}

   cout << "--- Starting Density calculations ---" <<endl;
   cout << "Density plots requested for " << config->rho_atoms.size() << " atom types: " << vec2str(config->rho_atoms) << endl;
   //We need to use unwrapped coordinates.  Unwrap if not already unwrapped.
   //vasprun->unwrap(); 
   
   //Start a simple bash script which calls GNUplot to plot the msd data
   ofstream of2;
   of2.open("output/plot_rho.sh");
   of2 << "#!/bin/bash" << "\n" 
       << "gnuplot -persist << GNUPLOTINPUT" << "\n"
       << "set title \"Density vs. Time\"\n"
       << "set term pdf\n"
       << "set output \"rho.pdf\"\n"
       << "set xlabel \"Time, picoseconds\"\n"
       << "set yrange [0:]\n"
       << "set ylabel \"Density, atoms per cubic Angstrom\"\n"
//       << "unset key" << "\n"
       << "plot ";

   //For each atom in the requested atom types
   for (int atomname=0; atomname < config->rho_atoms.size(); atomname++) {
      //this pointer will point to the atomType object for this type of atom 
      atomType* atomobject = vasprun->GetAtom(config->rho_atoms[atomname]);
      
      cout << "Beginning Density calculation for " << atomobject->element << ".\n";   
      //for each timestep which we have positions for
      for (int t=1; t < atomobject->timesteps.size()-2; t++ ) {
         //for each atom in the position vector of vectors

         vector<threevector>* positions= &atomobject->timesteps[t].ppp;
         
         double x0,y0,z0 = 1000000000;
         double x1,y1,z1 = -1000000000;
         for (int a=0; a<positions->size(); a++) {
            if ((*positions)[a][0] < x0) { x0 = (*positions)[a][0]; }
            if ((*positions)[a][0] > x1) { x1 = (*positions)[a][0]; }
            if ((*positions)[a][1] < y0) { y0 = (*positions)[a][1]; }
            if ((*positions)[a][1] > y1) { y1 = (*positions)[a][1]; }
            if ((*positions)[a][2] < z0) { z0 = (*positions)[a][2]; }
            if ((*positions)[a][2] > z1) { z1 = (*positions)[a][2]; }         
         }

         int atomcount = count_atoms_in_box(*positions, x0, x1, y0, y1, z0, z1);
         atomobject->timesteps[t].density = atomcount / ((x1-x0)*(y1-y0)*(z1-z0));   
      }

   //write out the data for this element to an element-specific file
   ofstream of;
   of.open("output/" + config->rho_data_prefix + atomobject->element + ".data");

   //write out the gnuplot command, scaling the x-axis increment by the timestep to get it in picoseconds
   of2 << "'" << config->rho_data_prefix + atomobject->element << ".data' using (\\$0*" << vasprun->dt << "*0.001):1 with lines title '" << atomobject->element << " (effective)' , ";
   of2 << "'" << config->rho_data_prefix + atomobject->element << ".data' using (\\$0*" << vasprun->dt << "*0.001):2 with lines title '" << atomobject->element << " (real)' , ";

   //write each timestep to a file
   for (int t=0; t < atomobject->timesteps.size(); t++) {
      if (atomobject->timesteps[t].MSD>=0) {
         of << (  atomobject->timesteps[t].density * pow(converter.angstroms_per_cm,3)*
                  atomobject->mass *
                  converter.moles_per_gram )  << 
               "\t" << 
               (  atomobject->atomspertype  / ((vasprun->latt[0][0]-1) * (vasprun->latt[1][1]-1) * (vasprun->latt[2][2]-1))
                  * pow(converter.angstroms_per_cm,3)
                  * atomobject->mass
                  * converter.moles_per_gram ) <<
               "\n" ;

         //cout << atomobject->timesteps[t].density << "\n";
      }
   }
   of.close();   
}



   //Close off the GNUPlot bash script
   of2 << "\nGNUPLOTINPUT\n";
   of2.close();

   //add a command to the global plot script to make the msd plots
   config->script_wrapper << "\nbash plot_rho.sh \n" ;   

   return 0;

}

#endif
