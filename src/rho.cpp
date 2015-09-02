#ifndef RHO_H
#define RHO_H

#include "data_structure.h"
#include "utility_functions.h"
#include <climits>

Conversions converter;

int global_density(VasprunXML *vasprun, Configuration *config) {
   if (!config->rho) {cout << "\nDensity function called but not requested in configuration. Exiting"; return 1;}

   screen.status << "Density calculations" ;
  
   //Make a gnuplot object.  It takes care of writing the data to a script. 
   GnuPlotScript gnuplot ;
   gnuplot.initialise("rho","Density vs. time","Time, [picoseconds]","Density","rho.pdf");
   gnuplot.command("set yrange [0:]");
   gnuplot.command("plot ",false);

   //For each atom in the requested atom types
   //for (int atomname=0; atomname < config->rho_atoms.size(); atomname++) {
      //this pointer will point to the atomType object for this type of atom 
      atomType* atomobject = &(config->atomfilters[config->rho_filter].atoms.atoms);
      
      screen.step << "Beginning Density calculation for filter " + config->rho_filter;   
      double x0,y0,z0;
      double x1,y1,z1;
      int atomcount;
      //for each timestep which we have positions for
      for (int t=0; t < atomobject->timesteps.size(); t++ ) {
         //for each atom in the position vector of vectors

         vector<threevector>* positions= &atomobject->timesteps[t].ppp;
         
         x0 = 10000000.0;
         y0 = 10000000.0;
         z0 = 10000000.0;
         x1 = -10000000.0;
         y1 = -10000000.0;
         z1 = -10000000.0;
         for (int a=0; a<positions->size(); a++) {
            if ((*positions)[a][0] < x0) { x0 = (*positions)[a][0]; }
            if ((*positions)[a][0] > x1) { x1 = (*positions)[a][0]; }
            if ((*positions)[a][1] < y0) { y0 = (*positions)[a][1]; }
            if ((*positions)[a][1] > y1) { y1 = (*positions)[a][1]; }
            if ((*positions)[a][2] < z0) { z0 = (*positions)[a][2]; }
            if ((*positions)[a][2] > z1) { z1 = (*positions)[a][2]; }         
         }
         
//         int atomcount = count_atoms_in_box(*positions, x0, x1, y0, y1, z0, z1);
         atomcount = atomobject->timesteps[0].ppp.size(); 
         atomobject->timesteps[t].density = atomcount / ((x1-x0)*(y1-y0)*(z1-z0));
      }
      screen.step << x1-x0;
      screen.step << y1-y0;
      screen.step << z1-z0;

   //write out the data for this element to an element-specific file
   ofstream of;
   of.open("output/" + config->rho_data_prefix + atomobject->element + ".data");

   //write out the gnuplot command, scaling the x-axis increment by the timestep to get it in picoseconds
   gnuplot.command("'" 
      + config->rho_data_prefix 
      + atomobject->element 
      + ".data' using ($0*" 
      + to_string(vasprun->dt) 
      + "*0.001):1 with lines title '" 
      + atomobject->element 
      + " (effective)' ls " 
      + gnuplot.style() 
      + " lw 3 , "
      ,false);

   gnuplot.command("'" 
      + config->rho_data_prefix 
      + atomobject->element 
      + ".data' using ($0*" 
      + to_string(vasprun->dt) 
      + "*0.001):2 with lines title '" 
      + atomobject->element 
      + " (real)' ls " 
      + gnuplot.style() 
      + " lw 3, "
      , false);





   //write each timestep to a file
   for (int t=0; t < atomobject->timesteps.size(); t++) {
      if (atomobject->timesteps[t].density>0) {
         of << (  atomobject->timesteps[t].density/atomobject->timesteps[0].ppp.size() * pow(converter.angstroms_per_cm,3)*
                  atomobject->mass *
                  converter.moles_per_gram )  << "\t" <<  config->rho_experimental << " \n" ;

         //cout << atomobject->timesteps[t].density << "\n";
      }
   }
   of.close();   
//}



   //Close off the GNUPlot
   gnuplot.close();

   //add a command to the global plot script to make the msd plots
   config->script_wrapper << "\ngnuplot plot_rho.gnu \n" ;   

   return 0;

}

#endif
