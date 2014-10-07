#ifndef SPATIAL_DIST_H
#define SPATIAL_DIST_H

#include "data_structure.h"

int spatial_distribution_projection(FileInfo *vasprun, Configuration *config) {
   if (!config->spatial_distribution_projection) {cout << "\nSpatial Distribution called but not requested in configuration. Exiting"; return 1;}

   cout << "---Spatial Distribution ---" <<endl;
   
   //Start a simple bash script which calls GNUplot to plot the msd data
   GnuPlotScript gnuplot;
   gnuplot.initialise("sdp","Spatial Distribution","x [Angstrom]","y [Angstrom]","sdp.pdf", "YlOrRd");



   int nbins_x = config->nbins_x;
   int nbins_y = config->nbins_y;
   int ix;
   int iy;
   int dim1 = 0;
   int dim2 = 1;
   
   //set the dimension references based on the configuration file
   if (config->collapse_dimension==0){
      dim1 = 1;
      dim2 = 2;
   } else if (config->collapse_dimension==1){
      dim1 = 0;
      dim2 = 2;
   } else {
      dim1 = 0;
      dim2 = 1;
   }
    
 
   cout << "Plotting spatial distribution heat map on axis " << dim1 << " vs. " << dim2 << ".\n";

  
   gnuplot.command("set xrange [0:" + to_string(vasprun->latt[dim1][dim1]) + "] ");
   gnuplot.command("set yrange [0:" + to_string(vasprun->latt[dim2][dim2]) + "] ");
   
   gnuplot.command("set pm3d map");
   gnuplot.command("set pm3d interpolate 0,0");
   gnuplot.command("set size ratio -1");
   gnuplot.command("set key off");



   //make empty vector of vectors of correct size to hold counts
   vector <vector <int>> bins;
   vector<int> dummy_row;
   for (int i=0; i<=nbins_x; i++){
      dummy_row.push_back(0);
   }
   for (int j=0; j<=nbins_y; j++) {
      bins.push_back(dummy_row);
   }




   //For each atom in the requested atom types
   for (int atomname=0; atomname < config->liquid_atoms.size(); atomname++) {
      //this pointer will point to the atomType object for this type of atom 
      atomType* atomobject = vasprun->GetAtom(config->msd_atoms[atomname]);
      cout << "Starting to bin the " << atomobject->element << " atoms for spatial distribution.\n";   
      
      //for each timestep
      for (int t=0; t < atomobject->timesteps.size()-1; t++ ) {
         //for each atom in the timestep
         for (int a=0; a<atomobject->timesteps[t].ppp.size(); a++) {
            //get the bin index in which the atom falls
            //   ie: take the x coordinate and divide by the bin width, round to nearest integer.
            ix = nint(atomobject->timesteps[t].ppp[a][dim1] / (vasprun->latt[dim1][dim1] / nbins_x));
            iy = nint(atomobject->timesteps[t].ppp[a][dim2] / (vasprun->latt[dim2][dim2] / nbins_y));
            bins[ix][iy]++;

         }


       }
   } 
 
   ofstream of;
   of.open("output/sdp.data");

   for (int ix=0; ix < nbins_x-1; ix++){
   //for each bin in the y direction
      for (int iy=0; iy < nbins_y-1; iy++) {
         of << bins[ix][iy] << "\t" ;
      }
      of << "\n";
   }
   of.close();

      //write out the gnuplot command, scaling the x-axis increment by the timestep to get it in picoseconds
   gnuplot.command("splot 'sdp.data' using ($1 * " 
         + to_string(vasprun->latt[dim1][dim1] / nbins_x)
         + "):($2 * "
         + to_string(vasprun->latt[dim2][dim2] / nbins_y)
         + "):3  matrix "         
         );


   //Close off the GNUPlot script
   gnuplot.close();

   //add a command to the global plot script to make the msd plots
   config->script_wrapper << "\ngnuplot plot_sdp.gnu \n" ;

   return 0;

}

#endif
