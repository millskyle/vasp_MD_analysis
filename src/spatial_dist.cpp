#ifndef SPATIAL_DIST_H
#define SPATIAL_DIST_H

#include "data_structure.h"

int spatial_distribution_projection(FileInfo *vasprun, Configuration *config) {
   if (!config->spatial_distribution_projection) {cout << "\nSpatial Distribution called but not requested in configuration. Exiting"; return 1;}

   cout << "---Spatial Distribution ---" <<endl;
   
   //Start a simple bash script which calls GNUplot to plot the msd data
   GnuPlotScript gnuplot;
   gnuplot.initialise("sdp","Spatial Distribution","x, Angstroms","y, Angstroms","sdp.pdf", "Reds");


   int nbins_x = config->nbins_x;
   int nbins_y = config->nbins_y;

  
  /*
   //make empty array of correct size to hold counts
   vector <vector <int>> bins;
   vector<int> dummy_row;
   for (int i=0; i<nbins_x; i++){
      dummy_row.push_back(0);
   }
   for (int j=0; j< nbins_y; j++) {
      bins.push_back(dummy_row);
   }
*/


   unsigned bins[200][200] = {};





   //For each atom in the requested atom types
   for (int atomname=0; atomname < config->liquid_atoms.size(); atomname++) {
      //this pointer will point to the atomType object for this type of atom 
      atomType* atomobject = vasprun->GetAtom(config->msd_atoms[atomname]);
      cout << "Starting to bin the " << atomobject->element << " atoms for spatial distribution.\n";   
      
   
   
   //fill an array with all positions from all timesteps to speed up the iteration:
   vector<threevector> pppp;
   cout << "Filling array with all timesteps\n";
   for (int t=0; t<atomobject->timesteps.size(); t++) {
      for (int a=0; a<atomobject->timesteps[t].ppp.size(); a++) {
         pppp.push_back(atomobject->timesteps[t].ppp[a]);
   
      }
   }

   


      cout << "Binning" << "\n";
      double xbin_min, xbin_max, ybin_min, ybin_max = 0;

      int dim1 = 0; //dimension 1 (ie: 0 for x)
      int dim2 = 1; //dimension 2 (ie: 1 for x)

      //for each timestep which we have positions for
//      for (int t=0; t < atomobject->timesteps.size(); t++ ) {
            //for each bin in the x direction
            for (int ix=0; ix < nbins_x; ix++){
               //for each bin in the y direction
               for (int iy=0; iy < nbins_y; iy++) {
                  xbin_min = vasprun->latt[dim1][dim1] * ix / nbins_x;
                  xbin_max = vasprun->latt[dim1][dim1] * (ix+1) / nbins_x; 
                  ybin_min = vasprun->latt[dim2][dim2] * iy / nbins_y;
                  ybin_max = vasprun->latt[dim2][dim2] * (iy+1) / nbins_y;
                  cout << "bin (" << ix << "," << iy << ")\n";               
                  bins[ix][iy] = bins[ix][iy] 
                     + count_atoms_in_box_pointer(&pppp, 
                        xbin_min, xbin_max, 
                        ybin_min, ybin_max, 
                        -9999999,9999999);
               
               }
            } 
  //       }
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
   gnuplot.command("plot 'sdp.data' matrix with image");


   //Close off the GNUPlot script
   gnuplot.close();

   //add a command to the global plot script to make the msd plots
   config->script_wrapper << "\ngnuplot plot_sdp.gnu \n" ;   

   return 0;

}

#endif
