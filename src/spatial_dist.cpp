#ifndef SPATIAL_DIST_H
#define SPATIAL_DIST_H

#include "data_structure.h"
#include "utility_functions.h"

int spatial_distribution_projection(VasprunXML *vasprun, Configuration *config) {
   if (!config->spatial_distribution_projection) {cout << "\nSpatial Distribution called but not requested in configuration. Exiting"; return 1;}

   screen.status <<  "Spatial Distribution projected onto plane";
   
   //Start a simple bash script which calls GNUplot to plot the msd data
   GnuPlotScript gnuplot;
   gnuplot.initialise("sdp","Spatial Distribution","x [Angstrom]","y [Angstrom]","sdp.pdf", "YlOrRd");
   gnuplot.term(" png font \"Times,0.1\" size 3840,2160 ","sdp.png");

   int nbins_h = config->sd_nbins_h;
   int nbins_v = config->sd_nbins_v;
   int ix;
   int iy;
   int dim1 = 0;
   int dim2 = 1;
   
   //set the dimension references based on the configuration file
   if (config->sd_collapse_dimension==0){
      dim1 = 1;
      dim2 = 2;
   } else if (config->sd_collapse_dimension==1){
      dim1 = 0;
      dim2 = 2;
   } else {
      dim1 = 0;
      dim2 = 1;
   }
    
 
 screen.step << "Plotting spatial distribution heat map on axis " + to_string(dim1) + " vs. " + to_string(dim2);

  
   gnuplot.command("set xrange [0:" + to_string(vasprun->latt[dim1][dim1]) + "] ");
   gnuplot.command("set yrange [0:" + to_string(vasprun->latt[dim2][dim2]) + "] ");
   
   gnuplot.command("set pm3d map");
   gnuplot.command("set pm3d interpolate 0,0");
   gnuplot.command("set size ratio -1");
   gnuplot.command("set key off");



   //make empty vector of vectors of correct size to hold counts
   vector <vector <int>> bins;
   vector<int> dummy_row;
   for (int i=0; i<=nbins_h; i++){
      dummy_row.push_back(0);
   }
   for (int j=0; j<=nbins_v; j++) {
      bins.push_back(dummy_row);
   }

     




     atomType* staticatoms = &(vasprun->atomfilters[config->sd_static_filter].atoms.atoms);
     double static_atom_coords_1[staticatoms->timesteps[0].ppp.size()];
     double static_atom_coords_2[staticatoms->timesteps[0].ppp.size()];


//     for (int t=0; t < staticatoms->timesteps.size(); t++) {
        int t=0;
        for (int a=0; a < staticatoms->timesteps[0].ppp.size(); a++) {
           static_atom_coords_1[a]+=staticatoms->timesteps[t].ppp[a][dim1];
           static_atom_coords_2[a]+=staticatoms->timesteps[t].ppp[a][dim2];
        }
//     }
     //divide by number of timesteps to get average
//     for (int a=0; a < staticatoms->timesteps[0].ppp.size(); a++) {
//        static_atom_coords_1[a]=static_atom_coords_1[a] / staticatoms->timesteps.size();
//        static_atom_coords_2[a]=static_atom_coords_2[a] / staticatoms->timesteps.size();
//     }
    
   //For each atom in the requested atom types
//   for (int atomname=0; atomname < config->liquid_atoms.size(); atomname++) {
      //this pointer will point to the atomType object for this type of atom 
      atomType* dynamicatoms = &(vasprun->atomfilters[config->sd_dynamic_filter].atoms.atoms); 



      screen.step << "Binning atoms in the filter " + config->sd_dynamic_filter;   
      
      //for each timestep
      for (int t=0; t < dynamicatoms->timesteps.size(); t++ ) {
         //get lattice centre of mass:
         //for each atom in the timestepi
         for (int a=0; a<dynamicatoms->timesteps[t].ppp.size(); a++) {
//         cout << " a = " << a << endl;
            for (int xx=0; xx<config->sd_h_wrap_times; xx++) {
              
               for (int yy=0; yy<config->sd_v_wrap_times; yy++) {
                  ix = nint((dynamicatoms->timesteps[t].ppp[a][dim1] + xx*config->sd_h_wrap_dist) / (vasprun->latt[dim1][dim1] / nbins_h));
                  iy = nint((dynamicatoms->timesteps[t].ppp[a][dim2] + yy*config->sd_v_wrap_dist) / (vasprun->latt[dim2][dim2] / nbins_v));
                  while (ix >= nbins_h) {
                     ix = ix - nbins_h;
                  }
                  while (iy >= nbins_v) {
                     iy = iy - nbins_v;
                  }
               bins[ix][iy]++;
                  
                  
                  ix = nint((dynamicatoms->timesteps[t].ppp[a][dim1] + (xx+0.5)*config->sd_h_wrap_dist) / (vasprun->latt[dim1][dim1] / nbins_h));
                  iy = nint((dynamicatoms->timesteps[t].ppp[a][dim2] + (yy + 1.0 + sqrt(2) )*config->sd_v_wrap_dist) / (vasprun->latt[dim2][dim2] / nbins_v));
                  while (ix >= nbins_h) {
                     ix = ix - nbins_h;
                  }
                  while (iy >= nbins_v) {
                     iy = iy - nbins_v;
                  }
               bins[ix][iy]++;
               
               }
            }
            //get the bin index in which the atom falls
            //   ie: take the x coordinate and divide by the bin width, round to nearest integer.
            //ix = nint(dynamicatoms->timesteps[t].ppp[a][dim1] / (vasprun->latt[dim1][dim1] / nbins_h));
//            iy = nint(dynamicatoms->timesteps[t].ppp[a][dim2] / (vasprun->latt[dim2][dim2] / nbins_v));
//            bins[ix][iy]++;
         }
       }
   //} 

   
   double norm;
   norm = (vasprun->latt[dim1][dim1] / nbins_h) 
      * (vasprun->latt[dim2][dim2] / nbins_v) 
      * vasprun->latt[config->sd_collapse_dimension][config->sd_collapse_dimension] 
      * config->sd_h_wrap_times 
      * 2 * config->sd_v_wrap_times 
      * dynamicatoms->timesteps.size() 
      * vasprun->dt;
   ofstream of;
   of.open("output/sdp_dyn.data");

   double imax=-1000000;
   double imin=10000000;
   double k = 0;
   //for each bin in the y direction
   for (int iy=0; iy < nbins_v; iy++) {
      for (int ix=0; ix < nbins_h; ix++){

//            if (  nint(floor( (iy * vasprun->latt[dim1][dim1] / nbins_h)) ) % 2 == 0 ) {
//               of << "100" << "\t";
//            } else {
//               of << "0" << "\t";
//            }
            
           k = pow(bins[ix][iy]/norm,2);
           imax = max(imax, k);
           imin = min(imin, k);
           of << k << "\t" ;
//         }

      }
      of << "\n";
   }
   of.close();

//   imin = 0;
//   imax = 100;

   gnuplot.command("set cbrange [" + to_string(imax/2.0) + ":" + to_string(imax) + "]  " );

   of.open("output/sdp_sta.data");
   int halftime = staticatoms->timesteps.size()/2;
   t=halftime;

//   for (int t=0; t<staticatoms->timesteps.size(); t+=1000) {
      for (int i=0; i<staticatoms->timesteps[t].ppp.size(); i++) {
         of << staticatoms->timesteps[t].ppp_uw[i][dim1] 
            << "\t" 
            << staticatoms->timesteps[t].ppp_uw[i][dim2]   
            << "\t" 
            <<  atomic_symbol_to_number(staticatoms->symbols[i])/2. 
            << "\n";
      }
//   }

   of.close();

//   gnuplot.command("set palette model RGB defined ( 5 '#FFC4EB', 12 '#222222', 13 '#CECECE', 7 '#0028FF'  )");

      //write out the gnuplot command, scaling the x-axis increment by the timestep to get it in picoseconds
   gnuplot.command("splot 'sdp_dyn.data' using ($1 * " 
         + to_string(vasprun->latt[dim1][dim1] / nbins_h)
         + "):($2 * "
         + to_string(vasprun->latt[dim2][dim2] / nbins_v)
         + "):3  matrix, 'sdp_sta.data' using 1:2:(0):3 w points pt 6 lt 1 ps variable  "         
         );


   //Close off the GNUPlot script
   gnuplot.close();

   //add a command to the global plot script to make the msd plots
   config->script_wrapper << "\ngnuplot plot_sdp.gnu \n" ;

   return 0;



}

#endif
