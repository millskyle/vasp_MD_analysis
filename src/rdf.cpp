#ifndef RDF_H
#define RDF_H

#include "data_structure.h"
#include <map>

/* int nint(float x) {
   return floor(x+0.5);
} */


GnuPlotScript gnuplot;



int radial_distribution_function(VasprunXML *vasprun, Configuration *config) {

   //find the minimum cell dimension...we can't plot the RDF past half of this.
   double minimum_dimension = 10000000000000000;
   for (int i=0; i<3;i++) {
      if (vasprun->latt[i][i] < minimum_dimension) {
         minimum_dimension=vasprun->latt[i][i];
      }
   }

   

   //For each atom in the requested atom types
//   for (int atomname=0; atomname < config->rdf_atoms.size(); atomname++) {
      //this pointer will point to the atomType object for this type of atom 
      atomType* atomobject0 = &(vasprun->atomfilters[config->rdf_filter_1].atoms.atoms);
      atomType* atomobject1 = &(vasprun->atomfilters[config->rdf_filter_2].atoms.atoms);
      screen.step << "Beginning RDF calculation between "  + config->rdf_filter_1 + " and " + config->rdf_filter_2 ;
      //screen.data("      Number of data points",atomobject0->timesteps[0].ppp.size()*atomobject1->timesteps[0].ppp.size()*atomobject0->timesteps.size());
      
      double rdf_sum=0; //the aggregate sum of of the displacements in the  timestep
      double xdist; // distances that the atom moved in x,y,z
      double ydist;
      double zdist;
      double xa,xb,ya,yb,za,zb;
      double distance;

      int allcount = 0;

      vector<double> bin_cutoff;
      vector<int> bin_count;
      
      int nbins = config->rdf_bins;
      vector<int> bins;
      //fill the bins vector with zeros
      for (int i=0;i<nbins;i++) {
         bins.push_back(0);
         bin_cutoff.push_back(i*minimum_dimension*0.5/nbins);
      }
      double bin_width = minimum_dimension*0.5/nbins;

      double xmin,ymin,zmin;
      double xmax,ymax,zmax;
      double sum_volume = 0;
      //for each timestep which we have positions for
      cout << "         t=0            ";
      for (int t=0; t < atomobject0->timesteps.size(); t++ ) {
         cout << "\r         t="<< t << "            ";
         cout.flush();

         xmin = 1e10;
         ymin = 1e10;
         zmin = 1e10;
         xmax = 0;
         ymax = 0;
         zmax = 0;
         //for each atom in the position vector of vectors
         for (int a=0; a<atomobject0->timesteps[0].ppp.size(); a++) {
            xa = atomobject0->timesteps[t].ppp[a][0];
            ya = atomobject0->timesteps[t].ppp[a][1];
            za = atomobject0->timesteps[t].ppp[a][2];
            //find the min/max coordinates for volume (density) calculation.
            if (xa > xmax) {xmax = xa;}
            if (ya > ymax) {ymax = ya;}
            if (za > zmax) {zmax = za;}
            if (xa < xmin) {xmin = xa;}
            if (ya < ymin) {ymin = ya;}
            if (za < zmin) {zmin = za;}

            for (int b=0; b<atomobject1->timesteps[0].ppp.size(); b++ ) { 
               //don't want to double-count as that wastes resources.
               //But if we don't double count, then we need to multiply by 2 later
               xb = atomobject1->timesteps[t].ppp[b][0];
               yb = atomobject1->timesteps[t].ppp[b][1];
               zb = atomobject1->timesteps[t].ppp[b][2];
               
               xdist = xb-xa;
               xdist = xdist - nint(xdist / vasprun->latt[0][0]) * vasprun->latt[0][0];
               ydist = yb-ya;
               ydist = ydist - nint(ydist / vasprun->latt[1][1]) * vasprun->latt[1][1];
               zdist = zb-za;
               zdist = zdist - nint(zdist / vasprun->latt[2][2]) * vasprun->latt[2][2];


               distance = sqrt(pow(xdist,2) + pow(ydist,2) + pow(zdist,2) );
               if (not(distance < 10e-10)) {
                  allcount++;
                  if (distance < (minimum_dimension) / 2.0) {
                     int binn = floor(distance * nbins / (minimum_dimension*0.5));
                     bins[binn]++;
                  }
               }
            }
         } 
         sum_volume+=(xmax-xmin)*(ymax-ymin)*(zmax-zmin); //volume taken by atoms this timestep
      }
      cout << endl;
  
      //write out the data for this element to an element-specific file
      ofstream of;
      of.open("output/rdf_" + vasprun->label + ".data");     
       
      //write out the gnuplot command, scaling the x-axis increment by the timestep to get it in picoseconds
      gnuplot.command("'rdf_" + vasprun->label + ".data' with lines title 'g(r) between " + config->rdf_filter_1 + " and " + config->rdf_filter_2 + " for " + vasprun->label + "' ls "  + gnuplot.style()  + " lw 3 , "
         , false);

//      double atomic_density = (atomobject0->timesteps[0].ppp.size() + atomobject1->timesteps[0].ppp.size()) / (sum_volume/atomobject0->timesteps.size());
      double atomic_density = allcount / (sum_volume) ;
      //write each timestep to a file
      double normalization;

      for (int n=0; n < nbins; n++) {
         //normalization = atomic_density * (4*3.14159264*bin_cutoff[n]*bin_cutoff[n]*bin_width*atomobject->timesteps.size() );
         normalization = atomic_density * 4.0 * 3.14159 * pow((bin_cutoff[n]),2) * bin_width  * atomobject0->timesteps.size() ; //* (atomobject0->timesteps[0].ppp.size() + atomobject1->timesteps[0].ppp.size() );
         //of << bin_cutoff[n] << "\t" << bins[n]/normalization/(atomobject0->timesteps[0].ppp.size() + atomobject1->timesteps[0].ppp.size() ) << "\n" ;
         of << bin_cutoff[n] << "\t" << bins[n]/normalization << "\n" ;
      }
      of.close();

  // }


   return 0;

}


int radial_distribution_function_wrapper(Configuration *config) {
   if (!config->rdf) {cout << "\nRDF called but not requested in configuration. Exiting"; return 1;}
   screen.status << "Radial Distribution Function";

   //Start a simple bash script which calls GNUplot to plot the msd data
   gnuplot.initialise("rdf","Radial distribution function","r, Angstroms","g(r)","rdf.pdf");
   gnuplot.command("set xrange[0:9]");
   gnuplot.command("set yrange[0:3]");

   gnuplot.command("plot ", false);

   for (int i=0; i<config->vaspruns.size(); i++) {
      radial_distribution_function(&config->vaspruns[i], config);
   }

   cout << "HERE" << endl;

   //Close off the GNUPlot script
   gnuplot.close();

   //add a command to the global plot script to make the msd plots
   config->script_wrapper << "\ngnuplot plot_rdf.gnu \n" ;   

   cout << "HERE" << endl;
   return 0;

}




#endif
