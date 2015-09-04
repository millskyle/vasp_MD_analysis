#ifndef MSD_H
#define MSD_H

#include "data_structure.h"
#include "utility_functions.h"

Conversions convert;



int mean_square_displacement(VasprunXML *vasprun, Configuration *config, GnuPlotScript *gnuplot) {
   


   //For each atom in the requested atom types
      //this pointer will point to the atomType object for this type of atom 
      atomType* atomobject = &(vasprun->atomfilters[config->msd_filter].atoms.atoms);
   

      //Calculate the center of mass for this atom_type. It'll be stored in the object.
      screen.step << "Calculating center of mass for " + atomobject->element;
      int garbageint = calculate_COM(atomobject,vasprun);
      screen.step << "Beginning MSD calculation for " + atomobject->element;
      int msd_count=0; //the integer number of data points that go into the aggregate sum
      double msd_sum=0; //the aggregate sum of of the displacements in the  timestep
      double xdist; // distances that the atom moved in x,y,z
      double ydist;
      double zdist;
      //for each timestep which we have positions for
      for (int t=1; t < atomobject->timesteps.size()-2; t++ ) {
         //for each atom in the position vector of vectors
         for (int a=0; a<atomobject->timesteps[0].ppp.size(); a++) {
               int t2=0; //the timestep to calculate displacement from
               // xdist = r(t) - r(0) - [  COM(t)  - COM(t) ]
               xdist = atomobject->timesteps[t].ppp_uw[a][0] - atomobject->timesteps[t2].ppp_uw[a][0]
                       - atomobject->timesteps[t].COM[0] + atomobject->timesteps[t2].COM[0];
               ydist = atomobject->timesteps[t].ppp_uw[a][1] - atomobject->timesteps[t2].ppp_uw[a][1]
                       - atomobject->timesteps[t].COM[1] + atomobject->timesteps[t2].COM[1];
               zdist = atomobject->timesteps[t].ppp_uw[a][2] - atomobject->timesteps[t2].ppp_uw[a][2]
                       - atomobject->timesteps[t].COM[2] + atomobject->timesteps[t2].COM[2];
               //the sum is the modulus of this vector w/o square root.
               msd_sum+=pow(xdist,2) + pow(ydist,2) + pow(zdist,2);
               //increment the sample counter (we'll use this to average later).
               msd_count++;
         }
         //push the data in to the object
         atomobject->timesteps[t].MSD= ( msd_sum / msd_count );
         //reset the counters
         msd_sum=0;
         msd_count=0;
      }

   //write out the data for this element to an element-specific file
   ofstream of;
   of.open("output/msd_" + vasprun->label + ".data");     
       
   //write out the gnuplot command
   gnuplot->command("'msd_" 
      + vasprun->label 
      + ".data' using 1:2 with lines title '" 
      + vasprun->label 
      + "' ls "
      + gnuplot->style()
      + " lw 3 , "
//      + " '"
//      + "msd_" + vasprun->label 
//      + ".data' using 1:3"
//      + " with lines title '"
//      + "Theoretical"
//      + "' ls 16 "
//      + " lw 3 ,"
      ,false);

   //write each timestep to a file
   for (int t=0; t < atomobject->timesteps.size(); t++) {
      if (atomobject->timesteps[t].MSD>=0) {
         of << vasprun->dt/1000.0 * t 
            << "\t" 
            << atomobject->timesteps[t].MSD 
            << "\t" 
            << 6*stod(config->msd_reference_D) * (t*vasprun->dt*convert.femto2pico)   * ( convert.m2A * convert.m2A / convert.to_pico   )
            << "\n" ;

      }
   }
   of.close();

   return 0;

}


int msd_wrapper(Configuration *config) {
   if (!config->msd) {cout << "\nMSD called but not requested in configuration. Exiting"; return 1;}
   screen.status << "Mean Squared Displacement";
  
   GnuPlotScript gnuplot;
   //Make a gnuplot object.  It takes care of writing the data to a script.
   gnuplot.initialise("msd","Mean Square displacement","Time, [picoseconds]","Distance squared [angstroms^2]","msd.pdf");
   gnuplot.command("plot ",false);

   for (int i=0; i<config->vaspruns.size(); i++ ) {
      mean_square_displacement(&config->vaspruns[i], config, &gnuplot);
   }
   
   //write command to plot the experimental mean square displacement
   gnuplot.command(" 6*x*" + to_string((convert.m2A*convert.m2A/convert.to_pico)) + "*" + config->msd_reference_D + " w l ls " + gnuplot.style() + " title 'Experimental (D=" + config->msd_reference_D + ")' lw 3 " );

   gnuplot.close();
   //add a command to the global plot script to make the msd plots
   config->script_wrapper << "\ngnuplot plot_msd.gnu \n" ;   

   return 0;

}



#endif
