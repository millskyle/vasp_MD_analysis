#ifndef MSD_H
#define MSD_H

#include "data_structure.h"
#include "utility_functions.h"

Conversions convert;

int mean_square_displacement(FileInfo *vasprun, Configuration *config) {
   if (!config->msd) {cout << "\nMSD called but not requested in configuration. Exiting"; return 1;}
   
   screen.status << "Mean Squared Displacement";
   screen.step << "Requested for " + to_string(config->msd_atoms.size())+ " atom types: " + vec2str(config->msd_atoms);
   //We need to use unwrapped coordinates.  Unwrap if not already unwrapped.
   vasprun->unwrap(); 
   
   //Make a gnuplot object.  It takes care of writing the data to a script.
   GnuPlotScript gnuplot ;
   gnuplot.initialise("msd","Mean Square displacement","Time, [picoseconds]","Disance squared, angstroms^2","msd.pdf");
   gnuplot.command("plot ",false);


   //For each atom in the requested atom types
   for (int atomname=0; atomname < config->msd_atoms.size(); atomname++) {
      //this pointer will point to the atomType object for this type of atom 
      atomType* atomobject = vasprun->GetAtom(config->msd_atoms[atomname]);
      //Calculate the center of mass for this atom_type. It'll be stored in the object.
      screen.step << "Calculating center of mass for " + atomobject->element;
      vasprun->calculate_COM(atomobject);
      screen.step << "Beginning MSD calculation for " + atomobject->element;
      int msd_count=0; //the integer number of data points that go into the aggregate sum
      double msd_sum=0; //the aggregate sum of of the displacements in the  timestep
      double xdist; // distances that the atom moved in x,y,z
      double ydist;
      double zdist;
      //for each timestep which we have positions for
      for (int t=1; t < atomobject->timesteps.size()-2; t++ ) {
         //for each atom in the position vector of vectors
         for (int a=0; a<atomobject->atomspertype-1; a++) {
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
   of.open("output/" + config->msd_data_prefix + atomobject->element + ".data");     
       
   //write out the gnuplot command, scaling the x-axis increment by the timestep to get it in picoseconds
   gnuplot.command("'" 
      + config->msd_data_prefix 
      + atomobject->element 
      + ".data' using ($0*" 
      + to_string(vasprun->dt)
      + "*0.001):1 with lines title '" 
      + atomobject->element 
      + "' ls "
      + gnuplot.style()
      + " lw 3 , "
      + " '"
      + config->msd_data_prefix
      + atomobject->element
      + ".data' using ($0*"
      + to_string(vasprun->dt)
      + "*0.001):2 with lines title '"
      + "Theoretical"
      + "' ls "
      + gnuplot.style()
      + " lw 3 ,"
      ,false);

   double experiment_diffusion_coefficient = 7.8e-9; //in m^2/s 
   //write each timestep to a file
   for (int t=0; t < atomobject->timesteps.size(); t++) {
      if (atomobject->timesteps[t].MSD>=0) {
         of << atomobject->timesteps[t].MSD 
            << "\t" 
            << 6*experiment_diffusion_coefficient * (t*vasprun->dt*convert.femto2pico)   * ( convert.m2A * convert.m2A / convert.to_pico   )
            << "\n" ;
//         cout << atomobject->timesteps[t].MSD << "\n" ;
      }
   }
   of.close();

   
}

   gnuplot.close();
   //add a command to the global plot script to make the msd plots
   config->script_wrapper << "\ngnuplot plot_msd.gnu \n" ;   

   return 0;

}

#endif
