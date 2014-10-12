#ifndef FORCES_CPP
#define FORCES_CPP

#include "data_structure.h"
#include "utility_functions.h"

int force_projections(FileInfo *vasprun, Configuration *config) {
   if (!config->forces) {cout << "\nForce projection called but not requested in configuration. Exiting"; return 1;}
   
   screen.status << "Force Projections";
   screen.step << "Force projections requested from " + config->forces_from_atom + " atoms to " + vec2str(config->msd_atoms) + " atoms ";
   //We need to use unwrapped coordinates.  Unwrap if not already unwrapped.
   vasprun->unwrap(); 
   
   //Make a gnuplot object.  It takes care of writing the data to a script.
   GnuPlotScript gnuplot ;
   gnuplot.initialise("forces","Force","horizontal","vertical","force.pdf");
   gnuplot.command("plot ",false);



   atomType* atomobject0 = vasprun->GetAtom(config->forces_from_atom);

   
   double dx,dy,dz,fx0,fy0,fz0, fx1,fy1,fz1; //components of the vector between atoms
   double proj0, proj1;
   double distance, total_proj;
   vector<vector<double>> all_data;
   vector<double> thisdata;
   double max_distance=0;
   double min_distance=100000000;

   //For each atom in the requested atom types
   for (int atomname=0; atomname < config->forces_to_atoms.size(); atomname++) {
      //this pointer will point to the atomType object for this type of atom 
      atomType* atomobject1 = vasprun->GetAtom(config->forces_to_atoms[atomname]);
      screen.step << "Beginning force calculation for " + atomobject0->element + " to " + atomobject1->element;
      //for each timestep which we have positions for
      for (int t=0; t < atomobject1->timesteps.size()-2; t++ ) {
         //for each atom in the position vector of vectors
         for (int a=0; a<atomobject1->atomspertype-1; a++) {
            //for each atom of the second type
            for (int b=0; b<atomobject0->atomspertype-1; b++) {
               //vector between the two atoms
               dx = atomobject0->timesteps[t].ppp[b][0] - atomobject1->timesteps[t].ppp[a][0];
               dy = atomobject0->timesteps[t].ppp[b][1] - atomobject1->timesteps[t].ppp[a][1];
               dz = atomobject0->timesteps[t].ppp[b][2] - atomobject1->timesteps[t].ppp[a][2];
               //minimum image correction:
               dx = dx - nint(dx / vasprun->latt[0][0])*vasprun->latt[0][0];
               dy = dy - nint(dy / vasprun->latt[1][1])*vasprun->latt[1][1];
               dz = dz - nint(dz / vasprun->latt[2][2])*vasprun->latt[2][2];

               distance = sqrt(dx*dx + dy*dy + dz*dz);
               if (distance > max_distance) { max_distance = distance; }
               if (distance < min_distance) { min_distance = distance; }

               //get the unit vector by dividing by magnitude
               dx = dx/(distance);
               dy = dy/(distance);
               dz = dz/(distance);
               //get the forces
               fx0 = atomobject0->timesteps[t].fff[b][0];
               fy0 = atomobject0->timesteps[t].fff[b][1];
               fz0 = atomobject0->timesteps[t].fff[b][2];

               fx1 = atomobject1->timesteps[t].fff[a][0];
               fy1 = atomobject1->timesteps[t].fff[a][1];
               fz1 = atomobject1->timesteps[t].fff[a][2];
               //project it
               proj0 = fx0*dx + fy0*dy + fz0*dz;
               proj1 = fx1*dx + fy1*dy + fz1*dz;
              
//               if (distance > 6.0) { proj1=0; proj0=0;  }
               
               thisdata.clear();
               thisdata.push_back(distance);
               thisdata.push_back(proj1);
               all_data.push_back(thisdata);
               
               thisdata.clear();
               thisdata.push_back(distance);
               thisdata.push_back(-proj0);
               all_data.push_back(thisdata);


            } 
         }
      }
   
   screen.step << "There are no atoms within less than " + to_string(min_distance) + " angstroms of each other.";

   //make vectors that are the correct size (ie: number of bins)
   vector<int> bins_count;
   vector<double> bins_sum;
   for (int i=0; i<config->forces_bins; i++) {
      bins_count.push_back(0);
      bins_sum.push_back(0.0);
   }
   
   double bin_width = max_distance / config->forces_bins;

   for (int i=0; i<all_data.size(); i++) { 
      bins_count[floor(all_data[i][0] / bin_width)]++;
      bins_sum[floor(all_data[i][0] / bin_width)]+= all_data[i][1];
   }



   
   
   
   //write out the data for this element to an element-specific file
   ofstream of;
   of.open("output/forces_" + atomobject1->element + ".data");     
       
   //write out the gnuplot command, scaling the x-axis increment by the timestep to get it in picoseconds
   gnuplot.command("'forces_" 
      + atomobject1->element 
      + ".data' using 1:2 with lines title '" 
      + atomobject1->element 
      + "' ls "
      + gnuplot.style()
      + " lw 3 , "
      ,false);

   //write each timestep to a file
   for (int i=0; i < bins_sum.size(); i++) {
         of << i*bin_width << "\t" << bins_sum[i]/bins_count[i] << "\n";
        //cout << bins_sum[i] << " / " << bins_count[i] << " = " << bins_sum[i] / bins_count[i] << "\n";
   }
   of.close();

   
}

   gnuplot.close();
   //add a command to the global plot script to make the msd plots
   config->script_wrapper << "\ngnuplot plot_forces.gnu \n" ;   

   return 0;

}

#endif
