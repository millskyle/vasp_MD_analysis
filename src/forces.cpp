#ifndef FORCES_CPP
#define FORCES_CPP

#include "data_structure.h"
#include "utility_functions.h"

int force_field(FileInfo *vasprun, Configuration *config) {
   if (!config->forces) {cout << "\nForce projection called but not requested in configuration. Exiting"; return 1;}   
   screen.status << "Force-field calculations";
   screen.step << "force field requested for " + vec2str(config->force_field_atoms) + " atoms";
   //We need to use unwrapped coordinates.  Unwrap if not already unwrapped.
   
   //Make a gnuplot object.  It takes care of writing the data to a script.
   GnuPlotScript gnuplot ;
   gnuplot.initialise("forcefield","Force","horizontal","vertical","forcefield.pdf");
   gnuplot.command("set yrange [0:" + to_string(vasprun->latt[2][2]) + "]");
   gnuplot.command("set xrange [0:" + to_string(vasprun->latt[1][1]) + "]");
   gnuplot.command("plot ",false);

   double dx,dy,dz,fx,fy,fz; //components of the vectors between atoms
   vector<vector<double>> all_data;
   vector<double> thisdata;
   double max_distance=0;
   double min_distance=100000000;
   
   vector<double> bin_width;
   bin_width.push_back(vasprun->latt[1][1] / config->force_field_resolution);
   bin_width.push_back(vasprun->latt[2][2] / config->force_field_resolution);
   
   vector<vector<vector<double>>> alldata;

   vector<double> onebin;
   onebin.push_back(0.0);
   onebin.push_back(0.0);
   onebin.push_back(0.0);
   onebin.push_back(0.0);
   onebin.push_back(0.0);

   vector<vector<double>> onerow;

   //make a "row vector" by pushing n bins into a vector
   for (int i=0; i<config->force_field_resolution; i++) {
      onerow.push_back(onebin);
   }
   //push N rows into the alldata vector
   for (int j=0; j<config->force_field_resolution; j++) {
      alldata.push_back(onerow);
   }

   //fill the first two elements of every bin with the y and z coordinates of the centre of that bin
   for (int i=0; i<config->force_field_resolution; i++) {
      for (int j=0; j< config->force_field_resolution; j++) {
         alldata[i][j][0] = i*bin_width[0] + 0.5*bin_width[0];
         alldata[i][j][1] = j*bin_width[1] + 0.5*bin_width[1];
      }
   }


   int thisbiny;
   int thisbinz;

   //For each atom in the requested atom types
   for (int atomname=0; atomname < config->force_field_atoms.size(); atomname++) {
      //this pointer will point to the atomType object for this type of atom 
      atomType* atomobject = vasprun->GetAtom(config->force_field_atoms[atomname]);
      screen.step << "Beginning force field calculation for " + atomobject->element; 
      //for each valid timestep
      for (int t=0; t < atomobject->timesteps.size()-2; t++ ) {
         //for each atom in the position vector of vectors
         for (int a=0; a<atomobject->atomspertype-1; a++) {
            //for each atom of the second type

                  //get the forces from the atom objects
                  fy = atomobject->timesteps[t].fff[a][1];
                  fz = atomobject->timesteps[t].fff[a][2];
                  
                  thisbiny = floor(atomobject->timesteps[t].ppp[a][1] / bin_width[0]);
                  thisbinz = floor(atomobject->timesteps[t].ppp[a][2] / bin_width[1]);

                  alldata[thisbiny][thisbinz][4]++;
                  alldata[thisbiny][thisbinz][2]+=fy;
                  alldata[thisbiny][thisbinz][3]+=fz;

         } 
      }
   


   
   //write out the gnuplot command
   gnuplot.command(
      "'force_field_" 
      + atomobject->element 
      + ".data' using 1:2:3:4 with vectors head filled title '" 
      + atomobject->element 
      + "' "
      ,false);

   
   double norm;
   //write out the data for this element to an element-specific file
   ofstream of;
   of.open("output/force_field_" + atomobject->element + ".data");     
   for (int i=0; i<config->force_field_resolution; i++) {
      for (int j=0; j< config->force_field_resolution; j++) {
         //norm=  sqrt(x**2 + y**2) * count
         norm = sqrt(alldata[i][j][2]*alldata[i][j][2] + alldata[i][j][3]*alldata[i][j][3])  ;
         if (norm==0) {
            norm=0.0000001;
         }
         of << alldata[i][j][0] 
            << "\t" 
            << alldata[i][j][1] 
            << "\t" 
            << bin_width[0]*alldata[i][j][2]/(4*norm)
            << "\t"
            << bin_width[0]*alldata[i][j][3]/(4*norm)
            <<"\n";
      }
   }
//for (int i=0; i < bins_sum.size(); i++) {
//      of << i*bin_width << "\t" << bins_sum[i]/bins_count[i] << "\t" << bins_std_dev[i] << "\n";
//   }
   of.close();

   
}

   gnuplot.close();
   //add a command to the global plot script to make the msd plots
   config->script_wrapper << "\ngnuplot plot_forcefield.gnu \n" ;   

   return 0;

}






























































int force_bond_projections(FileInfo *vasprun, Configuration *config) {
   if (!config->forces) {cout << "\nForce projection called but not requested in configuration. Exiting"; return 1;}   
   screen.status << "Force Projections";
   screen.step << "Bond force projections requested from " + config->forces_from_atom + " atoms to " + vec2str(config->forces_to_atoms) + " atoms ";
   //We need to use unwrapped coordinates.  Unwrap if not already unwrapped.
   vasprun->unwrap(); 
   
   //Make a gnuplot object.  It takes care of writing the data to a script.
   GnuPlotScript gnuplot ;
   gnuplot.initialise("forces","Force","distance","Average force","force.pdf");
//   gnuplot.command("set yrange [-1:1]");
   gnuplot.command("set xrange [0:]");
   gnuplot.command("plot ",false);

   //first atom_type object
   atomType* atomobject0 = vasprun->GetAtom(config->forces_from_atom);
   
   double dx,dy,dz,fx0,fy0,fz0,fx1,fy1,fz1; //components of the vectors between atoms
   double proj0, proj1; //scalar projections
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
      //for each valid timestep
      for (int t=0; t < atomobject1->timesteps.size()-2; t++ ) {
//         cout << "t=" << t << "\n";
         //for each atom in the position vector of vectors
         for (int a=0; a<atomobject1->atomspertype; a++) {
            //for each atom of the second type
            for (int b=0; b<atomobject0->atomspertype; b++) {

               //vector between the two atoms
               //if (rand()%10==1 || 1==1) {  //only actually calculate 10% of the force projections. (for speed).
                  dx = atomobject0->timesteps[t].ppp[b][0] - atomobject1->timesteps[t].ppp[a][0];
                  dy = atomobject0->timesteps[t].ppp[b][1] - atomobject1->timesteps[t].ppp[a][1];
                  dz = atomobject0->timesteps[t].ppp[b][2] - atomobject1->timesteps[t].ppp[a][2];
                  //minimum image correction:
                  dx = dx - nint(dx / vasprun->latt[0][0])*vasprun->latt[0][0];
                  dy = dy - nint(dy / vasprun->latt[1][1])*vasprun->latt[1][1];
                  dz = dz - nint(dz / vasprun->latt[2][2])*vasprun->latt[2][2];
                  //distance between both atoms:
                  distance = sqrt(dx*dx + dy*dy + dz*dz);
                  if (distance > max_distance) { max_distance = distance; }
                  if (distance < min_distance) { min_distance = distance; }
                  //get the forces from the atom objects
                  fx0 = atomobject0->timesteps[t].fff[b][0];
                  fy0 = atomobject0->timesteps[t].fff[b][1];
                  fz0 = atomobject0->timesteps[t].fff[b][2];
                  fx1 = atomobject1->timesteps[t].fff[a][0];
                  fy1 = atomobject1->timesteps[t].fff[a][1];
                  fz1 = atomobject1->timesteps[t].fff[a][2];
                 
         //         cout << "forces:  atom1: " << fx0 << "," << fy0 << "," << fz0 << ".\n" ; 
           //       cout << "forces:  atom2: " << fx1 << "," << fy1 << "," << fz1 << ".\n" ; 

                  //project them
                  // (F dot r) / |r|
                  proj0 = (fx0*dx + fy0*dy + fz0*dz)/distance;
                  proj1 = (fx1*dx + fy1*dy + fz1*dz)/distance;

         //         cout << "Projection: atom1 onto r: " << proj0 << "\n";
         //         cout << "Projection: atom2 onto r: " << proj1 << "\n";
                   
                  //push the first projection into the all_data vector
                  thisdata.clear();
                  thisdata.push_back(distance);
                  //one projection has to be negative because dx,dy,dz vector is in wrong direction for it. 
                  thisdata.push_back(-proj1);
                  all_data.push_back(thisdata);
                  
                  //push the second projection
                  thisdata.clear();
                  thisdata.push_back(distance);
                  thisdata.push_back(proj0);
                  all_data.push_back(thisdata);
//               } //endif
            } 
         }
      }
   
   //make vectors that are the correct size (ie: number of bins)
   //to hold the data
   vector<int> bins_count;
   vector<double> bins_sum;
   vector<double> bins_std_dev;
   for (int i=0; i<config->forces_bins; i++) {
      bins_std_dev.push_back(0);
      bins_count.push_back(0);
      bins_sum.push_back(0.0);
   }
  
   //calculate the width of each bin 
   double bin_width = max_distance / config->forces_bins;

   //for each data point, put it in the correct bin and increment the bin counter
   for (int i=0; i<all_data.size(); i++) { 
      //cout << "distance=" << all_data[i][0] << " and force=" << all_data[i][1] << "\n";
      bins_count[floor(all_data[i][0] / bin_width)]++;
      bins_sum[floor(all_data[i][0] / bin_width)]+= all_data[i][1];
   }

   //calculate the standard deviation for each bin
   int thisbin;
   for (int i=0; i<all_data.size(); i++) {
      thisbin = floor(all_data[i][0] / bin_width);
      bins_std_dev[thisbin] += pow((all_data[i][1] - (bins_sum[thisbin] / bins_count[thisbin])),2);
   }

   //normalize/sqrt the standard deviation
   for (int i=0; i<bins_std_dev.size(); i++) {
      bins_std_dev[i]=bins_std_dev[i] / bins_count[i];
      bins_std_dev[i] = sqrt(bins_std_dev[i]);
   }


   
   //write out the gnuplot command
   gnuplot.command(
      "'forces_" 
      + atomobject1->element 
      + ".data' using 1:2 with lines title '" 
      + atomobject1->element 
      + "' ls "
      + gnuplot.style()
      + " lw 3 , "
      + " 'forces_" 
      + atomobject1->element 
      + ".data' using 1:($2-$3) with lines ls " 
      + gnuplot.style() 
      + " lw 0.5, "
      + " 'forces_" 
      + atomobject1->element 
      + ".data' using 1:($2+$3) with lines ls " 
      + gnuplot.style() 
      + " lw 0.5, "
      ,false);

   float smallest_dimension;
   smallest_dimension = vasprun->latt[0][0];
   if (vasprun->latt[1][1] < smallest_dimension) {
      smallest_dimension = vasprun->latt[1][1];
   }
   if (vasprun->latt[2][2] < smallest_dimension) {
      smallest_dimension = vasprun->latt[2][2];
   }
   
   //write out the data for this element to an element-specific file
   ofstream of;
   of.open("output/forces_" + atomobject1->element + ".data");     
   for (int i=0; i <= bins_sum.size(); i++) {
      if (bins_count[i]>0) {
         if (i*bin_width < smallest_dimension/2.0) {
            of << i*bin_width << "\t" << bins_sum[i]/bins_count[i] << "\t" << bins_std_dev[i] << "\n";
         }
      }
   }
   of.close();

   
}

   gnuplot.close();
   //add a command to the global plot script to make the msd plots
   config->script_wrapper << "\ngnuplot plot_forces.gnu \n" ;   

   return 0;

}

#endif
