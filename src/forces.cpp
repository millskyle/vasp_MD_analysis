#ifndef FORCES_CPP
#define FORCES_CPP

#include "data_structure.h"
#include "utility_functions.h"




int force_bond_projections(VasprunXML *vasprun, Configuration *config, GnuPlotScript* gnuplot) {
   if (!config->forces) {cout << "\nForce projection called but not requested in configuration. Exiting"; return 1;}   

   cout << "HERE" << endl;

   atomType *atomobject0;
   atomType *atomobject1;

   atomobject0 = &(config->atomfilters[config->forces_set_1.c_str()].atoms.atoms);
   atomobject1 = &(config->atomfilters[config->forces_set_2.c_str()].atoms.atoms);
   

   double dx,dy,dz,fx0,fy0,fz0,fx1,fy1,fz1; //components of the vectors between atoms
   double proj0, proj1; //scalar projections
   double distance, total_proj;
   vector<vector<double>> all_data;
   vector<double> thisdata;
   double max_distance=0;
 
 //2015-09-11
  //make vectors that are the correct size (ie: number of bins)
   //to hold the data
   vector<int> bins_count;
   vector<double> bins_sum ;
   vector<double> bins_std_dev;
   for (int i=0; i<config->forces_bins; i++) {
      bins_std_dev.push_back(0);
      bins_count.push_back(0);
      bins_sum.push_back(0.0);
   }



   //calculate the width of each bin 
   double bin_width = max_distance / config->forces_bins;

   double min_distance=100000000;

   cout << "HERE" << endl;

   cout << atomobject1->timesteps[0].ppp.size() << endl;
   cout << atomobject0->timesteps[0].ppp.size() << endl;
   cout << (atomobject1->timesteps.size()-2) << endl ;

   int number_of_projections = 2*atomobject1->timesteps[0].ppp.size()*atomobject0->timesteps[0].ppp.size()*(atomobject1->timesteps.size()-2);
   cout << "HERE" << endl;



//   double all_data_array[number_of_projections][2];

   cout << "array declared" << endl;

   screen.step << "Beginning force projection."; 

   screen.data("Projections"," 2 orientations x " 
         + to_string(atomobject0->timesteps[0].ppp.size())  
         + " atoms x "  
         + to_string(atomobject1->timesteps[0].ppp.size()) 
         +  " atoms x " 
         + to_string((atomobject1->timesteps.size()-2)) 
         +  " timesteps = "
         + to_string(number_of_projections)
         + " projections"
   );
      //for each valid timestep
   int counter=0;
   int thisBin;
   cout << "        t = 0         0.00%     " ;
   for (int t=0; t < atomobject1->timesteps.size()-2; t++ ) {   // -2 as the last timestep could be incomplete
      cout << "\r        t = " << t << "         " << t*100.0/(atomobject1->timesteps.size()-2) << "%       " ;
      cout.flush();
      for (int a=0; a<atomobject1->timesteps[0].ppp.size(); a++) {
         for (int b=0; b<atomobject0->timesteps[0].ppp.size(); b++) {
            //vector between the two atoms
            cout << "1" << endl;
            dx = atomobject0->timesteps[t].ppp[b][0] - atomobject1->timesteps[t].ppp[a][0];
            dy = atomobject0->timesteps[t].ppp[b][1] - atomobject1->timesteps[t].ppp[a][1];
            dz = atomobject0->timesteps[t].ppp[b][2] - atomobject1->timesteps[t].ppp[a][2];
            //minimum image correction:
            dx = dx - nint(dx / vasprun->latt[0][0])*vasprun->latt[0][0];
            dy = dy - nint(dy / vasprun->latt[1][1])*vasprun->latt[1][1];
            dz = dz - nint(dz / vasprun->latt[2][2])*vasprun->latt[2][2];
            //distance between both atoms:
            cout << "2" << endl;
            distance = sqrt(dx*dx + dy*dy + dz*dz);
            max_distance = max(distance, max_distance);
            min_distance = min(distance, min_distance);

            //if (distance > max_distance) { max_distance = distance; }
            //if (distance < min_distance) { min_distance = distance; }

            //get the forces from the atom objects
            fx0 = atomobject0->timesteps[t].fff[b][0];
            fy0 = atomobject0->timesteps[t].fff[b][1];
            fz0 = atomobject0->timesteps[t].fff[b][2];
            fx1 = atomobject1->timesteps[t].fff[a][0];
            fy1 = atomobject1->timesteps[t].fff[a][1];
            fz1 = atomobject1->timesteps[t].fff[a][2];
           
            cout << "3" << endl;
            //project them
            // (F dot r) / |r|
            proj0 = (fx0*dx + fy0*dy + fz0*dz)/distance;
            proj1 = (fx1*dx + fy1*dy + fz1*dz)/distance;

//            all_data_array[counter][0] = distance;
//            all_data_array[counter][1] = -proj1;
            counter++;   //Increment the counter. The reverse projection needs to go into another element
//            all_data_array[counter][0] = distance;
//            all_data_array[counter][1] = proj0;
            counter ++; 
      //for each data point, put it in the correct bin and increment the bin counter
   //for (int i=0; i<number_of_projections-2; i++) { 
            cout << "4" << endl;
            thisBin = floor(distance / bin_width);
            bins_count[thisBin]+=2;
            bins_sum[thisBin]+= -proj1;
            bins_sum[thisBin]+= proj0;

   //}

         


         } 
      }
   }
   cout << "\r";

   screen.data("Separation distances","");
   screen.data("   Minimum",min_distance);
   screen.data("   Maximum",max_distance);

   
   //gnuplot.command("set xrange [" + to_string(min_distance) + ":]"); 
//   gnuplot.command("set xrange [0:]"); // + to_string(min_distance) + ":]"); 

   
//   //calculate the standard deviation for each bin
//   int thisbin;
//   for (int i=0; i<number_of_projections-2; i++) {
//      thisbin = floor(all_data_array[i][0] / bin_width);
//      bins_std_dev[thisbin] += pow((all_data_array[i][1] - bins_sum[thisbin]/bins_count[thisbin]),2);
//   }


//   //normalize/sqrt the standard deviation
//   for (int i=0; i<bins_std_dev.size(); i++) {
////      bins_std_dev[i]=bins_std_dev[i] / (bins_count[i],2);
//      bins_std_dev[i] = sqrt(bins_std_dev[i]) / bins_count[i];
//   }

   //write out the gnuplot command
   gnuplot->command(  
      + " 'forces_" + vasprun->label + ".data' using 1:($2-$3):($2+$3) with filledcurves lc rgb '#D1D1D1' lw 0 title '" + vasprun->label + "' ,  '' using 1:2 w l ls " 
      + gnuplot->style() 
      + " lw 2 "
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
   of.open("output/forces_" + vasprun->label + ".data");     
   for (int i=0; i <= bins_sum.size(); i++) {
      if (bins_count[i]>0) {
         if (i*bin_width < smallest_dimension/2.0) {
//            cout << bins_count[i] << endl;
            of << i*bin_width << "\t" << bins_sum[i]/bins_count[i] << "\t" << 0 << "\n";
            //of << i*bin_width << "\t" << bins_sum[i]/bins_count[i] << "\t" << bins_std_dev[i] << "\n";
         }
      }
   }
   of.close();

   
   return 0;

}

int force_bond_projections_wrapper(Configuration *config) {
   screen.status << "Force Projections";
   screen.step << "Bond force projections requested."; 
   
   //Make a gnuplot object.  It takes care of writing the data to a script.
   GnuPlotScript gnuplot ;
   gnuplot.initialise("forces","Force","Separation distance [Angstrom]","Average force [eV/Angstrom]","force.pdf");
   gnuplot.command("set style fill transparent solid 0.40 noborder");
   gnuplot.command("set yrange [-10:10]");
   gnuplot.command("plot 0 lc rgb '#CCCCCC',  ", false );


   for (int i=0; i<config->vaspruns.size(); i++) {
      force_bond_projections(&config->vaspruns[i], config, &gnuplot);
   }

   gnuplot.close();

   config->script_wrapper << "\ngnuplot plot_forces.gnu \n ";

   return 0;

}



#endif
