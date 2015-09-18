#ifndef FORCES_CPP
#define FORCES_CPP

#include "data_structure.h"
#include "utility_functions.h"


int force_bond_projections(VasprunXML *vasprun, Configuration *config, GnuPlotScript* gnuplot) {
   if (!config->forces) {cout << "\nForce projection called but not requested in configuration. Exiting"; return 1;}   

   atomType *atomobject0;
   atomType *atomobject1;

   atomobject0 = &(vasprun->atomfilters[config->forces_set_1.c_str()].atoms.atoms);
   atomobject1 = &(vasprun->atomfilters[config->forces_set_2.c_str()].atoms.atoms);
   
   double dx,dy,dz,fx0,fy0,fz0,fx1,fy1,fz1; //components of the vectors between atoms
   double proj0, proj1; //scalar projections
   double distance, total_proj;
   vector<vector<double>> all_data;
   vector<double> thisdata;

   double max_distance = 
     max(
          (pow(vasprun->latt[0][0],2) + pow(vasprun->latt[0][1],2) + pow(vasprun->latt[0][2],2) ), 
          max(
               (pow(vasprun->latt[1][0],2) + pow(vasprun->latt[1][1],2) + pow(vasprun->latt[1][2],2) ), 
               (pow(vasprun->latt[2][0],2) + pow(vasprun->latt[2][1],2) + pow(vasprun->latt[2][2],2) )
             )
        );
     max_distance = sqrt(max_distance)/1.9;

  int total_timesteps = floor(atomobject1->timesteps.size()/config->time_skip) ;

  //make vectors that are the correct size (ie: number of bins) to hold the data
   vector<int> bins_count;
   vector<double> bins_sum ;
   vector<double> bins_square_sum ;
   vector<double> bins_std_dev;
   vector<double> cum_sum; //cumulative sum for integration of force curve
   //calculate the width of each bin 
   double bin_width = max_distance / config->forces_bins;

   for (int i=0; i<config->forces_bins; i++) {
      bins_std_dev.push_back(0);
      bins_count.push_back(0);
      bins_sum.push_back(0.0);
      bins_square_sum.push_back(0.0);
      cum_sum.push_back(0.0);
   }

/* Generate Lennard-Jones force for testing *//*
   for (int i=0; i<config->forces_bins; i++) {
      double x=i*bin_width;      
      bins_sum[i] =1e8*6*(pow(x,6) - 2) / ( pow(x,13));
      bins_count[i] = 1e8;
      bins_square_sum[i] = pow(bins_sum[i],2);
   }

*/

   double min_distance=100000000;

   int number_of_projections = 2*atomobject1->timesteps[0].ppp.size()*atomobject0->timesteps[0].ppp.size()*(atomobject1->timesteps.size());

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
   //if ( 0 == 1 ) {  /* If testing with Lennard Jones, make sure this is uncommented to avoid calculating projections */
      int thisBin;
      cout << "        t = 0         0.00%     " ;
      for (int t=0; t < (atomobject1->timesteps.size() - config->time_skip); t+=config->time_skip ) {   // -2 as the last timestep could be incomplete
         cout << "\r        t = " << t << "         " << t*100.0/(atomobject1->timesteps.size()) << "%       " ;
         cout.flush();
         for (int a=0; a<atomobject1->timesteps[0].ppp.size(); a++) {
   //         if ((atomobject0->timesteps[t].ppp[a][2]) > vasprun->latt[2][2] /2.0) {
            for (int b=0; b<atomobject0->timesteps[0].ppp.size(); b++) {
   //          if ((atomobject0->timesteps[t].ppp[b][2]) > vasprun->latt[2][2] /2.0) {
               
               //vector between the two atoms
               dx = atomobject0->timesteps[t].ppp[b][0] - atomobject1->timesteps[t].ppp[a][0];
               dy = atomobject0->timesteps[t].ppp[b][1] - atomobject1->timesteps[t].ppp[a][1];
               dz = atomobject0->timesteps[t].ppp[b][2] - atomobject1->timesteps[t].ppp[a][2];
               //minimum image correction:
               dx = dx - nint(dx / vasprun->latt[0][0])*vasprun->latt[0][0];
               dy = dy - nint(dy / vasprun->latt[1][1])*vasprun->latt[1][1];
               dz = dz - nint(dz / vasprun->latt[2][2])*vasprun->latt[2][2];
               //distance between both atoms:
               distance = sqrt(dx*dx + dy*dy + dz*dz);

               //get the forces from the atom objects
               fx0 = atomobject0->timesteps[t].fff[b][0];
               fy0 = atomobject0->timesteps[t].fff[b][1];
               fz0 = atomobject0->timesteps[t].fff[b][2];
               fx1 = atomobject1->timesteps[t].fff[a][0];
               fy1 = atomobject1->timesteps[t].fff[a][1];
               fz1 = atomobject1->timesteps[t].fff[a][2];
           
               //project them
               // (F dot r) / |r|
               proj0 = (fx0*dx + fy0*dy + fz0*dz)/distance;
               proj1 = (fx1*dx + fy1*dy + fz1*dz)/distance;

               thisBin = floor(distance / bin_width);

            
               bins_sum[thisBin]+= -proj1;
               bins_sum[thisBin]+= proj0;

               bins_square_sum[thisBin]+= pow(proj1,2);
               bins_square_sum[thisBin]+= pow(proj0,2);
            
               bins_count[thisBin]+=2;

         }
      }
      cout << "\r" << endl;
   }


   for (int i=0; i<bins_square_sum.size(); i++) {
   //   bins_std_dev[i] = sqrt( bins_square_sum[i]/bins_count[i] - pow(bins_sum[i]/bins_count[i],2) )  ;
      bins_std_dev[i] =  sqrt (   ( bins_square_sum[i] - (bins_sum[i]*bins_sum[i])/bins_count[i] ) / (bins_count[i] - 1 )   ) / sqrt(bins_count[i]);
   }


/* go backward through the data and integrate it from infinity to the point */
   for (int i=bins_square_sum.size()-1; i>0; i-=1) {
      if (not(bins_count[i-1]==0)) {
         cum_sum[i-1] = (bins_sum[i-1]/bins_count[i-1])*bin_width + cum_sum[i] ;
      }
   }
   
// gnuplot->command(  
//      + " 'forces_" 
//      + vasprun->label 
//      + ".data' using 1:2 w l ls " 
//      + gnuplot->style() 
//      + " lw 2  title '" 
//      + vasprun->label 
//      + "', "
//      ,false);
   
  //write out the gnuplot command
   gnuplot->command(  
      + " 'forces_" 
      + vasprun->label 
      + ".data' using 1:($2-$3):($2+$3) with filledcurves lc rgb '#FF4444' lw 0 title 'standard error',  " 
      + " '' using 1:2 w l ls " + gnuplot->style() + " lw 1 title '" + vasprun->label + " ("+ config->forces_set_1.c_str() + " - " + config->forces_set_2.c_str() + ") ',"
      + " '' u 1:4 w l ls 14 lw 2 title 'Potential V(r)' ; "
      ,false);

   cout << "HERE" << endl;

   float smallest_dimension;
   smallest_dimension = vasprun->latt[0][0];
   if (vasprun->latt[1][1] < smallest_dimension) {
      smallest_dimension = vasprun->latt[1][1];
   }
   if (vasprun->latt[2][2] < smallest_dimension) {
      smallest_dimension = vasprun->latt[2][2];
   }
  // cout << "HERE" << bins_sum.size() << endl;
  
   
   //write out the data for this element to an element-specific file
   ofstream of;

  // cout << "ofstream of;" << endl;


   of.open("output/forces_" + vasprun->label + ".data");     
   //cout << "FILE OPENED" << endl;
   for (int i=0; i < bins_sum.size(); i++) {
//      cout << "i=" << i <<  endl;
      if (bins_count[i]>0) {
         if (i*bin_width < smallest_dimension/2.0) {
  //          cout << "E" << endl;
//            cout << bins_count[i] << endl;
            //of << i*bin_width << "\t" << bins_sum[i]/bins_count[i] << "\t" << 0 << "\n";
            of << i*bin_width << "\t" << bins_sum[i]/bins_count[i] << "\t" << bins_std_dev[i] << "\t" << -cum_sum[i] << "\n";
         }
      }
   }
   of.close();

//   cout << "HERE" << endl;
   
   return 0;

}

int force_bond_projections_wrapper(Configuration *config) {
   screen.status << "Force Projections";
   screen.step << "Bond force projections requested."; 
   
   //Make a gnuplot object.  It takes care of writing the data to a script.
   GnuPlotScript gnuplot ;
   gnuplot.initialise("forces","Force","Separation distance [Angstrom]","Average force [eV/Angstrom]","force.pdf");
   gnuplot.command("set style fill transparent solid 0.40 noborder");
//   gnuplot.command("set yrange [-10:10]");
   gnuplot.command("plot 0 lc rgb '#CCCCCC' notitle,  ", false );


   for (int i=0; i<config->vaspruns.size(); i++) {
      screen.step << config->vaspruns[i].label; 
      force_bond_projections(&config->vaspruns[i], config, &gnuplot);
   }

   gnuplot.close();

   config->script_wrapper << "\ngnuplot plot_forces.gnu \n ";

   return 0;

}



#endif
