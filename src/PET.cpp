#ifndef PET_H
#define PET_H

#include "data_structure.h"
#include "utility_functions.h"



//Make a gnuplot object.  It takes care of writing the data to a script.
//GnuPlotScript gnuplot ;
GnuPlotScript gnuplotE ;
















int PET_plots(VasprunXML *vasprun, Configuration *config) {
   

   //For each atom in the requested atom types
//2015-09-01   for (int atomname=0; atomname < config->msd_atoms.size(); atomname++) {
      //this pointer will point to the atomType object for this type of atom 
   screen.step << "Pressure plot vertical range set to range of values for t > " + to_string((config->p_t_start_timestep * vasprun->dt)/1000.0) + " picoseconds (timestep " + to_string(config->p_t_start_timestep) + ")." ;
   screen.step << "Energy plot vertical range set to range of values for t > " + to_string((config->e_t_start_timestep * vasprun->dt)/1000.0) + " picoseconds (timestep " + to_string(config->e_t_start_timestep) + ")." ;

   ofstream et;
   et.open("output/e-t_" + vasprun->label + ".data");

   ofstream pt;
   pt.open("output/p-t_" + vasprun->label + ".data");


   double Emin = 100000000.0;
   double Emax = -10000000.0;
   double E = 0;

   double Pmin = 100000000.0;
   double Pmax = -10000000.0;
   double P = 0;

   for (int t=0; t<vasprun->timesteps.size(); t++) {
      P = vasprun->timesteps[t].pressure();
      E = vasprun->timesteps[t].energy["e_fr_energy"]; 

      if (t > config->p_t_start_timestep) {
         Pmin = min(Pmin, P);
         Pmax = max(Pmax, P);
      }
      if (t > config->e_t_start_timestep) {
         Emin = min(Emin, E);
         Emax = max(Emax, E);
      }
         
      //write out the data for this element to an element-specific file
      
      pt << (t*vasprun->dt)/1000.0 << "\t" << P << "\n" ;
      et << (t*vasprun->dt)/1000.0 << "\t" << E << "\n" ;

   }

   pt.close();
   et.close();

   //write out the gnuplot command, scaling the x-axis increment by the timestep to get it in picoseconds
   gnuplotE.command("'e-t_" + vasprun->label + ".data' using 1:2 with lines ls " + gnuplotE.style() + " lw 3  , " ,false);
   gnuplot.command("'p-t_" + vasprun->label + ".data' using 1:2 with lines ls " + gnuplot.style() + " lw 3  , " ,false);
   return 0;

}

int PET_plots_wrapper(Configuration *config) {
   
   screen.status << "Plotting pressure, energy as a function of time";
  
   gnuplotE.initialise("E-t","Energy vs. Time","Time, [picoseconds]","Total Energy [eV]","E-t.pdf");
   gnuplot.initialise("P-t","Pressure vs. Time","Time, [picoseconds]","Pressure [kB]","P-t.pdf");

   gnuplotE.command("set xrange [0:]");
   gnuplot.command("set xrange [0:]");

   gnuplotE.command("plot ",false);
   gnuplot.command("plot ",false);

   for (int i=0; i<config->vaspruns.size(); i++) {
      PET_plots( &config->vaspruns[i], config);
   }


   gnuplot.close();
   gnuplotE.close();


   //add a command to the global plot script to make the msd plots
   config->script_wrapper << "\ngnuplot plot_P-t.gnu \n" ;   
   config->script_wrapper << "\ngnuplot plot_E-t.gnu \n" ;   
   

   return 0;

}


#endif
