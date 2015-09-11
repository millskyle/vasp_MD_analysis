#include "readXML.cpp"
#include "read_config.cpp"
#include "utility_functions.h"
#include "screen.h"
#include "msd.cpp"
#include "rdf.cpp"
#include "rho.cpp"
#include "spatial_dist.cpp"
#include "forces.cpp"
#include "PET.cpp"


using namespace std;

//extern ScreenOutput screen;



int main() {
   cout << BOLDRED << "WELCOME!!!" << RESET << "\n";
   Configuration config;
   read_configfile(config);



   VasprunXML* v;
   
   for (int f=0; f<config.vaspruns_labels.size(); f++) {
      screen.section << "READING FILE " + config.vaspruns_labels[f];
      v = &config.vaspruns[f];
      if (readXML(&config.vaspruns[f], &config) == 1)  {
         screen.error << "READING XML FILE FAILED. EXITING\n";
         return 0;
      } else {
         screen.finished << "Successfully read XML file";
      }



      screen.section << "CREATING ATOM SUBSETS";
      for (int i=0; i<config.filter_name_list.size(); i++) {
         screen.step << "Creating set: " + config.filter_name_list[i]; 
         config.atomfilters[config.filter_name_list[i]].execute_filter(v);
         screen.data("   atoms",config.atomfilters[config.filter_name_list[i]].atoms.atoms.timesteps[0].ppp.size() );
      }


   v->atomfilters = config.atomfilters;


}




   screen.section << "PLOTS";
   if (config.rdf) {
      radial_distribution_function_wrapper(&config);
   }

   PET_plots_wrapper(&config);

   if (config.msd) {
      msd_wrapper(&config);
   }
   
   if (config.forces) {
      force_bond_projections_wrapper(&config);
   }
/*
   if (config.rho) {
      global_density(&v, &config);
   }

   */
   if (config.spatial_distribution_projection) {
      spatial_distribution_projection(&config.vaspruns[0], &config);
   }



   config.script_wrapper.close();
   return 0;
   

}


