#include "readXML.cpp"
#include "read_config.cpp"
#include "utility_functions.h"
#include "screen.h"
#include "msd.cpp"
#include "rdf.cpp"
#include "rho.cpp"
#include "spatial_dist.cpp"
#include "forces.cpp"


using namespace std;

//extern ScreenOutput screen;


int main() {
   cout << BOLDRED << "WELCOME!!!" << RESET << "\n";
   VasprunXML v;
   Configuration config;
   read_configfile(config);

   v.input_filename="vasprun.xml";

   screen.section << "READING FILE";
   if (readXML(&v)==1)  {
      screen.error << "READING XML FILE FAILED. EXITING\n";
      return 0;
   } else {
      screen.finished << "Successfully read XML file";
   }

   screen.section << "CREATING ATOM SUBSETS";
   for (int i=0; i<config.filter_name_list.size(); i++) {
      screen.step << "Creating set: " + config.filter_name_list[i]; 
      config.atomfilters[config.filter_name_list[i]].execute_filter(&v);
      screen.data("      atoms",config.atomfilters[config.filter_name_list[i]].atoms.atoms.timesteps[0].ppp.size() );
   }



   screen.section << "PLOTS";
   if (config.msd) {
      mean_square_displacement(&v, &config);
   }
   if (config.rdf) {
      radial_distribution_function(&v, &config);
   }
   if (config.rho) {
      global_density(&v, &config);
   }
   if (config.spatial_distribution_projection) {
      spatial_distribution_projection(&v, &config);
   }
   if (config.forces) {
      force_bond_projections(&v, &config);
   }


   config.script_wrapper.close();
   return 0;
   

}


