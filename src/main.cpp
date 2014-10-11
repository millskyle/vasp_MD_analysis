#include "readXML.cpp"
#include "read_config.cpp"
#include "utility_functions.h"
#include "screen.h"
#include "msd.cpp"
#include "rdf.cpp"
#include "rho.cpp"
#include "spatial_dist.cpp"


using namespace std;

//extern ScreenOutput screen;


int main() {
   cout << BOLDRED << "WELCOME!!!" << RESET << "\n";
   FileInfo v;
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

   config.script_wrapper.close();
   return 0;
   

}


