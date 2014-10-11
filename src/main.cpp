#include "readXML.cpp"
#include "read_config.cpp"
#include "utility_functions.h"


#include "msd.cpp"
#include "rdf.cpp"
#include "rho.cpp"
#include "spatial_dist.cpp"

using namespace std;


int main() {
   cout << "\n Starting XML read"<<endl;
   FileInfo v;
   Configuration config;
   read_configfile(config);
   cout << "\n done"<<endl;

   v.input_filename="vasprun.xml";
   readXML(&v);

   //if (readXML(&v)==0)  {
//      cout << "\nERROR READING XML FILE. EXITING\n";
//      return 0;
//   }


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


