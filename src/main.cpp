#include "readXML.cpp"
#include "read_config.cpp"
#include "utility_functions.h"


#include "msd.cpp"
#include "rdf.cpp"

using namespace std;





int main() {
   cout << "\n Starting XML read"<<endl;
   FileInfo v;
   Configuration config;
   read_configfile(config);
   cout << "\n done"<<endl;

   v.input_filename="vasprun.xml";
   readXML(&v);


   if (config.msd) {
      mean_square_displacement(&v, &config);
   }

 if (config.rdf) {
      radial_distribution_function(&v, &config);
   }

   config.script_wrapper.close();
   return 0;
}


