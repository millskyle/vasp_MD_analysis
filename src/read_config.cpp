#include <fstream>
#include <string>
#include "yaml-cpp/yaml.h"
#include <iostream>
#include "utility_functions.h"
#include "data_structure.h"

using namespace std;
using namespace YAML;

template<class T>
static void parse(const Node& node, const char* key, T& value) {
    if (node.FindValue(key)) {
        node[key] >> value;
    } else {
        printf("CONFIG: '%s' was not found!\n", key);
    }
}


void parse_inputfile(Configuration& config, const Node& node) {
   const Node& files = node["files"];
   parse(files, "log",config.log_file_location);
     config.log.open(config.log_file_location);
   parse(files, "script_wrapper", config.script_wrapper_location);
     config.script_wrapper.open(config.script_wrapper_location);
     config.script_wrapper << "#!/bin/sh\n";


   const Node& plots = node["plots"];
   parse(plots, "msd", config.msd);
   parse(plots, "msd_data_prefix", config.msd_data_prefix);
   parse(plots, "msd_atoms", config.tempstr);
      config.msd_atoms = str2vec(config.tempstr);
    
   parse(plots, "rdf", config.rdf);
   parse(plots, "rdf_data_prefix", config.rdf_data_prefix);
   parse(plots, "rdf_bins",config.rdf_bins);
   parse(plots, "rdf_cut_half_lv",config.rdf_cut_half_lv);
   parse(plots, "rdf_plot_type",config.rdf_plot_type);
   parse(plots, "rdf_atoms", config.tempstr);
      config.rdf_atoms = str2vec(config.tempstr);

   parse(plots, "rho",config.rho);
   parse(plots,"rho_atoms",config.tempstr);
      config.rho_atoms = str2vec(config.tempstr);
   parse(plots, "rho_data_prefix", config.rho_data_prefix);

   parse(plots, "spatial_distribution", config.spatial_distribution);
   parse(plots, "collapse_dimension", config.collapse_dimension);
   parse(plots, "lattice_atoms", config.tempstr);
      config.lattice_atoms = str2vec(config.tempstr);
   parse(plots, "liquid_atoms", config.tempstr);
      config.liquid_atoms = str2vec(config.tempstr);
   parse(plots, "nbins_x", config.nbins_x);
   parse(plots, "nbins_y", config.nbins_y);



}




bool read_configfile(Configuration& config) {
    fstream file("config.yaml", fstream::in);
    if (!file.is_open()) {
        return false;
    }
    Parser parser(file);
    Node root;
    parser.GetNextDocument(root);
    parse_inputfile(config, root);
    cout << "Configuration file \"config.yaml\" read.\n\n";
    return true;
}
