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

static map<string, atomfilter> parse_atom_selections(Configuration& config, const Node& node) {
   const Node& atom_selections = node["atom_sets"];
   map<string, atomfilter> allfilters;
   for (int i=0; i < atom_selections.size(); i++) {
      atomfilter thisFilter;
      parse(atom_selections[i], "name", thisFilter.name);
      parse(atom_selections[i], "filter_type", thisFilter.filter_type);
      parse(atom_selections[i], "criteria", thisFilter.criteria);
      allfilters[thisFilter.name] = thisFilter;
//      allfilters[thisFilter.name].execute_filter();
      config.filter_name_list.push_back(thisFilter.name);
   }
   return allfilters;
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
   parse(plots, "msd_filter", config.msd_filter);
   parse(plots, "msd_reference_D", config.msd_reference_D);
    
   parse(plots, "rdf", config.rdf);
   parse(plots, "rdf_filter_1", config.rdf_filter_1);
   parse(plots, "rdf_filter_2", config.rdf_filter_2);
   parse(plots, "rdf_data_prefix", config.rdf_data_prefix);
   parse(plots, "rdf_bins",config.rdf_bins);
   parse(plots, "rdf_plot_type",config.rdf_plot_type);

   parse(plots, "rho",config.rho);
   parse(plots, "rho_data_prefix", config.rho_data_prefix);
   parse(plots, "rho_filter", config.rho_filter);
   parse(plots, "rho_experimental", config.rho_experimental);

   parse(plots, "sd_projection", config.spatial_distribution_projection);
   parse(plots, "sd_collapse_dimension", config.sd_collapse_dimension);
   parse(plots, "sd_dynamic_filter", config.sd_dynamic_filter);
   parse(plots, "sd_static_filter", config.sd_static_filter);
   parse(plots, "sd_nbins_v", config.sd_nbins_v);
   parse(plots, "sd_nbins_h", config.sd_nbins_h);

   parse(plots, "forces", config.forces);
   parse(plots, "forces_bins", config.forces_bins);
   parse(plots, "forces_set_1", config.forces_set_1);
   parse(plots, "forces_set_2", config.forces_set_2);


   config.atomfilters = parse_atom_selections(config, node );
   


}


bool read_configfile(Configuration& config) {
    screen.section << "Reading configuration file \"config.yaml\" ";
    fstream file("config.yaml", fstream::in);
    if (!file.is_open()) {
        return false;
    }
    Parser parser(file);
    Node root;
    parser.GetNextDocument(root);
    parse_inputfile(config, root);
    screen.finished << "Configuration file \"config.yaml\" read.";
    return true;
}
