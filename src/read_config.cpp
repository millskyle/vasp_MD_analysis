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



static int parse_input_file_choices(Configuration &config, const Node& node) {
   const Node& vaspruns = node["vaspruns"];
   string thisLabel;
   string thisLocation;
   for (int i=0; i<vaspruns.size(); i++) {
      parse(vaspruns[i], "label", thisLabel);
      parse(vaspruns[i], "file", thisLocation);
      config.vaspruns_labels.push_back(thisLabel);
      VasprunXML v;
      v.input_filename = thisLocation;
      v.label = thisLabel;
      config.vaspruns.push_back(v);
   }
   return 0;
}




void parse_inputfile(Configuration& config, const Node& node) {
   const Node& files = node["files"];
   parse(files, "log",config.log_file_location);
     config.log.open(config.log_file_location);
   parse(files, "script_wrapper", config.script_wrapper_location);
     config.script_wrapper.open(config.script_wrapper_location);
     config.script_wrapper << "#!/bin/sh\n";


   const Node& plots = node["plots"];
   parse(plots, "p_t_start_timestep", config.p_t_start_timestep);
   parse(plots, "e_t_start_timestep", config.e_t_start_timestep);


   parse(plots, "msd", config.msd);
   parse(plots, "msd_data_prefix", config.msd_data_prefix);
   parse(plots, "msd_filter", config.msd_filter);
   parse(plots, "msd_reference_D", config.msd_reference_D);
   parse(plots, "msd_reference_shift", config.msd_reference_shift);
    
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

   parse(plots, "sd_wrap", config.sd_wrap);
   parse(plots, "sd_h_wrap_dist", config.sd_h_wrap_dist);
   parse(plots, "sd_v_wrap_dist", config.sd_v_wrap_dist);
   parse(plots, "sd_h_wrap_times", config.sd_h_wrap_times);
   parse(plots, "sd_v_wrap_times", config.sd_v_wrap_times);




   parse(plots, "forces", config.forces);
   parse(plots, "forces_bins", config.forces_bins);
   parse(plots, "forces_set_1", config.forces_set_1);
   parse(plots, "forces_set_2", config.forces_set_2);
   parse(plots, "forces_yrange", config.forces_yrange);


   parse(files, "time_start", config.time_start);
   parse(files, "time_end", config.time_end);
   parse(files, "time_skip", config.time_skip);
   if (config.time_skip <=0 ) {
      config.time_skip = 1;
   }
   
   cout << config.time_start << "   " << config.time_end << endl;


   config.atomfilters = parse_atom_selections(config, node );

   parse_input_file_choices(config, node);
   

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
