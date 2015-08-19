

main: readXML.cpp data_structure.h msd.cpp rdf.cpp rho.cpp utility_functions.h read_config.cpp spatial_dist.cpp
	$(CXX) -std=c++11 main.cpp pugixml-1.6/src/pugixml.cpp -Ipugixml-1.6/src -Iyaml-cpp/include/ -Lyaml-cpp/ -lyaml-cpp -Ltinyxml/ -orun

clang: readXML.cpp data_structure.h  
	clang++ -std=c++11 main.cpp pugixml-1.6/src/pugixml.cpp -Ipugixml-1.6/src  -Iyaml-cpp/include/ -Lyaml-cpp/ -lyaml-cpp -Ltinyxml/ -orun

clang_test: readXML.cpp data_structure.h  
	clang++ -std=c++11 main.cpp pugixml-1.6/src/pugixml.cpp -Ipugixml-1.6/src  -Iyaml-cpp/include/ -Lyaml-cpp/ -lyaml-cpp -Ltinyxml/ -orun && ./run && (cd output ; sh make_all_plots.sh ; cd .. )


clean:
	rm run



