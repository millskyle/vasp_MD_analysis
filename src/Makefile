

main: readXML.cpp data_structure.h msd.cpp rdf.cpp rho.cpp utility_functions.h read_config.cpp spatial_dist.cpp
	$(CXX) -std=c++11 main.cpp tinyxml/tinyxml.cpp tinyxml/tinyxmlerror.cpp tinyxml/tinyxmlparser.cpp tinyxml/tinystr.cpp -Iyaml-cpp/include/ -Lyaml-cpp/ -lyaml-cpp -Ltinyxml/ -orun

clang: readXML.cpp data_structure.h  
	clang++ -std=c++11 main.cpp tinyxml/tinyxml.cpp tinyxml/tinyxmlerror.cpp tinyxml/tinyxmlparser.cpp tinyxml/tinystr.cpp -Iyaml-cpp/include/ -Lyaml-cpp/ -lyaml-cpp -Ltinyxml/ -orun

clean:
	rm run



