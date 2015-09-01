#ifndef UTILITY_FUNCTIONS_H
#define UTILITY_FUNCTIONS_H
#include <string>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>


using namespace std;

typedef vector<float> threevector;




int count_atoms_in_box_pointer(vector<threevector>* positions, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) { 
   int total=0;
   for (int a=0; a<positions->size();a++) {
      if ((*positions)[a][0] <= xmax && (*positions)[a][0] > xmin &&  (*positions)[a][1] <= ymax &&
          (*positions)[a][1] > ymin &&  (*positions)[a][2] <= zmax && (*positions)[a][2] > zmin)
       {
          total +=1;
       }
   }

return total;
}

int count_atoms_in_box(vector<threevector> positions, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) { 
   int total=0;
   for (int a=0; a<positions.size();a++) {
      if ((positions)[a][0] <= xmax && (positions)[a][0] > xmin &&  (positions)[a][1] <= ymax &&
          (positions)[a][1] > ymin &&  (positions)[a][2] <= zmax && (positions)[a][2] > zmin)
       {
          total +=1;
       }
   }

return total;
}



string string_replace( string old, string nnew, string subject) {
   subject.replace(subject.find(old), old.length(), nnew);
   return subject;
}


string vec2str(vector<string> v, string delimiter=" ") {
   string s;
   for (int i=0; i<v.size(); i++) {
      s.append(v[i]);
      s.append(" ");
   }
   return s;
}
      



vector<string> str2vec(string s, string delimiter=" ,/;") {
   vector<string> sv;
   stringstream stringStream(s);
   string line;
   while(getline(stringStream, line)) 
   {
      size_t prev = 0, pos;
      while ((pos = line.find_first_of(delimiter, prev)) != string::npos)
      {
         if (pos > prev) {
            sv.push_back(line.substr(prev, pos-prev));
         }
        prev = pos+1;
      }
      if (prev < line.length()) {
         sv.push_back(line.substr(prev, std::string::npos));
      }
return sv;
}
return sv;
}

int nint(float x) {
   return floor(x+0.5);
}





// trim from start
static inline std::string ltrim(std::string s) {
           s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
                   return s;
}

// trim from end
static inline std::string rtrim(std::string s) {
           s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
                   return s;
}

// trim from both ends
static inline std::string trim(std::string s) {
           return ltrim(rtrim(s));
}






  















/*
vector<string> string_split_to_vector2(string s) {
   vector<string> sv;
   string delimiter = " ";
   size_t pos = 0;
   string token;
   while ((pos = s.find(delimiter)) != std::string::npos) {
      token = s.substr(0, pos);
      if (0!=strcmp(delimiter.c_str(), token.c_str())) {
         sv.push_back(token); //std::cout << token << std::endl;
      
      s.erase(0, pos + delimiter.length());
   }
   sv.push_back(s); //std::cout << s << std::endl; 
   return sv;
}

}


*/



struct Conversions {
      double angstroms_per_cm = 1e8;
      double cm_per_angstrom = 1e-8;

      double m2A = 1e10;
      double A2m = 1/1e10;
      
      double to_micro=1e6;
      double to_nano= 1e9;
      double to_pico=1e12;
      double to_femto=1e15;
   
      double atoms_per_mole = 6.0221412927e23;
      double moles_per_gram = 1.0 / atoms_per_mole ;
      double pico2femto = 1000;
      double femto2pico = 1.0/pico2femto;

  Conversions(){
  }
};





     








#endif
