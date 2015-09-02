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






int atomic_symbol_to_number(string symbol) {
map <string, int> m;
m["H"]=1;
m["He"]=2;
m["Li"]=3;
m["Be"]=4;
m["B"]=5;
m["C"]=6;
m["N"]=7;
m["O"]=8;
m["F"]=9;
m["Ne"]=10;
m["Na"]=11;
m["Mg"]=12;
m["Al"]=13;
m["Si"]=14;
m["P"]=15;
m["S"]=16;
m["Cl"]=17;
m["Ar"]=18;
m["K"]=19;
m["Ca"]=20;
m["Sc"]=21;
m["Ti"]=22;
m["V"]=23;
m["Cr"]=24;
m["Mn"]=25;
m["Fe"]=26;
m["Co"]=27;
m["Ni"]=28;
m["Cu"]=29;
m["Zn"]=30;
m["Ga"]=31;
m["Ge"]=32;
m["As"]=33;
m["Se"]=34;
m["Br"]=35;
m["Kr"]=36;
m["Rb"]=37;
m["Sr"]=38;
m["Y"]=39;
m["Zr"]=40;
m["Nb"]=41;
m["Mo"]=42;
m["Tc"]=43;
m["Ru"]=44;
m["Rh"]=45;
m["Pd"]=46;
m["Ag"]=47;
m["Cd"]=48;
m["In"]=49;
m["Sn"]=50;
m["Sb"]=51;
m["Te"]=52;
m["I"]=53;
m["Xe"]=54;
m["Cs"]=55;
m["Ba"]=56;
m["La"]=57;
m["Ce"]=58;
m["Pr"]=59;
m["Nd"]=60;
m["Pm"]=61;
m["Sm"]=62;
m["Eu"]=63;
m["Gd"]=64;
m["Tb"]=65;
m["Dy"]=66;
m["Ho"]=67;
m["Er"]=68;
m["Tm"]=69;
m["Yb"]=70;
m["Lu"]=71;
m["Hf"]=72;
m["Ta"]=73;
m["W"]=74;
m["Re"]=75;
m["Os"]=76;
m["Ir"]=77;
m["Pt"]=78;
m["Au"]=79;
m["Hg"]=80;
m["Tl"]=81;
m["Pb"]=82;
m["Bi"]=83;
m["Po"]=84;
m["At"]=85;
m["Rn"]=86;
m["Fr"]=87;
m["Ra"]=88;
m["Ac"]=89;
m["Th"]=90;
m["Pa"]=91;
m["U"]=92;
m["Np"]=93;
m["Pu"]=94;
m["Am"]=95;
m["Cm"]=96;
m["Bk"]=97;
m["Cf"]=98;
m["Es"]=99;
m["Fm"]=100;
m["Md"]=101;
m["No"]=102;
m["Lr"]=103;
m["Rf"]=104;
m["Db"]=105;
m["Sg"]=106;
m["Bh"]=107;
m["Hs"]=108;
m["Mt"]=109;
m["Ds"]=110;
m["Rg"]=111;
m["Cp"]=112;
m["Uut"]=113;
m["Uuq"]=114;
m["Uup"]=115;
m["Uuh"]=116;
m["Uus"]=117;
m["Uuo"]=118;

return m[symbol];


}



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
