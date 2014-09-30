#ifndef UTILITY_FUNCTIONS_H
#define UTILITY_FUNCTIONS_H
#include <string>
#include <vector>
#include <stdio.h>
#include <string.h>

using namespace std;


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
}

vector<string> string_split_to_vector2(string s) {
   vector<string> sv;
   string delimiter = " ";
   size_t pos = 0;
   string token;
   while ((pos = s.find(delimiter)) != std::string::npos) {
      token = s.substr(0, pos);
      if (!0==strcmp(delimiter.c_str(), token.c_str())) {
         sv.push_back(token); //std::cout << token << std::endl;
      
      s.erase(0, pos + delimiter.length());
   }
   sv.push_back(s); //std::cout << s << std::endl; 
   return sv;



}

}



#endif
