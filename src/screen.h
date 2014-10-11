#ifndef SCREEN_CPP
#define SCREEN_CPP
#include <string>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <iostream>


#define RESET   "\033[0m"
#define BLACK   "\033[90m"      /* Black */
#define RED     "\033[91m"      /* Red */
#define GREEN   "\033[92m"      /* Green */
#define YELLOW  "\033[93m"      /* Yellow */
#define BLUE    "\033[94m"      /* Blue */
#define MAGENTA "\033[95m"      /* Magenta */
#define CYAN    "\033[96m"      /* Cyan */
#define WHITE   "\033[97m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */




struct bWHITE {
    int status_count = 1;
    template<typename T> bWHITE operator << (T&& s)  {
              std::cout << BOLDWHITE << status_count << ". "<<  s << RESET << "\n";
              status_count++;
              return *this;
      };
};

struct bRED {
    template<typename T> bRED operator << (T&& s)  {
              std::cout << BOLDRED << "ERROR:\t" << s << RESET << "\n";
              return *this;
      };
};


struct oYELLOW {
    template<typename T> oYELLOW operator << (T&& s)  {
              std::cout << YELLOW << "   " << s << RESET << "\n";
              return *this;
      };
};


struct bMAGENTA {
    template<typename T> bMAGENTA operator << (T&& s)  {
              std::cout << BOLDMAGENTA << "------------------------------------------\n";
              std::cout << "   "  << s << "\n";
              std::cout << "------------------------------------------" <<RESET << "\n";
              return *this;
      };
};

struct bGREEN {
    template<typename T> bGREEN operator << (T&& s)  {
              std::cout << BOLDGREEN << s << RESET << "\n";
              return *this;
      };
};


struct ScreenOutput {
   bRED error;
   bWHITE status;
   oYELLOW step;
   bGREEN finished;
   bMAGENTA section;


   template <typename T>
   void data(string name, T&& dd) {
      std::cout << "     " << BOLDCYAN << name << ":\t" << CYAN << dd << RESET << "\n";
   };

} screen ;















#endif
