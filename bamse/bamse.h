#define the_max_riskset 15
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string>
#include<cstdlib>
#include <math.h>
#include <ctime>
#include<stdio.h>
#define Kvadrat(X) ((X) * (X))
#define MAX(X,Y) (X > Y ? X : Y)
#define MIN(X,Y) (X < Y ? X : Y)
#include <algorithm>
using namespace std;
const char *the_filename;
const char *the_optionfile;
const char *the_missingfile;


int min_riskset=0; 
int the_max_active=0;
int the_min_active=0; //not implemented
//int relaxed=1; //1: starts with all ind in riskset 1
float max_u=0.01;
float active_prior_p=0; //the parameters for geometric distributed number of active parameteres
float riskset_prior_p=0; //the parameters for geometric distributed number of risk sets
//int antal_variable=antal_snp+antal_env+antal_trait;
float the_max_varians=0;
float var_change=0.05;
float kappa_multi=0;
bool use_midrange =0; //if 1 used the midrange. else used the average mean
int print=0;
float start_varians=0; //sd

const double Pi = 3.14159265358979323846;



    int sample_rate=100;
    int runtime=50000;
    float over=0.05;                     //default=0.05
    int burnin=0;
    int das_burnin=0;
    int mean_update_ratio=10;
    int flip_rate;                       //default=0.01
    int env_change_rate;                 //default=0.20
////////////////////////////////////////////////////////////////////////////


#ifndef RANDOMC_H
#define RANDOMC_H

#include <math.h>
#include <assert.h>
#include <stdio.h>

// Define 32 bit signed and unsigned integers.
// Change these definitions, if necessary, on 64 bit computers
typedef   signed long int32;     
typedef unsigned long uint32;     

class TRandomMersenne {                // encapsulate random number generator
  #if 0
    // define constants for MT11213A:
    // (32 bit constants cannot be defined as enum in 16-bit compilers)
    #define MERS_N   351
    #define MERS_M   175
    #define MERS_R   19
    #define MERS_U   11
    #define MERS_S   7
    #define MERS_T   15
    #define MERS_L   17
    #define MERS_A   0xE4BD75F5
    #define MERS_B   0x655E5280
    #define MERS_C   0xFFD58000
  #else    
    // or constants for MT19937:
    #define MERS_N   624
    #define MERS_M   397
    #define MERS_R   31
    #define MERS_U   11
    #define MERS_S   7
    #define MERS_T   15
    #define MERS_L   18
    #define MERS_A   0x9908B0DF
    #define MERS_B   0x9D2C5680
    #define MERS_C   0xEFC60000
  #endif
  public:
  TRandomMersenne(uint32 seed) {       // constructor
    RandomInit(seed);}
  void RandomInit(uint32 seed);        // re-seed
  void RandomInitByArray(uint32 seeds[], int length); // seed by more than 32 bits
  int IRandom(int min, int max);       // output random integer
  double Random();                     // output random float
  uint32 BRandom();                    // output random bits
  private:
  uint32 mt[MERS_N];                   // state vector
  int mti;                             // index into mt
  enum TArch {LITTLEENDIAN, BIGENDIAN, NONIEEE};
  TArch Architecture;                  // conversion to float depends on computer architecture
  };    

class TRanrotBGenerator {              // encapsulate random number generator
  enum constants {                     // define parameters
    KK = 17, JJ = 10, R1 = 13, R2 =  9};
  public:
  void RandomInit(uint32 seed);        // initialization
  int IRandom(int min, int max);       // get integer random number in desired interval
  double Random();                     // get floating point random number
  TRanrotBGenerator(uint32 seed);      // constructor
  protected:
  int p1, p2;                          // indexes into buffer
  uint32 randbuffer[KK];               // history buffer
  uint32 randbufcopy[KK*2];            // used for self-test
  enum TArch {LITTLEENDIAN, BIGENDIAN, NONIEEE};
  TArch Architecture;                  // conversion to float depends on computer architecture
};


class TRanrotWGenerator {              // encapsulate random number generator
  enum constants {                     // define parameters
    KK = 17, JJ = 10, R1 = 19, R2 =  27};
  public:
  void RandomInit(uint32 seed);        // initialization
  int IRandom(int min, int max);       // get integer random number in desired interval
  long double Random();                // get floating point random number
  uint32 BRandom();                    // output random bits  
  TRanrotWGenerator(uint32 seed);      // constructor
  TRanrotWGenerator();      // constructor
  protected:
  int p1, p2;                          // indexes into buffer
  union {                              // used for conversion to float
    long double randp1;
    uint32 randbits[3];};  
  uint32 randbuffer[KK][2];            // history buffer
  uint32 randbufcopy[KK*2][2];         // used for self-test
  enum TArch {LITTLEENDIAN, BIGENDIAN, NONIEEE, EXTENDEDPRECISIONLITTLEENDIAN};
  TArch Architecture;                  // conversion to float depends on computer architecture
};

class TRandomMotherOfAll {             // encapsulate random number generator
  public:
  void RandomInit(uint32 seed);        // initialization
  int IRandom(int min, int max);       // get integer random number in desired interval
  double Random();                     // get floating point random number
  TRandomMotherOfAll(uint32 seed);     // constructor
  protected:
  double x[5];                         // history buffer
  };

#endif




//int theseed=1159777120;
int theseed=time(0);
int32 seed = theseed; //sets seed random using the clock
TRanrotWGenerator rg(seed); //generates the object rg
#include"functions.h"
#define RAND_UNIFORME (double)random()/(double)RAND_MAX
#include"bamseclass.cpp"
