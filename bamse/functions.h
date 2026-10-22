#ifndef _ERRSTOP_H
#define _ERRSTOP_H
#include <stdio.h>
#include <stdlib.h>

#define err_stop(STR) { printf(__FILE__ " %5d:" STR "\n", __LINE__); exit(-1); }
#define err_start(STR) { printf(__FILE__ " %5d:" STR "\n", __LINE__);}

#endif
#include <string.h> // some compilers require <mem.h> instead

// If your system doesn't have a rotate function for 32 bits integers,
// then use the definition below. If your system has the _lrotl function 
// then remove this.
uint32 _lrotl (uint32 x, int r) {
  return (x << r) | (x >> (sizeof(x)*8-r));}


// constructor:
TRanrotWGenerator::TRanrotWGenerator(uint32 seed) {
  RandomInit(seed);
  // detect computer architecture
  randbits[2] = 0; randp1 = 1.0;
  if (randbits[2] == 0x3FFF) Architecture = EXTENDEDPRECISIONLITTLEENDIAN;
  else if (randbits[1] == 0x3FF00000) Architecture = LITTLEENDIAN;
  else if (randbits[0] == 0x3FF00000) Architecture = BIGENDIAN;
  else Architecture = NONIEEE;}

// constructor:
TRanrotWGenerator::TRanrotWGenerator() {
  uint32 seed=time(0);
 // uint32 seed=44444;
  RandomInit(seed);
  // detect computer architecture
  randbits[2] = 0; randp1 = 1.0;
  if (randbits[2] == 0x3FFF) Architecture = EXTENDEDPRECISIONLITTLEENDIAN;
  else if (randbits[1] == 0x3FF00000) Architecture = LITTLEENDIAN;
  else if (randbits[0] == 0x3FF00000) Architecture = BIGENDIAN;
  else Architecture = NONIEEE;}

uint32 TRanrotWGenerator::BRandom() {
  // generate next random number
  uint32 y, z;
  // generate next number
  z = _lrotl(randbuffer[p1][0], R1) + randbuffer[p2][0];
  y = _lrotl(randbuffer[p1][1], R2) + randbuffer[p2][1];
  randbuffer[p1][0] = y; randbuffer[p1][1] = z;
  // rotate list pointers
  if (--p1 < 0) p1 = KK - 1;
  if (--p2 < 0) p2 = KK - 1;
  // perform self-test
  if (randbuffer[p1][0] == randbufcopy[0][0] &&
    memcmp(randbuffer, randbufcopy[KK-p1], 2*KK*sizeof(uint32)) == 0) {
      // self-test failed
      if ((p2 + KK - p1) % KK != JJ) {
        // note: the way of printing error messages depends on system
        // In Windows you may use FatalAppExit
        printf("Random number generator not initialized");}
      else {
        printf("Random number generator returned to initial state");}
      exit(1);}
  randbits[0] = y;
  randbits[1] = z;
  return y;
  }


long double TRanrotWGenerator::Random() {
  // returns a random number between 0 and 1.
  uint32 z = BRandom();  // generate 64 random bits
  switch (Architecture) {
  case EXTENDEDPRECISIONLITTLEENDIAN:
    // 80 bits floats = 63 bits resolution
    randbits[1] = z | 0x80000000;  break;
  case LITTLEENDIAN:
    // 64 bits floats = 52 bits resolution
    randbits[1] = (z & 0x000FFFFF) | 0x3FF00000;  break;
  case BIGENDIAN:
    // 64 bits floats = 52 bits resolution
    randbits[0] = (randbits[0] & 0x000FFFFF) | 0x3FF00000;  break;
  case NONIEEE: default:
    // not a recognized floating point format. 32 bits resolution
    return (double)z * (1./((double)(uint32)(-1L)+1.));}
  return randp1 - 1.0;}


int TRanrotWGenerator::IRandom(int min, int max) {
  // get integer random number in desired interval
  int iinterval = max - min + 1;
  if (iinterval <= 0) return 0x80000000;  // error
  int i = int(iinterval * Random());      // truncate
  if (i >= iinterval) i = iinterval-1;
  return min + i;}


void TRanrotWGenerator::RandomInit (uint32 seed) {
  // this function initializes the random number generator.
  int i, j;

  // make random numbers and put them into the buffer
  for (i=0; i<KK; i++) {
    for (j=0; j<2; j++) {
      seed = seed * 2891336453UL + 1;
      randbuffer[i][j] = seed;}}
  // set exponent of randp1
  randbits[2] = 0; randp1 = 1.0;
  // initialize pointers to circular buffer
  p1 = 0;  p2 = JJ;
  // store state for self-test
  memcpy (randbufcopy, randbuffer, 2*KK*sizeof(uint32));
  memcpy (randbufcopy[KK], randbuffer, 2*KK*sizeof(uint32));
  // randomize some more
  for (i=0; i<31; i++) BRandom();
}


///////////////////
/*returns the a random normail distributet number*/
double random_normal(double esp,double var)
{
  double u1=rg.Random();
  double u2=rg.Random();
  double z=sqrt(-2.0*log(u1))*cos(2*M_PI*u2);
  return esp+var*z;
}

//////////////////////////////////////////
/*returns the binomial coefficient*/
int binom(int n, int k) {
   if ( k==0 || k==n )
      return 1;
   else 
      return binom(n-1,k-1) 
           + binom(n-1,k) ;
}

//float random_normal_min_max
#include <cmath>
#ifndef PI
#define PI 3.141592653589793238462643
#endif

double n(double z){//normal distribution function
  return(1.0/sqrt(2.0*PI))*exp(-0.5*z*z);
}



float FRandom(float themin,float themax){
      return ((themax-themin)*rg.Random()+themin);
}

float intNormal(float current_risk_mean, float new_risk_mean,float prior_mean, float var){
  float x=(Kvadrat(new_risk_mean)-Kvadrat(current_risk_mean)-2*prior_mean*(new_risk_mean-current_risk_mean))/(2*Kvadrat(var));
  return((-x));     
}

int is_in_array(int item,int length,int array[]){
  for(int tal=0;tal<length;tal++){
    if(item==array[tal]) return(1);
  }
return(0);
}
      
void sampler(int size,int number,int chosen[]){
  int current=0;
  for(int tal=0;tal<number;tal++){
    while(1){
      current=rg.IRandom(0,size-1);
      chosen[tal]=current;
      if(!is_in_array(current,tal,chosen)) break;
    }
    //cout<<current<<endl;
  }
}

void divider(int number,int chosen[]){
  int current=0;
  for(int tal=0;tal<number;tal++){
      current=rg.IRandom(1,3);
      if(current<=2)
      chosen[tal]=current;
      else
      chosen[tal]=3;
  }
}

void print_array(int size,int array[]){
  for(int tal=0;tal<size;tal++){
    cout<<array[tal]<<" ";
  }
  cout<<endl;
}

 void print_array_float(int size,float array[]){
  for(int tal=0;tal<size;tal++){
    cout<<array[tal]<<" ";
  }
  cout<<endl;
}


void SortInt(int x[], int n) {
    for (int pass=0; pass<n-1; pass++) {
        int potentialSmallest = pass;  // assume this is smallest

        //--- Look over remaining elements to find smallest.
        for (int i=pass+1; i<n; i++) {
            if (x[i] < x[potentialSmallest]) {
                //--- Remember index for latter swap.
                potentialSmallest = i;
            }
        }
        
        //--- Swap smallest remaining element
        int temp = x[pass];
        cout<<"emp "<<potentialSmallest<<endl;
        x[pass] = x[potentialSmallest];
        x[potentialSmallest] = temp;
    }
}

void SortFloat(float x[], int n) {
    for (int pass=0; pass<n-1; pass++) {
        int potentialSmallest = pass;  // assume this is smallest

        //--- Look over remaining elements to find smallest.
        for (int i=pass+1; i<n; i++) {
            if (x[i] < x[potentialSmallest]) {
                //--- Remember index for latter swap.
                potentialSmallest = i;
            }
        }
        
        //--- Swap smallest remaining element
        float temp = x[pass];
        cout<<"emp "<<potentialSmallest<<endl;
        x[pass] = x[potentialSmallest];
        x[potentialSmallest] = temp;
    }
}

void OrderFloat2(float x[], int n,int order[]) {
    for (int pass=0; pass<n; pass++)
      order[pass]=-1;
    for (int pass=0; pass<n; pass++) {
        int number=0; 
        for (int i=0; i<n; i++) {
            if (x[i] < x[pass]) {
                number++;
            }
        }
        if(is_in_array(number,n,order))
          order[pass]=number+1;
        else
        order[pass]=number;
    }
    for(int tal=0;tal<n-1;tal++){
      if(order[tal]==order[tal+1]){
        order[tal+1]=order[tal+1]-1;   
      }
    }
}

void OrderFloat(float x[], int n,int order[]) {
  int *y;
  y = new int[n];
  OrderFloat2(x,n,y);
  for (int pass=0; pass<n; pass++)
  order[pass]=-1;
    for (int pass=0; pass<n; pass++) {
        int number=0; 
        for (int i=0; i<n; i++) {
            if (y[i] > y[pass]) {
                number++;
            }
        }
        order[number]=pass;
    }
    delete[] y;
}




//******************************************************************************

double d_uniform_01 ( int *seed )

//******************************************************************************
//
//  Purpose:


//
//    D_UNIFORM_01 is a portable pseudorandom number generator.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Modified:
//
//    11 August 2004
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, L E Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double D_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}

//*********************************************************************

void i_swap ( int *i, int *j )

//*********************************************************************
//
//  Purpose:
//
//    I_SWAP switches two integer values.
//
//  Modified:
//
//    07 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *I, *J.  On output, the values of I and
//    J have been interchanged.
//
{
  int k;

  k = *i;
  *i = *j;
  *j = k;
 
  return;
}

//********************************************************************

int i_uniform ( int b, int c, int *seed )

//********************************************************************
//
//  Purpose:
//
//    I_UNIFORM returns a random integer in a given range.
//
//  Modified:
//
//    03 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int B, C, the minimum and maximum values acceptable
//    for I.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int I_UNIFORM, the randomly chosen integer.
//
{
  int i;
//
  i = b + ( int ) ( d_uniform_01 ( seed ) * ( double ) ( c + 1 - b ) );
//
//  In case of oddball events at the boundary, enforce the limits.
//
  if ( i < b )
  {
    i = b;
  }

  if ( c < i )
  {
    i = c;
  }

  return i;
}

//********************************************************************

void perm_uniform ( int n, int *seed, int p[])
//********************************************************************
//    PERM_UNIFORM selects a random permutation of N objects.
//    Input, int N, the number of objects to be permuted.
//    Input/output, int *SEED, a seed for the random number generator.
//    Output, int PERM_UNIFORM[N], a permutation of (1,, 1, ..., N).
{
  int i;
  int j;
//  int *p;
//  p = new int[n];
  for ( i = 1; i <= n; i++ )
  {
    p[i-1] = i;
  }
  for ( i = 1; i <= n; i++ )
  {
    j = i_uniform ( i, n, seed );
    i_swap ( &p[i-1], &p[j-1] );
  }
 // return p;
}
//**************************************
float max_array(float cc[],int L){
  float tmax=cc[0];
  for(int tal =1;tal<L;tal++)
    if(cc[tal]>tmax)
      tmax=cc[tal];
  return(tmax);

}

//**********************************************************
int gibbs_pick(float like[], int size){
  float sum=0;
  float max_like=max_array(like,size);
  for(int num=0;num<size;num++){
    if(like[num]>max_like-20){
      like[num]=exp(like[num]-max_like);
      sum+=like[num];
    }
    else
      like[num]=0;
    if(like[num]>1)
      cout<<"like="<<like[num]<<" num="<<num<<endl;
  }
  float ran=FRandom(0,1);
  float cumsum=0;


  for(int num=0;num<size;num++){
    cumsum+=like[num]/sum;

    if(cumsum>1.00001){
      cout<<"cumsum="<<cumsum<<" size="<<size<<" like[0]="<<like[0]<<" sum="<<sum<<" chose="<<num<<endl;
      err_stop("error in opick gibbs")
    }
    if(cumsum>ran)
      return(num);
  }
  cout<<"cumsum="<<cumsum<<" size="<<size<<" like[0]="<<like[0]<<" sum="<<sum<<endl;
  err_stop("gibbs pick does dot sum to one");
}
