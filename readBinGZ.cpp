
//  g++ -o readBinGz  -lz -O3 readBinGZ.cpp 

// large data set
// ./readBinGz -f /home/albrecht/gitfuck/test/ceu8gl2.glf.gz 

#include <cstdio> // NULL
#include <cstring> // strcmp
#include <zlib.h> //gzFile ect. 
#include <cstdlib> //exit
#include <ctime> //time


int bufSize;


/*
double read-in
first get size (first read)
then read in data
NB!: does not work for files larger than can be represented by unsigned int :-(

*/
size_t readGZbin(const char* fname,double * &data){
  time_t t2=time(NULL);

  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  int temp=0;
  size_t nSites=0;
  double *tmp_gls = new double[bufSize];

  //getFile size
  while( temp = gzread(fp,tmp_gls,sizeof(double)*bufSize) ){
    nSites += temp/sizeof(double); 
  }

  gzclose(fp);
  delete[] tmp_gls;


  data = new double[nSites];
  fp=gzopen(fname,"r");
  size_t nSites2 = gzread(fp,data,sizeof(double)*nSites);
  gzclose(fp);
  nSites2 = nSites2/sizeof(double);

  if(nSites!=nSites2){
    fprintf(stdout,"reading failed nSites%lu != nSites2 %lu\n",nSites,nSites2);
  }
fprintf(stdout, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  

  return nSites2;
}



/*
double read-in
first get size (first read)
then read in data

Works for large files;

*/
size_t readGZbin2(const char* fname,double * &data){
  time_t t2=time(NULL);

  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  int temp=0;
  size_t nSites=0;
  double *tmp_gls = new double[bufSize];

  //getFile size
  while( temp = gzread(fp,tmp_gls,sizeof(double)*bufSize) ){
    nSites += temp/sizeof(double); 
  }

  gzclose(fp);
  delete[] tmp_gls;


  data = new double[nSites];
  fp=gzopen(fname,"r");
  size_t nSites2=0;
  size_t offset=0;
  int tmpBufSize=bufSize;
  while( temp = gzread(fp,data+offset,sizeof(double)*tmpBufSize) ){
    nSites2 += temp/sizeof(double); 
    offset+=tmpBufSize;
  }

  gzclose(fp);

  if(nSites!=nSites2){
    fprintf(stdout,"reading failed nSites%lu != nSites2 %lu\n",nSites,nSites2);
  }
fprintf(stdout, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  

  return nSites2;
}



/*
single read-in
realloc every chunk
maximum ram waste: bufSize doubles

Works for large files;

*/
size_t readGZbin3(const char* fname,double * &data){
  time_t t2=time(NULL);

  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  int temp=0;
  size_t nSites=0;


  data = (double *) malloc(sizeof(double)*bufSize);
  while( temp = gzread(fp,data+nSites,sizeof(double)*bufSize) ){
    nSites += temp/sizeof(double); 
    if(temp/sizeof(double)==bufSize)
      data = (double *) realloc(data,sizeof(double)*(nSites+bufSize));
  }

  gzclose(fp);



 fprintf(stdout, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  return nSites;
}

/*
single read-in
realloc doubles size
maximum ram waste: size of file

Works for large files;

*/
size_t readGZbin4(const char* fname,double * &data){
  time_t t2=time(NULL);

  size_t reBufSize=bufSize;
  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  int temp=0;
  size_t nSites=0;

  
  data = (double *) malloc(sizeof(double)*bufSize);
  while( temp = gzread(fp,data+nSites,sizeof(double)*bufSize) ){
    nSites += temp/sizeof(double); 
    if(nSites>=reBufSize){ // double the size
      reBufSize *=2;
      data = (double *) realloc(data,sizeof(double)*(reBufSize));
    }
  }

  gzclose(fp);



 fprintf(stdout, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  return nSites;
}

int main(int argc, char **argv){

  const char* gzBinFileName = NULL;
 bufSize=100000;

  argv++;
  while(*argv){
    if(strcmp(*argv,"-gzFile")==0 || strcmp(*argv,"-f")==0) gzBinFileName=*++argv; 
    else{
      fprintf(stderr,"Unknown arg:%s\n",*argv);
      return 0;
    }
    ++argv;
  }

  if(gzBinFileName==NULL){
      fprintf(stderr,"Please supply input fine: -f");
      return 0;
  }



  //   /*
  fprintf(stdout,"\nreadGZbin\t");
  double *data=NULL;
  size_t totalSites=readGZbin(gzBinFileName,data);
  fprintf(stdout,"read %lu sites\n\n",totalSites);
  //    */


  //  /*
  fprintf(stdout,"readGZbin2\t");
  double *data2=NULL;
  size_t totalSites2=readGZbin2(gzBinFileName,data2);
  fprintf(stdout,"read %lu sites\n\n",totalSites2);
  //    */

  //   /*
  fprintf(stdout,"readGZbin3\t");
  double *data3=NULL;
  size_t totalSites3=readGZbin3(gzBinFileName,data3);
  fprintf(stdout,"read %lu sites\n\n",totalSites3);
  //    */

  //   /*
  fprintf(stdout,"readGZbin4\t");
  double *data4=NULL;
  size_t totalSites4=readGZbin3(gzBinFileName,data4);
  fprintf(stdout,"read %lu sites\n\n",totalSites4);
  //    */



    
    for(size_t i=0;i<totalSites2; i++){
      if(data3[i]-data2[i] != 0 || data4[i]-data2[i] != 0){
//	fprintf(stdout,"diff\t%f\t%f\tsite: %lu\n",data3[i],data2[i],i);
	fprintf(stdout,"diff\t%f\t%f\t%f\t%f\tsite: %lu\t%lu\n",data[i],data2[i],data3[i],data4[i],i,totalSites2);
	exit(0);
      }

    }
 


 return 0;

}
