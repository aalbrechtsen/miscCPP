#include "bamse.h"

int main (int argc, char *argv[]){
////////////////////////open and read options file///////////////////////////
  the_optionfile = "options.txt";	//default opotions file name
  the_missingfile = "missinggtprobs.txt";	//default opotions file name
  int nr_files = 1;		//the current comman argument

  if (!strcmp ("-o", argv[1])) {	//change the option file name
    the_optionfile = argv[2];
    nr_files = nr_files + 2;
  }
   if (!strcmp ("-m", argv[nr_files])) {	//change the option file name
    the_missingfile = argv[nr_files+1];
    nr_files = nr_files + 2;
  }

  ///open options file///
  float options[28];		//array of option arguuments
  ifstream exsamplefile (the_optionfile);	//opens option file
  if (!exsamplefile.is_open ())
    err_stop ("Error opening options file");
  int tal = 0;			// read option file
  char buffer[256];
  while (!exsamplefile.eof ()) {
    exsamplefile.getline (buffer, 1000);
    char *buffer_ar[1000] = { 0 };
    buffer_ar[0] = strtok (buffer, " ");
    options[tal] = atof (buffer_ar[0]);
    tal++;
    if (tal == 28)
      break;
  }
  exsamplefile.close ();
  int antal_ind = (int) options[0];
  int antal_snp = (int) options[1];
  int antal_env = (int) options[2];
  int nr_mean_risk = (int) options[3] + 1;	//start value (min 2)
  int das_max_riskset = (int) options[4] + 1;	//max_riskset
  min_riskset = (int) options[5] + 1;
  the_max_active = (int) options[6];
  the_min_active = (int) options[7];	//not implemented
  active_prior_p = options[8];	//the parameters for geometric distributed number of active parameteres
  riskset_prior_p = options[9];	//the parameters for geometric distributed number of risk sets
  the_max_varians = options[10];
  kappa_multi = options[11];
  use_midrange = options[12];	//if 1 used the midrange. else used the average mean
  print = (int) options[14];
  start_varians = options[13];	//sd
  sample_rate = (int) options[17];
  runtime = (int) options[15];
  das_burnin = (int) options[16];
  mean_update_ratio = (int) options[18];
  int binary = (int) options[19];
  int permutations=(int) options[20];
  int the_higher_mean=(int) options[21];
  int the_adj=(int) options[22];
  float the_min_range= options[23];
  float the_max_range= options[24];
  float the_range_frac= options[25];
  float missing_fP= (int)options[26];
  float missing_cutoff= options[27];
  float likelihood2;
///////////////////////////////////////////////////

  if (argc < 2) {
    cout << "please input filename(s)";
    return 0;
  }
  int error_check=0;
  float *perm_post;
  float *true_perm_post;
  perm_post = new float[permutations+1]; //array for the posterior value from each permutation
  true_perm_post = new float[permutations+1]; //array for the posterior value from each permutation
  
  while (nr_files < argc) {	//loops over the input files
    the_filename = argv[nr_files];
    nr_files++;
    cout << "# file: " << the_filename << endl;
    int the_four_allowed=0;
    int *true_number_of_risksets;
    true_number_of_risksets=new int[(int)(runtime/sample_rate)];
    int *number_of_risksets;
    number_of_risksets=new int[(int)(runtime/sample_rate)];
    Matrix *data = new Matrix (the_four_allowed);	//generates an matrix object
   data->make ( antal_ind, antal_env, antal_snp, runtime, mean_update_ratio, nr_mean_risk, the_max_active, permutations,the_higher_mean,das_max_riskset,the_adj,the_min_range,the_max_range,the_range_frac,binary);	//opens the file filename and initiates a matrix
   if(missing_fP)
     data->open_missingfile(missing_cutoff);
   data->generate_parameters();	//generates parametres (default=0 meaning random config. 1 means relaxed)
   data->print_var ();
   
   if (print)
     cout << "# start SD " << data->varians << endl;
   float like=data->likelihood;
   for(int nr_perm=0;nr_perm<=permutations;nr_perm++){ ///////////loop over permutations////////
     int tal3 = 1;
     if(nr_perm>0){
       data->perm_matrix();
       data->generate_parameters();	//generates parametres (default=0 meaning random config. 1 means relaxed)
     }
     for (int tal = 0; tal <= runtime; tal++) {//starts the MCMC algo
       
       if(data->check_env()){
         err_stop("env error 1");
       }
       ///////the algorithm/////////
       data->iteration = tal;
       if(error_check&das_burnin<tal)
	 cout<<"snp_one"<<tal<<" "<<data->varians<<" like="<<data->likelihood<<endl;
       data->update_snp_one ();	//with test
       if(error_check&das_burnin<tal)
	 cout<<"adjust"<<tal<<" "<<data->varians<<" like="<<data->likelihood<<endl;
       data->update_adjust ();	//with test
       if(error_check&das_burnin<tal){
	 cout<<"env"<<tal<<" "<<data->varians<<" like="<<data->likelihood<<endl;   
	 if(data->check_env())
	   err_stop("env error 2");
       }
       data->update_env_one ();	//with test
       if(error_check&das_burnin<tal)
	 cout<<"mean"<<tal<<" "<<data->varians<<" like="<<data->likelihood<<endl;    
       data->update_mean_normal ();//with test
       if(error_check&das_burnin<tal)
	 cout<<"move"<<tal<<" "<<data->varians<<" like="<<data->likelihood<<endl;      
       data->update_active_risk_parameters ();	//with test
       data->move_active_risk_parameters ();	//with test
       if(missing_fP)
	 data->update_missing_fP();//with test
       else 
	 data->update_missing();//with test
       //         data->gibbs_move_active_risk_parameters ();//with test
       
       if(error_check&das_burnin<tal)
	 cout<<"mean"<<tal<<" "<<data->varians<<" like="<<data->likelihood<<endl;
       data->update_mean_normal ();//with test
       if (!data->test ()&0&error_check) {
	 data->generate_p_riskset ();
	 data->calc_p_likelihood ();
	 likelihood2 = data->p_likelihood;
	 cout << tal << "fejl " << " p like " << likelihood2 << " like "<< data->likelihood << endl;
	 cin.get ();
       }
       if(error_check&das_burnin<tal)
	 cout<<"riskset"<<tal<<" "<<data->varians<<" like="<<data->likelihood<<endl;
       if(data->check_env()){
	 err_stop("env error 3");
       }
       data->update_riskset_nr ();	//with test
       if(data->check_env()){
	 err_stop("env error 5");
       }
       if(error_check&das_burnin<tal)
	 cout<<"mean/varian"<<tal<<" "<<data->varians<<" like="<<data->likelihood<<endl;
       for (int tal2 = 0; tal2 < mean_update_ratio; tal2++) {
	 data->update_varians (var_change);//with test
	 data->update_mean_normal ();
       }
       if(error_check&das_burnin<tal)
	 cout<<"print"<<tal<<" "<<data->varians<<" like="<<data->likelihood<<endl;
       if (tal == tal3*sample_rate) {     
	 tal3++;
	 if(error_check&das_burnin<tal)
	   cout<<"print2"<<tal<<" "<<data->varians<<" like="<<data->likelihood<<endl;
	 if (tal <= das_burnin)
	   continue;
	 
	 
	 true_number_of_risksets[tal3-2]=data->true_nr_riskset ();
	 number_of_risksets[tal3-2]=data->nr_mean_risk;
	 if (print&nr_perm==0) {
	   cout << endl;
	   data->print_sample(tal);
	   if(data->nr_mean_risk>1){
	     data->print_active_risk_parameters (tal);
	     like=data->likelihood;
	   }
	 }
	 if (print == 3)
	   data->print_active_risk_parameters (3);
       }
       for (int risktallet = 1; risktallet < data->nr_mean_risk; risktallet++) {
	 if (data->mean_risk[0] > data->mean_risk[risktallet]&the_higher_mean) {
	   cout << "fejl med mean" << endl;
	   cin.get ();
	 }
       }
     }				//end loop over updates
     float *post;
     post = new float[das_max_riskset];
     float *true_post;
     true_post=new float[das_max_riskset];
     int samp=(int)(runtime/sample_rate);
     for(int tal=0;tal<das_max_riskset;tal++){
       post[tal]=0;
       true_post[tal]=0;
     }
     
     for(int tal=(int)(das_burnin/sample_rate);tal<samp;tal++){
      post[number_of_risksets[tal]-1]+=1;
      true_post[true_number_of_risksets[tal]-1]+=1;
    }
     if(nr_perm==0)
       cout<<"# Posterior distribution with and without empty risk sets"<<endl;
     else
       cout<<"permutation "<<nr_perm<<" ";
     for(int tal=0;tal<das_max_riskset;tal++){
       post[tal]=post[tal]/(samp-(int)(das_burnin/sample_rate));
       true_post[tal]=true_post[tal]/(samp-(int)(das_burnin/sample_rate));
       if((true_post[tal]!=0|post[tal]!=0)&nr_perm==0)
	 cout<<"# "<<tal<<" riskset "<<post[tal]<<" "<<true_post[tal]<<endl;
    }
     perm_post[nr_perm]=1-post[0];
     true_perm_post[nr_perm]=1-true_post[0];
     cout<<"# P(nm>0)= "<<perm_post[nr_perm]<<" "<<true_perm_post[nr_perm]<<endl;    
     
   }// end loop over permutations
   data->delete_matrix();
   delete data;
  }// end loop over files
  if(permutations>0){
    cout<<"unpermutated posterior="<<perm_post[0]<<endl;
    int num=0;
    for(int tal=1;tal<=permutations;tal++)
      if(perm_post[tal]>perm_post[0])                 
        num++;
    float p_value=num*1.0/permutations;
    cout<<"simulated p-value="<<p_value<<endl; 
    num=0;
    for(int tal=1;tal<=permutations;tal++)
      if(true_perm_post[tal]>true_perm_post[0])                 
        num++;
    p_value=num*1.0/permutations;
    cout<<"simulated p-value="<<p_value<<" without empty risk sets"<<endl; 
    delete[] perm_post;
    delete[] true_perm_post;
  }
  
  return (0);
}
