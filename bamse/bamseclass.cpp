using namespace std;

class Matrix {
private:
  int ind, variable, env, snp, hej, trait, mean_update_ratio,adj;	//numbers of INDividual,VARIABLE,ENViromental,SNPs and TRAIT data
  int thetrait;			//nummer trait som undersøges
  float (*matrix)[600];
  int (*miss_matrix)[600];
  float (*prior_miss_matrix)[600][3];
  int (*active_risk_parameters)[the_max_riskset];
  int (*p_active_risk_parameters)[the_max_riskset];
  int *fastbinom;
  int max_active;
  int include;
  int higher_mean;
  int the_trait;
  int max_riskset;
  float *trait_vektor;
  float *maf;
  int *miss;
  void calc_maf();
  void random_missing(int nr_snp,int nr_ind);
  int binary;
  //TRanrotWGenerator rg;   

public:
  //functioner
    Matrix (int thehej);	//constructor
  void make ( int antal_ind, int antal_env,
	     int antal_snp, int the_runtime,
	     int mean_update_ratio, int the_nr_mean_risk, int the_max_active, int include, int higher_meanint,int das_max_riskset,int the_adj,float the_min_range,float the_max_range,float the_range_frac,int the_binary);
  void new_snp (int riskset, int new_snp);
  void generate_max_min ();
  void generate_parameters (int exclude);
  void calc_p_likelihood ();
  bool generate_p_riskset ();
  void calc_mean (float das_mean[]);
  void calc_binom ();
  void update_mean_normal (int chosen_riskset);
  void update_varians (float varians_window);
  void update_snp_one ();
  void update_adjust (float adjust_window);
  void update_env_one (int chosen_riskset, int chosen_para);
  void update_active_risk_parameters ();
  void move_active_risk_parameters ();
  void gibbs_move_active_risk_parameters ();
  void accept ();
  void reject ();
  bool test (float prob_ratio);
  void print_riskset ();
  void print_risk_parameters ();
  void print_active_risk_parameters (int iteration);
  void print_sample(int iteration);
  void calc_average_mean ();
  void update_riskset_nr ();
  void create_riskset ();
  void delete_riskset ();
  void calc_prior_nr_riskset ();
  void calc_prior_nr_active ();
  float prior_nr_active_ratio(int A, int a);
  void calc_prior_nr_active_dist ();
  int true_nr_riskset ();
  void print_var ();
  int order(int nr_riskset);
  void delete_matrix();
  int check_env();
  void check_snp();
  void update_missing(int nr_snp);
  void update_missing_fP(int nr_snp);
  void perm_matrix();  
  void open_missingfile(float missing_cutoff);
  //variable
  float *max_env;
  float *min_env;
  float (*risk_parameters)[the_max_riskset];	//the parameters that defines the riskset  (SNP+env)
  float (*p_risk_parameters)[the_max_riskset];	//the proposed parameters that defines the riskset  (SNP+env)   
  float max_trait[1];
  float min_trait[1];
  float *mean_risk, *p_mean_risk, *sum_mean_risk, *mean_mean_risk;
  float likelihood, p_likelihood;
  float min_range, max_range,range_frac;
  int *p_riskset;
  int *riskset;
  int *nr_active;
  int *p_nr_active;
  int (*active)[the_max_riskset];
  int (*p_active)[the_max_riskset];
  float varians, p_varians;
  int nr_mean_risk, p_nr_mean_risk;
  float *adj_risk,*p_adj_risk;
  bool four_allowed;
  float *expected_mean;
  int *antal_in_riskset;
  float *riskset_prior;
  float *active_prior;
  float *active_dist_prior;
  float average_trait;
  float midrange;
  float kappa;
  int iteration;
};

//constructor
Matrix::Matrix (int four) {
    four_allowed = four;
}


void Matrix::make (int antal_ind,
		int antal_env, int antal_snp, int the_runtime,
		int the_mean_update_ratio, int the_nr_mean_risk,
		int the_max_active, int the_include, int the_higher_mean,int das_max_riskset,int the_adj,
    float the_min_range,float the_max_range,float the_range_frac,int the_binary)
{  
  max_riskset=das_max_riskset;   
  ind = antal_ind;
  env = antal_env;
  snp = antal_snp;
  adj=the_adj;  
  variable = antal_snp + antal_env + 1 + the_adj;


  trait = 1     ;
  the_trait=0;
  thetrait = the_trait;
  runtime = the_runtime;
  mean_update_ratio = the_mean_update_ratio;
  include=the_include;
  higher_mean=the_higher_mean;
  matrix = new float[ind][600];
  miss_matrix = new int[ind][600];
  prior_miss_matrix = new float[ind][600][3];
  adj_risk=new float[adj];
  p_adj_risk=new float[adj];
  //Psample = new float[int (runtime / sample_rate)][100][the_max_riskset];
  risk_parameters = new float[snp + env][the_max_riskset];
  p_risk_parameters = new float[snp + env][the_max_riskset];
  max_env = new float[env];
  min_env = new float[env];
  trait_vektor =new float[ind];
  riskset = new int[ind];
  p_riskset = new int[ind];
  maf = new float[snp];
  miss = new int[snp];
  mean_risk = new float[max_riskset];
  p_mean_risk = new float[max_riskset];
  sum_mean_risk = new float[max_riskset];
  mean_mean_risk = new float[max_riskset];
  nr_active = new int[the_max_riskset];
  p_nr_active = new int[the_max_riskset];
  active_risk_parameters = new int[snp + env][the_max_riskset];
  p_active_risk_parameters = new int[snp + env][the_max_riskset];
  max_active = the_max_active;
  fastbinom = new int[max_active+2];
  active = new int[max_active][the_max_riskset];
  p_active = new int[max_active][the_max_riskset];
  expected_mean = new float[max_riskset];
  antal_in_riskset = new int[max_riskset];
  riskset_prior = new float[max_riskset + 1];
  active_prior = new float[max_active + 1];
  active_dist_prior = new float[max_active + 1];
  char buffer[5000];
  min_range=the_min_range;
  max_range=the_max_range;
  range_frac=the_range_frac;
  binary=the_binary;
  ifstream exsamplefile (the_filename);
 if(!exsamplefile.is_open())
  err_stop("file is not open");
  int tal = 0;
  int done=0;
  while (!exsamplefile.eof ()) {
    exsamplefile.getline (buffer, 5000);
    char *buffer_ar[5000] = { 0 };
    for (int tal2 = 0; tal2 < variable+1; tal2++) {
      if (tal2)
	buffer_ar[tal2] = strtok (NULL, " ");
      else
	buffer_ar[tal2] = strtok (buffer, " ");

      if(tal2==variable){
        if(!buffer_ar[tal2])
         break;
	cout<<"line "<<tal+1<<" in the input file"<<endl;
	err_stop("The specified number of parameters does not match the number of colloms in the input file")
      }
      if(buffer_ar[tal2]==0&tal2<variable){
        done=1;
	break;
      }
      matrix[tal][tal2] = atof (buffer_ar[tal2]);
    }
    if(done)
      break;
    tal++;
  }
  if(tal!=ind){
   cout<<"specified number="<<ind<<" read numbers="<<tal<<endl;
    err_stop("the number of individuals does not match the number of lines in the file")

   }
  exsamplefile.close ();


  calc_binom();

  //check that the SNP parameters contain only 1 2 3
  //  if(check_env())
  // err_stop("Env error in make");
  generate_max_min ();
//calculates the riskset prior
  calc_prior_nr_riskset ();
//calculates the active prior
  calc_prior_nr_active ();
//calculates the active distribution prior
  calc_prior_nr_active_dist ();
  calc_average_mean ();

}


///////////////////////////////////////////////////////////////////////////////
//generate parameters
void Matrix::generate_parameters (int exclude = 0){
    calc_maf();
    for(int tal=0;tal<ind;tal++)
    trait_vektor[tal]=matrix[tal][thetrait];

  nr_mean_risk = 1;
  p_nr_mean_risk = 1;
  //generater de active parametre for det nye riskset (tal2)
  for (int tal11 = 1; tal11 < nr_mean_risk; tal11++) {
    if (max_active != 0) {
      for (int tal4 = 0; tal4 < snp + env; tal4++) {
    	p_active_risk_parameters[tal4][tal11] = 0;
	    active_risk_parameters[tal4][tal11] = 0;
      }
      int tal6 = 0;
      while (tal6 < rg.IRandom (1, max_active)) {
    	int tal5 = rg.IRandom (0, snp + env-1);
    	if (p_active_risk_parameters[tal5][tal11] == 0) {
    	  p_active_risk_parameters[tal5][tal11] = 1;
    	  active_risk_parameters[tal5][tal11] = 1;
    	  p_active[tal6][tal11] = tal5;
    	  active[tal6][tal11] = tal5;
    	  tal6++;
    	}
      }
      nr_active[tal11] = tal6;
      p_nr_active[tal11] = tal6;
    }
  }

//generates risk parameters
    for (int tal = 0; tal < snp; tal++) {	//random snp parameters
      for (int tal7 = 0; tal7 < nr_mean_risk; tal7++) {
        if(four_allowed)
  	      p_risk_parameters[tal][tal7] = risk_parameters[tal][tal7] = rg.IRandom (1, 7);	//all gets 4
        else {
          int temp = rg.IRandom (1, 6);
          if (temp > 3) 
            temp++;
          p_risk_parameters[tal][tal7] = risk_parameters[tal][tal7] = temp;
        }
      }
    }				//1=aa, 2=(Aa,aA), 3=AA, 4=(1,2,3), 5=(1,2), 6=(2,3), 7=(1,3)
    for (int tal2 = 0; tal2 < env; tal2++) {
      for (int tal6 = 1; tal6 < nr_mean_risk; tal6++) {
	    p_risk_parameters[tal2 + snp][tal6] = risk_parameters[tal2 + snp][tal6] = FRandom(min_env[tal2],max_env[tal2]);	//includes all
      }
    }
    if(check_env())
      err_stop("Error in generate parameters - env1");

    //adjust
  for(int tal=0;tal<adj;tal++){
    adj_risk[tal]=0;
    p_adj_risk[tal]=0;
  }
//generates mean_risks
  p_mean_risk[0] = mean_risk[0] = average_trait;
  for (int tal4 = 1; tal4 < nr_mean_risk; tal4++) {
    p_mean_risk[tal4] = mean_risk[tal4] = average_trait;
  }
//generates risksets:
  generate_p_riskset ();
  for (int nr_ind = 0; nr_ind < ind; nr_ind++) {
    riskset[nr_ind] = p_riskset[nr_ind];
    if (riskset[nr_ind] > nr_mean_risk) {
      cout << "fejl in generate" << endl;
      cin.get ();
    }
  }
// generates varians
  if (start_varians) {
    p_varians = varians = start_varians;
  }
  else {
    p_varians = varians = FRandom (0.5, 10);
  }
// calculates likelihood
  calc_p_likelihood ();
  likelihood = p_likelihood;

  if(check_env())
    err_stop("Error in generate parameters - env2");

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
  Matrix::generate_max_min ()
{
  for (int tal = 0; tal < env; tal++) {
    float
      themin = 3.4E+37;
    float
      themax = 0;
    for (int tal2 = 0; tal2 < ind; tal2++) {
      if (matrix[tal2][tal + trait + snp] < themin)
	themin = matrix[tal2][tal + trait + snp];
      if (matrix[tal2][tal + trait + snp] > themax)
	themax = matrix[tal2][tal + trait + snp];
    }
    max_env[tal] = themax;
    min_env[tal] = themin;
    //    cout << "#min env " << themin << " max env " << themax <<	" env nr " << tal << endl;
  }
  for (int tal = 0; tal < trait; tal++) {
    float
      themin = 3.4E+37;
    float
      themax = 0;

    for (int tal2 = 0; tal2 < ind; tal2++) {
      //cout<<themin<<" "<<themax<<endl;
      //cin.get();
      if (matrix[tal2][tal] < themin)
	themin = matrix[tal2][tal];
      if (matrix[tal2][tal] > themax)
	themax = matrix[tal2][tal];
    }
    max_trait[tal] = themax;
    min_trait[tal] = themin;
    if(min_range!=0|max_range!=0){
      max_trait[tal] = max_range;
      min_trait[tal] = min_range;                              
    }
    else if(range_frac!=1){
      float *trait_vec;
      trait_vec=new float[ind];
      for (int tal2 = 0; tal2 < ind; tal2++)
         trait_vec[tal2]=matrix[tal2][0]; 
      sort(trait_vec,trait_vec+ind);
      float maxx=(1-(1-range_frac)/2)*ind;
      float minn=(1-range_frac)/2*ind;
      max_trait[tal] = trait_vec[(int)maxx];
      min_trait[tal] = trait_vec[(int)minn];
      delete[] trait_vec;
    }
    if (print)
      cout << "# range of trait "<< min_trait[tal] << " - " << max_trait[tal]<< endl;
    //cin.get();
  }
  midrange =
    min_trait[the_trait] + (max_trait[thetrait] - min_trait[thetrait]) / 2.0;
  if (print)
    cout << "# midrange " << midrange << endl;
  kappa =
    1.0 * kappa_multi / Kvadrat (max_trait[thetrait] - min_trait[thetrait]);
  if (print)
    cout << "# kappa " << 1.0 / kappa << endl;
}

///////////////////////////////////////////////////////////////////////////////
//calculates the likelihood
void Matrix::calc_p_likelihood (){
/*
//////////////////////////////////////////    
  if(binary){  
    int i,q;
    double u,v,z_mean;
    for (int nr_ind=0;nr_ind<ind;nr_ind++){
      u=(((double)rand())/RAND_MAX);   	  
    	v=ndf(-p_mean_risk[p_riskset[nr_ind]]); 
    	if (matrix[nr_ind][0]==1.0){
    	  trait_vektor[nr_ind]=inv_ndf((v+u*(1.0-v)));	     	    
   	  }
    	else{
    	  trait_vektor[nr_ind]=inv_ndf((u*v));	     	    
  	  }
    } 
  }

///////////////////////////////////////////     
 */    
     
  float adjust=0;
  double  temp_likelihood = 0;
  if(!binary){
    if(adj>0){
      for (int nr_ind = 0; nr_ind < ind; nr_ind++) {
        adjust=0;
        for(int adj_nr=0;adj_nr<adj;adj_nr++)
          adjust=adjust+matrix[nr_ind][thetrait+snp+1+env]*p_adj_risk[adj_nr];
        temp_likelihood += Kvadrat (trait_vektor[nr_ind] - p_mean_risk[p_riskset[nr_ind]]-adjust);
      }
    }
    else{
      for (int nr_ind = 0; nr_ind < ind; nr_ind++) {
        temp_likelihood += Kvadrat (trait_vektor[nr_ind] - p_mean_risk[p_riskset[nr_ind]]-adjust);
      }
    }  
    p_likelihood =-log (2 * Pi) * ind / 2 - ind * log (p_varians) -    temp_likelihood / (2 * Kvadrat (p_varians));
  }
  else{
    for (int nr_ind=0;nr_ind<ind;nr_ind++)
    temp_likelihood +=trait_vektor[nr_ind]*log(p_mean_risk[p_riskset[nr_ind]])+(1-trait_vektor[nr_ind])*log(1-p_mean_risk[p_riskset[nr_ind]]);
    p_likelihood=temp_likelihood;
  }  
  if (p_likelihood < -100000000) {
    for (int tal = 0; tal < p_nr_mean_risk; tal++)
      cout << "mean" << tal << " " << p_mean_risk[tal] << endl;
    cout << " fejl i calc likelihood, p_like=" << p_likelihood <<
      " varians " << p_varians << " temp=" << temp_likelihood << endl;
    cin.get ();
  }
}



///////////////////////////////////////////////////////////////////////////////
void Matrix::update_varians (float varians_window = 0.2){
   p_varians = FRandom (varians - varians_window, varians + varians_window);
  if (varians < 0) {
    p_varians = -varians;
  }
  if (p_varians > the_max_varians)
    update_varians (0.2);
   bool boolean = test(0);
  if(boolean)
    varians = p_varians;
  else{
    p_varians = varians;
  }
}

///////////////////////////////////////////////////////////////////////////////
void  Matrix::update_snp_one (){

  if (nr_mean_risk == 1)
    return;
  int nr_riskset = rg.IRandom (1, nr_mean_risk - 1);	//select random riskset
  int * perm1;
  perm1 = new int[nr_active[nr_riskset]];
  if(nr_active[nr_riskset]==1){
     perm1[0]=1;
  }
  else{
    int new_seed = rg.IRandom (0, 1000000);
    perm_uniform (nr_active[nr_riskset], &new_seed,perm1);
  }
  int change_para = -1;
  bool any_snps = 0;
  for (int tal = 0; tal < nr_active[nr_riskset]; tal++) {
    change_para = active[perm1[tal] - 1][nr_riskset];
    if (change_para >= snp)
      continue;
    new_snp(nr_riskset,change_para);
    any_snps = 1;
    break;
  }
   delete[] perm1;
  if (!any_snps)
    return;
  generate_p_riskset ();
  test (0);
}

///////////////////////////////////////////////////////////////////////////////
void Matrix::update_env_one (int chosen_riskset = 0, int chosen_para = -1){
  if(!chosen_riskset){
    if(check_env())
      err_stop("error in env 1");
  }
  if (p_nr_mean_risk == 1)
    return;
  if (env == 0)
    return;
  int  nr_riskset;
  if (chosen_riskset)
    nr_riskset = chosen_riskset;
  else
    nr_riskset = rg.IRandom (1, p_nr_mean_risk - 1);
  if (chosen_para>-1) {
    int fortegn = 1;
    if (rg.IRandom (0, 1))
      fortegn *= -1;
    p_risk_parameters[chosen_para][nr_riskset] = FRandom (min_env[chosen_para - snp], max_env[chosen_para - snp]) * fortegn;
  }
  else {
    int new_seed = rg.IRandom (0, 1000000);
    int *perm;
    perm = new int[nr_active[nr_riskset]];
    perm_uniform (nr_active[nr_riskset], &new_seed,perm);
    
    for (int tal = 0; tal < nr_active[nr_riskset]; tal++) {
      int change_para = active[perm[tal] - 1][nr_riskset];
      if (change_para < snp)
      	continue;
      int fortegn = 1;
      if (rg.IRandom (0, 1))
      	fortegn *= -1;
       p_risk_parameters[change_para][nr_riskset] =FRandom (min_env[change_para - snp], max_env[change_para - snp]) * fortegn;
      generate_p_riskset ();
      test (0);
      break;

    }
    delete[] perm;
  }
  if(!chosen_riskset){
    if(check_env())
      err_stop("error in env 2");
  }
}

///////////////////////////////////////////////////////////////////////////////
void Matrix::update_active_risk_parameters (){
  if (nr_mean_risk == 1)
    return;
  for (int nr_ind = 0; nr_ind < ind; nr_ind++) {
    if (p_riskset[nr_ind] > nr_mean_risk) {
      cout << "fejl in update active " << endl;
      cin.get ();
    }
  }

  int error = 0;
  if (max_active == 1)
    return;
  int nr_riskset = rg.IRandom (1, nr_mean_risk - 1);	//the risk set
  /*   cout<<"før update active"<<endl;
     cout<<"nr active 1: "<<nr_active[1]<<" 2: "<<nr_active[2]<<" 3: "<<nr_active[3]<<endl;
     print_active_risk_parameters(0,0,-1);
     cout<<"riskset "<<nr_active[nr_riskset]<<endl;
     cin.get();
   */
  float prob = 1.0;			//pmM/pMm
  int random = rg.IRandom (0, 1);
  if (nr_active[nr_riskset] < 1|nr_active[nr_riskset] >max_active) {
    cout<<"fejl in active. current number active "<<nr_active[nr_riskset]<<endl;
    cin.get();
  }
  if (nr_active[nr_riskset] == max_active | nr_active[nr_riskset] >= snp + env) {

    prob = 1;
    if(random==1)
      return;
    random = 0;
  }
  if (nr_active[nr_riskset] == 1) {
    prob = 1; 
    if(random==0)
      return;
    random = 1;
  }
  generate_p_riskset ();
  if (random == 0) {		//deactive parameter
    while (1) {
      int tal2 = rg.IRandom (0, nr_active[nr_riskset] - 1);	//active parameter
      int old_active = active[tal2][nr_riskset];
      if (active_risk_parameters[old_active][nr_riskset] == 1) {
        p_active_risk_parameters[old_active][nr_riskset] = 0;
        p_nr_active[nr_riskset] = nr_active[nr_riskset] - 1;
	      for (int tal = tal2; tal < p_nr_active[nr_riskset]; tal++) {
	        p_active[tal][nr_riskset] = active[tal + 1][nr_riskset];
	      }
	      float  prior_active_ratio = prior_nr_active_ratio (nr_active[nr_riskset], p_nr_active[nr_riskset]);
	      float ratio =  log (prior_active_ratio * nr_active[nr_riskset] * 1.0 * prob /((env + snp - p_nr_active[nr_riskset])));
      	if (abs (ratio) < 0.000001)
      	  ratio=0;
//        cout<<ratio<<endl;
        if (abs (ratio) > 10000) {
      	  print_active_risk_parameters (0);
      	  cout << "error in deactivate ratio" << endl;
      	  cin.get ();
      	}
      	generate_p_riskset ();
      	bool  boolean = test (ratio);
      	if (boolean) {
      	  //  cout<<"deac"<<endl;
      	}
	      return;
      }
      error++;
      if (error > 1000 & error < 1002) {
      	print_active_risk_parameters (0);
      	cout << "stuck in update active 1" << endl;
      	cout << "old_active" << old_active << " nr active " << nr_active[nr_riskset] << endl;
      	cout << "riskset " << nr_riskset << " tal2 " << tal2 << endl;
      	cout << "active " << active[0][nr_riskset] << " " << active[1][nr_riskset] << " " << active[2][nr_riskset] << " "	  << active[3][nr_riskset] << " " << endl;
      }
    }
  }
  else {			//active parameter
    while (1) {
      int new_active = rg.IRandom (0, env + snp - 1);
      if (active_risk_parameters[new_active][nr_riskset] == 0) {
      	p_active_risk_parameters[new_active][nr_riskset] = 1;
      	p_nr_active[nr_riskset] = nr_active[nr_riskset] + 1;
      	if (new_active < snp) {
      	  if(four_allowed)
            p_risk_parameters[new_active][nr_riskset] = rg.IRandom (1, 7);
          else {
            int temp = rg.IRandom (1, 6);
            if (temp > 3) 
              temp++;
            p_risk_parameters[new_active][nr_riskset] = temp;
          }
      	}
      	else {
      	  update_env_one (nr_riskset, new_active);
      	}
      	p_active[p_nr_active[nr_riskset] - 1][nr_riskset] = new_active;
      	float prior_active_ratio = prior_nr_active_ratio (nr_active[nr_riskset],p_nr_active[nr_riskset]);
      	float ratio = log (prior_active_ratio * (env + snp - nr_active[nr_riskset]) * 1.0 * prob /(p_nr_active[nr_riskset]));
//          cout<<ratio<<endl;
        if (abs (ratio) > 10000) {
      	  cout << "error in activate ratio" << endl;
      	  cin.get ();
      	}
      	//cout<<"ratio "<<ratio<<endl;
      	generate_p_riskset ();
      	bool
      	  boolean = test (ratio);
      	if (boolean) {
      
      	}
      	return;
      }
      error++;
      if (error > 1000 & error < 1002)
      	cout << "stuck in update active 2" << endl;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
void
  Matrix::move_active_risk_parameters ()
{
  if (nr_mean_risk == 1)
    return;
  for (int nr_ind = 0; nr_ind < ind; nr_ind++) {
    if (p_riskset[nr_ind] > nr_mean_risk) {
      cout << "fejl in move" << endl;
      cin.get ();
    }
  }
  int nr_riskset = rg.IRandom (1, nr_mean_risk - 1);	//nr_riskset
  if (nr_active[nr_riskset] >= snp + env - 1)
    return;
  int tal2 = rg.IRandom (0, nr_active[nr_riskset] - 1);	//active parameter
  int old_active = active[tal2][nr_riskset];
  int error = 1;
  while (1) {
    int new_active = rg.IRandom (0, env + snp - 1);
    if (p_active_risk_parameters[new_active][nr_riskset] == 0) {
      p_active_risk_parameters[old_active][nr_riskset] = 0;
      p_active[tal2][nr_riskset] = new_active;
      p_active_risk_parameters[new_active][nr_riskset] = 1;
      if (new_active < snp){
         new_snp (nr_riskset,new_active); 
        // generate_p_riskset ();
      }
      else
      	update_env_one (nr_riskset, new_active);
         generate_p_riskset ();
      test (0);
      return;
    }
    error++;
    if (error > 1000 & error < 1002)
      cout << "stuck in move" << endl;
  }
}

///////////////////////////////////////////////////////////////////////////////
void Matrix::accept (){
  if(check_env()){
    cout<<"error in accept 1"<<endl;
  }
  for (int tal = 0; tal < snp + env; tal++) {
    for (int tal4 = 1; tal4 < p_nr_mean_risk; tal4++) {
      risk_parameters[tal][tal4] = p_risk_parameters[tal][tal4];
      active_risk_parameters[tal][tal4] = p_active_risk_parameters[tal][tal4];
    }
  }
  if(check_env()){
    cout<<"error in accept 12"<<endl;
  }
  for (int tal2 = 0; tal2 < ind; tal2++) {
    riskset[tal2] = p_riskset[tal2];
  }
  likelihood = p_likelihood;
  for (int tal3 = 0; tal3 < p_nr_mean_risk; tal3++) {
    mean_risk[tal3] = p_mean_risk[tal3];
  }
  for (int tal2 = 1; tal2 < p_nr_mean_risk; tal2++) {
    nr_active[tal2] = p_nr_active[tal2];
  }
         if(check_env()){
    cout<<"error in accept 14"<<endl;
  }
  for (int tal4 = 1; tal4 < p_nr_mean_risk; tal4++) {
    for (int tal = 0; tal < p_nr_active[tal4]; tal++) {
      active[tal][tal4] = p_active[tal][tal4];
    }
  }
           if(check_env()){
    cout<<"error in accept 15"<<endl;
  }
  nr_mean_risk = p_nr_mean_risk;
  if(check_env()){
    cout<<"error in accept 2"<<endl;
  }
}

///////////////////////////////////////////////////////////////////////////////
void  Matrix::reject (){
  for (int tal = 0; tal < snp + env; tal++) {
    for (int tal4 = 1; tal4 < nr_mean_risk; tal4++) {
      p_risk_parameters[tal][tal4] = risk_parameters[tal][tal4];
      p_active_risk_parameters[tal][tal4] = active_risk_parameters[tal][tal4];
    }
  }
  for (int nr_ind = 0; nr_ind < ind; nr_ind++) {
    if (riskset[nr_ind] > nr_mean_risk) {
      cout << "fejl i reject" << endl;
      cin.get ();
    }
  }
  for (int tal2 = 0; tal2 < ind; tal2++) {
    p_riskset[tal2] = riskset[tal2];
  }
  p_likelihood = likelihood;
  for (int tal3 = 0; tal3 < nr_mean_risk; tal3++) {
    p_mean_risk[tal3] = mean_risk[tal3];
  }
  for (int tal2 = 1; tal2 < nr_mean_risk; tal2++) {
    p_nr_active[tal2] = nr_active[tal2];
  }
  for (int tal4 = 1; tal4 < nr_mean_risk; tal4++) {
    for (int tal = 0; tal < nr_active[tal4]; tal++) {
      p_active[tal][tal4] = active[tal][tal4];
    }
  }
  p_nr_mean_risk = nr_mean_risk;
}

///////////////////////////////////////////////////////////////////////////////
bool Matrix::test (float prob_ratio = 0){
  if (abs (prob_ratio) > 10000) {
    cout << prob_ratio << endl;
    cout << "error in test ratio" << endl;
    cin.get ();
  }
  calc_p_likelihood ();
  if (prob_ratio + p_likelihood - likelihood >= 0) {
    accept ();
    return (1);
  }
  else if (prob_ratio + p_likelihood - likelihood < -20) {
    reject ();
    return (0);
  }
  else {
    float ratio = exp (prob_ratio + p_likelihood - likelihood);
    float random = rg.Random ();
     if (random < MIN (ratio, 1)) {
      accept ();
      return (1);
    }
    else {
      reject ();
      return (0);
    }
  }
  cout << "fejl i test" << endl;
  cin.get ();
}



///////////////////////////////////////////////////////////////////////////////
void
  Matrix::print_riskset ()
{
  for (int tal2 = 0; tal2 < ind; tal2++) {
    cout << tal2 << " riskset= " << riskset[tal2] << endl;
  }
}

///////////////////////////////////////////////////////////////////////////////
void
  Matrix::print_risk_parameters ()
{
  for (int tal2 = 0; tal2 < env + snp; tal2++) {
    for (int tal3 = 1; tal3 < nr_mean_risk; tal3++) {
      cout << tal2 << " risk_parameters= " << risk_parameters[tal2][tal3]
	<< endl;
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
void Matrix::print_sample(int iteration){
  int *total2;
  total2 = new int[nr_mean_risk];
  for (int tal = 0; tal < nr_mean_risk; tal++)
    total2[tal] = 0;
  for (int tal = 0; tal < ind; tal++)
    total2[riskset[tal]] += 1;
  
  cout <<"# "<< iteration << " like: " << likelihood << " SD: " << varians << " u0: " << mean_risk[0] << " " << " Tnr " << true_nr_riskset () <<" empty: "<<total2[0]<<" ";
  
  for (int tal2 = 0; tal2 < adj; tal2++)
    cout << adj_risk[tal2] << " ";
  cout<<endl;
  
  delete[] total2;

}

///////////////////////////////////////////////////////////////////////////////
void Matrix::print_active_risk_parameters (int iteraction)
{
int oneline = 0;
 int p =0;
 int nr_riskset = -1;
  float *das_mean;
  das_mean=new float[nr_mean_risk];
  calc_mean(das_mean);
  
  if (nr_mean_risk == 1) {
    return;
  }
  int *total2;
  total2 = new int[nr_mean_risk];
  if(nr_mean_risk-1<1)
    err_stop("tatal2 dumt");
  for (int tal = 0; tal < nr_mean_risk; tal++)
    total2[tal] = 0;
  for (int tal = 0; tal < ind; tal++)
    total2[riskset[tal]] += 1;


  if (oneline == 1) {
    for (int tal3 = 1; tal3 < nr_mean_risk; tal3++) {
   
      if (nr_riskset >= 0)
	tal3 = nr_riskset;
      for (int tal2 = 0; tal2 < env + snp; tal2++) {
	if (p)
	  cout << tal2 << " p_active_risk_parameters= " <<
	    p_active_risk_parameters[tal2][tal3] << endl;
	else
	  cout << tal2 << " active_risk_parameters= " <<
	    active_risk_parameters[tal2][tal3] << endl;
      }
      if (nr_riskset >= 0)
	break;
    }
  }
  else {
    int antal_risksets = nr_mean_risk;
    if (p)
      antal_risksets = p_nr_mean_risk;
    for (int tal3 = 1; tal3 < antal_risksets; tal3++) {
   cout<<"% "<<iteration<<" ";
      if (nr_riskset >= 0)
	tal3 = nr_riskset;
      for (int tal2 = 0; tal2 < env + snp; tal2++) {
	if (p) {
	  if (p_active_risk_parameters[tal2][tal3] == 1) {
	    cout << p_risk_parameters[tal2][tal3] << " ";
	  }
	  else {
	    cout << 0 << " ";
	  }
	}
	else {
	  if (active_risk_parameters[tal2][tal3] == 1) {
      float number=risk_parameters[tal2][tal3]*100;
      int numb=(int)number ;
      number=numb*1.0/100;
	    cout << number << " ";
	    if (tal2 <
		snp & (risk_parameters[tal2][tal3] <
		       1 | risk_parameters[tal2][tal3] > 7)) {
	      cout << "error in print" << endl;
	      cin.get ();
	    }
	  }
	  else {
	    cout << 0 << " ";
	  }
	}
      }

      
      cout << "nr " << total2[tal3] << " u" << tal3 << ": " <<
	mean_risk[tal3] <<" adjusted: "<<das_mean[tal3]<< endl;
      if (nr_riskset >= 0)
	break;
    }
  }
   delete[] total2;
   delete[] das_mean;
}


///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void Matrix::new_snp (int nr_riskset,int new_snp){
  if (four_allowed) {
    p_risk_parameters[new_snp][nr_riskset] = rg.IRandom (1, 7);
  }
  else {
    int temp = rg.IRandom (1, 6);
    if (temp > 3) 
      temp++;
    p_risk_parameters[new_snp][nr_riskset] = temp;
  }
}


////////////////////////////////////////////////////////////////
void  Matrix::calc_average_mean ()
{
  average_trait = 0;
  for (int nr_ind = 0; nr_ind < ind; nr_ind++) {
    average_trait += matrix[nr_ind][thetrait];
  }
  average_trait = average_trait / ind;
  //cout<<"average trait "<<average_trait<<endl;
}


//////////////////////////////////////////////////////////////
void Matrix::update_mean_normal (int chosen_riskset = -1){
  float log_ratio = 0;
  float the_prior_mean99 = 0;
  float var = 1.0 / kappa;
  int nr_riskset;
  int error = 0;

  if (chosen_riskset == -1) {
    nr_riskset = rg.IRandom (0, nr_mean_risk - 1);
    if (nr_mean_risk != p_nr_mean_risk) {
      cout << "error in update mean" << endl;
      cin.get ();
    }
  }
  else
    nr_riskset = chosen_riskset;

  if (nr_riskset) {
    if (use_midrange) {
      p_mean_risk[nr_riskset] = random_normal (midrange, var);
    }
    else {
      p_mean_risk[nr_riskset] = random_normal (average_trait, var);
    }
    bool correct_mean=1;
    while (correct_mean) {
      if (use_midrange) {
        p_mean_risk[nr_riskset] = random_normal (midrange, var);
      }
      else {
        p_mean_risk[nr_riskset] = random_normal (average_trait, var);
      }
      error++;
      if (error == 10 & chosen_riskset == -1) {	//mean proposal 2
        if (use_midrange) {
          the_prior_mean99 = midrange;
        }
	      else {
	        the_prior_mean99 = average_trait;
        }
        float min_others = p_mean_risk[0];
        float new_risk_mean = FRandom (min_others, max_trait[thetrait]);
        p_mean_risk[nr_riskset] = new_risk_mean;
        log_ratio =intNormal (mean_risk[nr_riskset], new_risk_mean,
        the_prior_mean99, var);
      }
      if (error == 10000902) {
      	cout << " stuck in update mean normal" << endl;
       	cout << "mean0 " << p_mean_risk[0] << endl;
        cout << "risk set nr " << nr_riskset << endl;
        cout << "chosen=" << chosen_riskset << endl;
      }
      correct_mean=0;
      if(higher_mean){
        if(p_mean_risk[nr_riskset] < min_trait[thetrait] | p_mean_risk[nr_riskset] > max_trait[thetrait] | p_mean_risk[nr_riskset] < p_mean_risk[0])
          correct_mean=1;  
      }
      else{
        if(p_mean_risk[nr_riskset] < min_trait[thetrait] | p_mean_risk[nr_riskset] > max_trait[thetrait])
          correct_mean=1;  
      }
    }
  }
  else if (nr_riskset == 0) {
    bool stop = 1;
    while (stop == 1) {
      int	stop = 0;
      if (use_midrange) {
        p_mean_risk[nr_riskset] = random_normal (midrange, var);
        the_prior_mean99 = midrange;
      }
      else {
        p_mean_risk[nr_riskset] = random_normal (average_trait, var);
        the_prior_mean99 = average_trait;
      }
      if(higher_mean)
        for (int tal = 1; tal < p_nr_mean_risk; tal++) {
          if (p_mean_risk[nr_riskset] > p_mean_risk[tal]) {
            stop = 1;
          }
        }
      if (p_mean_risk[nr_riskset] <  min_trait[thetrait] | p_mean_risk[nr_riskset] > max_trait[thetrait]) {
        stop = 1;
      }
      error++;
      if (error == 10 & chosen_riskset == -1) {	//mean proposal 2
        float min_others = max_trait[thetrait];
	      for (int tal = 1; tal < p_nr_mean_risk; tal++) {
	        if (min_others > p_mean_risk[tal]) {
	          min_others = p_mean_risk[tal];
          }
        }
        float new_risk_mean;
        if(higher_mean)
	        new_risk_mean = FRandom (min_trait[thetrait], min_others);
	      else
	        new_risk_mean = FRandom (min_trait[thetrait],max_trait[thetrait]);

        p_mean_risk[nr_riskset] = new_risk_mean;
        log_ratio = intNormal (mean_risk[nr_riskset], new_risk_mean,the_prior_mean99, var);
        stop = 0;
        break;
      }
      if (error == 11) {
        cout << " stuck in update mean normal2" << endl;
        for (int tal = 0; tal < p_nr_mean_risk; tal++) {
          cout << " mean" << tal << " " << p_mean_risk[tal];
        }
        cout << endl;
      }
      if (!stop) {
        break;
      }
    }
  }
  else
    cout << "error in update mean normal 3" << endl;

  if(higher_mean)
  if (p_mean_risk[nr_riskset] <
      min_trait[thetrait] | p_mean_risk[nr_riskset] >
      max_trait[thetrait] | p_mean_risk[nr_riskset] < p_mean_risk[0]) {
    cout << "wrong mean" << endl;
    cout << "p mean0 " << p_mean_risk[0] << "p mean" << nr_riskset << " " <<
      p_mean_risk[nr_riskset] << endl;
    cout << "  mean0 " << mean_risk[0] << "  mean" << nr_riskset << " " <<
      mean_risk[nr_riskset] << endl;
  }


  if (chosen_riskset == -1) {
    int b=0;
     b=order(nr_riskset);
    if(b&p_nr_mean_risk>2){
      generate_p_riskset ();
    }
    test (log_ratio);
  }

}


//////////////////////////////////////////////////////////////////////////////////////
void Matrix::update_riskset_nr (){
  if (rg.IRandom (0, 1)) {
    //   cout<<"create"<<endl;
    create_riskset ();
  }
  else {
    //    cout<<"delete"<<endl;
    delete_riskset ();
  }
}

/////////////////////////////////////////////////////////////////////////////
void Matrix::delete_riskset ()
{
//cout<<"delete"<<endl;
  if (nr_mean_risk == min_riskset)
    return;
  p_nr_mean_risk = nr_mean_risk - 1;

  int deleted = rg.IRandom (1, nr_mean_risk - 1);

  int total1 = 0;
  for (int tal = 0; tal < ind; tal++)
    if (riskset[tal] == deleted)
      total1++;

  float pmM = 1.0;
  float pMm = 1.0;			//(nr_mean_risk-1);
  pmM = pmM / (max_active);
  pmM = pmM / fastbinom[nr_active[deleted]];

  /////////////////////////////////
  for (int tal5 = deleted; tal5 < nr_mean_risk - 1; tal5++) {	// move old risk parameters
    for (int tal6 = 0; tal6 < snp + env; tal6++) {
      p_active_risk_parameters[tal6][tal5] = p_active_risk_parameters[tal6][tal5 + 1];
      p_risk_parameters[tal6][tal5] = p_risk_parameters[tal6][tal5 + 1];
    }
    p_nr_active[tal5] = p_nr_active[tal5 + 1];
    for (int tal = 0; tal < p_nr_active[tal5]; tal++) {
      p_active[tal][tal5] = p_active[tal][tal5 + 1];
    }
    p_mean_risk[tal5] = p_mean_risk[tal5 + 1];
  }
  /////////////////////////////////
  float prior_nr_of_active = active_prior[nr_active[deleted]];
  float prior_active_dist = active_dist_prior[nr_active[deleted]];
  float prior_nr_riskset = riskset_prior[nr_mean_risk] / riskset_prior[p_nr_mean_risk];
  //cout<<"p_nr_active "<<prior_nr_of_active<<endl;
  //cout<<"delete: "<<"pMm "<<pMm<<" pmM "<<pmM<<" p_nr_riskset "<<prior_nr_riskset<<endl;
  //float log_ratio=log(pmM*prior_nr_of_active*1.0/(pMm))  ;
  float log_ratio = log (pmM * 1.0 /(pMm * prior_active_dist * prior_nr_of_active * prior_nr_riskset));

  generate_p_riskset ();
  //cout<<log_ratio<<endl;
  bool boolean = test (log_ratio);
  if (!boolean) {
    p_mean_risk[0] = mean_risk[0];
    //   cout<<"reject"<<endl;
    if (total1 == 0 & log_ratio >= 0) {
      cout << "del " << deleted << "nr " << total1 << " of " << 	nr_mean_risk << " ";
      cout << "like " << likelihood << endl;
      cin.get ();
    }
  }
  else {
    //  cout<<"deleted"<<endl;
    mean_risk[0] = p_mean_risk[0];
  }
  if(check_env())
    err_stop("env error in delete risk sets");

}

////////////////////////////////////////////////////////////////////////////
void Matrix::create_riskset (){
  //cout<<"create"<<varians<<endl;
  if(check_env())
    err_stop("Error in create - env First");

  if (nr_mean_risk == max_riskset)
    return;
  float pMm = 1.0;
  p_nr_mean_risk = nr_mean_risk + 1;
  float pmM = 1.0;			//(p_nr_mean_risk-1);

  int new_riskset = p_nr_mean_risk - 1;
  for (int tal = 0; tal < env + snp; tal++) {
    p_active_risk_parameters[tal][new_riskset] = 0;
  }

  p_nr_active[new_riskset] = rg.IRandom (1, max_active);
  pMm = pMm / (max_active);
  int * perm2;
    perm2=new int[snp+env];
  if(snp+env==1){
    perm2[0]=1;
  }
  else{
    int new_seed = rg.IRandom (0, 1000000);
    perm_uniform (snp + env, &new_seed,perm2);
  }
  pMm = pMm / fastbinom[p_nr_active[new_riskset]];
  if(check_env())
    err_stop("Error in create - env 11");

  for (int tal = 0; tal < p_nr_active[new_riskset]; tal++) {
    p_active[tal][new_riskset] = perm2[tal] - 1;
    p_active_risk_parameters[p_active[tal][new_riskset]][new_riskset] = 1;
   if (p_active[tal][new_riskset] < snp) {
      if(four_allowed)
        p_risk_parameters[p_active[tal][new_riskset]][new_riskset] =rg.IRandom (1, 7);
      else {
        int temp = rg.IRandom (1, 6);
        if (temp > 3) 
          temp++;
        p_risk_parameters[p_active[tal][new_riskset]][new_riskset] = temp;
      }  
    }
    else{
      update_env_one (new_riskset, p_active[tal][new_riskset]);
    }
  }


  float prior_nr_of_active = active_prior[p_nr_active[new_riskset]];
  float prior_nr_riskset = riskset_prior[p_nr_mean_risk] / riskset_prior[nr_mean_risk];
  float prior_active_dist = active_dist_prior[p_nr_active[new_riskset]];
  float log_ratio = log (pmM * prior_active_dist * prior_nr_riskset * prior_nr_of_active * 1.0 / (pMm));
  p_mean_risk[new_riskset] = max_trait[thetrait];
  update_mean_normal (new_riskset);
  if (p_mean_risk[new_riskset] <
      min_trait[thetrait] | p_mean_risk[new_riskset] >
      max_trait[thetrait] | (p_mean_risk[new_riskset] < p_mean_risk[0]&higher_mean)) {
    cout << "wrong mean in create før" << endl;
    cout << "p mean0 " << p_mean_risk[0] << " pmean" << new_riskset << " "
	 << p_mean_risk[new_riskset] << endl;
    cout << "  mean0 " << mean_risk[0] << "  mean" << new_riskset << " " <<
      mean_risk[new_riskset] << endl;
  }
  generate_p_riskset ();

  bool boolean = test (log_ratio);
  if(check_env()){
    cout<<"error in riskset 52"<<endl;
  }
  if (!boolean) {
    p_nr_mean_risk = nr_mean_risk;
    p_mean_risk[0] = mean_risk[0];
    p_mean_risk[new_riskset] = mean_risk[new_riskset];
  }
  else {
    mean_risk[new_riskset] = p_mean_risk[new_riskset];
    nr_mean_risk = p_nr_mean_risk;
  }
  if(check_env()){
    cout<<"error in riskset 6"<<endl;
  }
  if ((p_mean_risk[new_riskset] <min_trait[thetrait] | p_mean_risk[new_riskset] > max_trait[thetrait] | (p_mean_risk[new_riskset] < p_mean_risk[0]&higher_mean)) & boolean) {
    cout << "test " << boolean << endl;
    cout << "wrong mean in create" << endl;
    cout << "p mean0 " << p_mean_risk[0] << " pmean" << new_riskset << " "
	 << p_mean_risk[new_riskset] << endl;
    cout << "  mean0 " << mean_risk[0] << "  mean" << new_riskset << " " <<
      mean_risk[new_riskset] << endl;
  }
  delete[] perm2;
  if(check_env())
    err_stop("env error in create risk sets last");

}
//////////////////////////////////////////////////////
//prior number of active parameters
void Matrix::calc_prior_nr_riskset (){
  float sum_prior = 0;
  for (int tal = min_riskset; tal <= max_riskset; tal++) {
    sum_prior = sum_prior + pow (1 - riskset_prior_p, tal);
  }
  for (int tal = min_riskset; tal <= max_riskset; tal++) {
    float prior = 1;
    prior = (pow (1 - riskset_prior_p, tal) * 1.0 / sum_prior);
    riskset_prior[tal] = prior;
  }
}
////////////////////////////////////////////////////////////////
void Matrix::calc_prior_nr_active_dist (){
  for (int tal = 1; tal <= max_active; tal++) {
    float prior = 0;
    prior = prior + binom (snp + env, tal);
    active_dist_prior[tal] = 1.0 / prior;
  }
}
///////////////////////////
void Matrix::calc_binom (){
  for (int tal = 0; tal <= max_active; tal++) {
    fastbinom[tal]=binom (snp + env,tal);
  }
}
///////////////////////////////////////////////////////////////
void Matrix::calc_prior_nr_active (){
  float sum_prior = 0;
  for (int tal = 1; tal <= max_active; tal++) {
    sum_prior = sum_prior + pow (1 - active_prior_p, tal);
  }
  for (int tal = 1; tal <= max_active; tal++) {
    float prior = 1;
    prior = (pow (1 - active_prior_p, tal) * 1.0 / sum_prior);
    active_prior[tal] = prior;
  }
}
///////////////////////////////////////////////////////////////
float Matrix::prior_nr_active_ratio(int A,int a){
  float prior;
  prior=(active_dist_prior[a]/active_dist_prior[A])*(active_prior[a]/active_prior[A]);
  return(prior);      
}
/////////////////////////////////////////////////
//calculated the true number of risksets (i.e. non empty)
/////////////////////////////////////////////////
int Matrix::true_nr_riskset (){
   int the_number = nr_mean_risk;
  int *total = new int[nr_mean_risk];
  if(total==NULL) exit(1);
  for (int tal = 0; tal < nr_mean_risk; tal++)
    total[tal] = 0;
  for (int tal = 0; tal < ind; tal++){
    if(riskset[tal]>nr_mean_risk-1){
      err_stop("fejl in true nr risk set");
    }
    else
      total[riskset[tal]] += 1;
  }
  for (int tal = 0; tal < nr_mean_risk; tal++)
    if (total[tal] == 0)
      the_number--;
    delete[] total;
   return (the_number);
}
//###################################
//prints the parameters
void
  Matrix::print_var ()
{
  cout << "# individuals " << ind << endl;
  cout << "# SNPs " << snp << endl;
  cout << "# environmental factors " << env << endl;
  cout << "# adjustment " << adj << endl;
  cout << "# maximum number of risk sets " << max_riskset << endl;
  cout << "# minimun number of risk sets " << min_riskset << endl;
  cout << "# maximum number of active parameters " << the_max_active << endl;
  cout << "# minimun number of active parameters " << the_min_active << endl;
  cout << "# prior for the number of active components " << active_prior_p << endl;
  cout << "# prior for the number of risk sets " << riskset_prior_p << endl;
  cout << "# the maximum sd " << the_max_varians << endl;
  cout << "# kappa multiplied by " << kappa_multi << endl;
  if (use_midrange)
    cout << "# midrange used as xi" << endl;
  else
    cout << "# the imperical mean used as xi" << endl;
  cout << "# the start varians " << start_varians << endl;
  cout << "# the burnin " << das_burnin << endl;
  cout << "# the run time " << runtime << endl;
  cout << "# the sample rate " << sample_rate << endl;
  cout << "# seed used " << theseed << endl;
  if(higher_mean)
    cout<<"# higher mean for risk sets"<<endl;
  else
    cout<<"# The mean are not assumed to be higher"<<endl;


}

//####################################################
//is the order of the risk sets changed?
int Matrix::order(int nr_riskset){
  if(higher_mean){
        if(nr_riskset==0)
        return(0);
    for(int tal=0;tal<p_nr_mean_risk;tal++){
      if(p_mean_risk[nr_riskset]<p_mean_risk[tal]){
        if(mean_risk[nr_riskset]>=p_mean_risk[tal]){
          //cout<<"ord"<<endl;
          return(1);
        }
      }
      else if(p_mean_risk[nr_riskset]>p_mean_risk[tal]){
        if(mean_risk[nr_riskset]<=p_mean_risk[tal]){
          //cout<<"ord"<<endl;
          return(1);        
        }                                        
      }    
            
    }
  }
  else{
    if(nr_riskset==0)
      return(1);
    for(int tal=0;tal<p_nr_mean_risk;tal++){
      if(abs(p_mean_risk[nr_riskset]-p_mean_risk[0])<abs(p_mean_risk[tal]-p_mean_risk[0])){
        if(abs(mean_risk[nr_riskset]-mean_risk[0])>=abs(p_mean_risk[tal]-p_mean_risk[0])){
          //cout<<"ord"<<endl;
          return(1);
        }
      }
      else if(abs(p_mean_risk[nr_riskset]-p_mean_risk[0])>abs(p_mean_risk[tal]-p_mean_risk[0])){
        if(abs(mean_risk[nr_riskset]-mean_risk[0])<=abs(p_mean_risk[tal]-p_mean_risk[0])){
          //cout<<"ord"<<endl;
          return(1);        
        }                                        
      }    
            
    }
  }
 return(0);   
}

///////////////////////////////////////////////////////////////////////////////
bool Matrix::generate_p_riskset (){
  if(p_nr_mean_risk==0)
    return(0);
  int nonempty = 0;
  for (int tal4 = 0; tal4 < ind; tal4++) {
    p_riskset[tal4] = 0;
  }
  int *order;
  order=new int[p_nr_mean_risk];
  if(higher_mean)
  OrderFloat(p_mean_risk,p_nr_mean_risk,order);
  else{
    float *temp;
    temp = new float[p_nr_mean_risk];
    for(int tal=0;tal<p_nr_mean_risk;tal++)  
      temp[tal]=abs(p_mean_risk[tal]-p_mean_risk[0]);
    OrderFloat(temp,p_nr_mean_risk,order);    
    delete[] temp;
  }
  //print_array(p_nr_mean_risk,order);
  //print_array_float(p_nr_mean_risk,p_mean_risk);

  for (int nr_riskset = 0; nr_riskset < p_nr_mean_risk-1; nr_riskset++) {	//nr_riskset is the riskset number
    for (int nr_ind = 0; nr_ind < ind; nr_ind++) {	//nr_ind is the ind number
      if(p_riskset[nr_ind])  //individuals is already in a riskset
        continue;
      for (int nr_snp = 0; nr_snp < p_nr_active[order[nr_riskset]]; nr_snp++) {	//nr_snp is the snp number
        if (p_active_risk_parameters[p_active[nr_snp][order[nr_riskset]]][order[nr_riskset]] == 0|p_active[nr_snp][order[nr_riskset]]>=snp) {
	        continue;		//only use the active parameters
	      }
//      	if (matrix[nr_ind][p_active[nr_snp][order[nr_riskset]] + trait] == 0&!include)	//missing data excludes
//          goto label;
//        if (matrix[nr_ind][p_active[nr_snp][order[nr_riskset]] + trait] != 0)	//missing data inkludes
      	  switch ((int) p_risk_parameters[p_active[nr_snp][order[nr_riskset]]][order[nr_riskset]]) {
          case 0:
            err_stop("missing error")
            break;
      	  case 1:
      	    if (matrix[nr_ind][p_active[nr_snp][order[nr_riskset]] + trait] != 1) {
              goto label;
      	    }
      	    break;
      	  case 2:
      	    if (matrix[nr_ind][p_active[nr_snp][order[nr_riskset]] + trait] != 2) {
      	      goto label;
      	    }
      	    break;
      	  case 3:
      	    if (matrix[nr_ind][p_active[nr_snp][order[nr_riskset]] + trait] != 3) {
      	      goto label;
      	    }
      	    break;
      	  case 4:
               if(!four_allowed){
                cout<<"error in generate risksets - four not allowed"<<endl;
                cin.get();
               }
      	    break;
      	  case 5:
      	    if (matrix[nr_ind][p_active[nr_snp][order[nr_riskset]] + trait] == 3) {
      	      goto label;
      	    }
      	    break;
      	  case 6:
      	    if (matrix[nr_ind][p_active[nr_snp][order[nr_riskset]] + trait] == 1) {
      	      goto label;
      	    }
      	    break;
      	  case 7:
      	    if (matrix[nr_ind][p_active[nr_snp][order[nr_riskset]] + trait] == 2) {
      	      goto label;
      	    }
      	    break;
      	  default:
      	    cout << "error in generate risk groups " << (int)
      	      p_risk_parameters[p_active[nr_snp][order[nr_riskset]]][order[nr_riskset]] << endl;
      	    break;
      	  }
      }				//1=aa, 2=(Aa,aA), 3=AA, 4=(1,2,3), 5=(1,2), 6=(2,3), 7=(1,3)

      for (int nr_snp = 0; nr_snp < p_nr_active[order[nr_riskset]]; nr_snp++) { 
	      if (p_active_risk_parameters[p_active[nr_snp][order[nr_riskset]]][order[nr_riskset]] == 0|p_active[nr_snp][order[nr_riskset]]<snp) {
	        continue;		//only use the active parameters
  	    }
//  	    cout<<p_active[nr_snp][order[nr_riskset]]<<endl;
      	if (p_risk_parameters[p_active[nr_snp][order[nr_riskset]]][order[nr_riskset]] > 0) {
      	  if (matrix[nr_ind][p_active[nr_snp][order[nr_riskset]] + trait] < p_risk_parameters[p_active[nr_snp][order[nr_riskset]]][order[nr_riskset]]) {
      	    goto label;
      	    break;
      	  }
      	}
      	else {
      	  if (matrix[nr_ind][p_active[nr_snp][order[nr_riskset]] + trait] >
      	      -p_risk_parameters[p_active[nr_snp][order[nr_riskset]]][order[nr_riskset]]) {
      	    goto label;
      	    break;
      	  }
      	}
      }

//      if (p_riskset[nr_ind] == 0) {
	    p_riskset[nr_ind] = order[nr_riskset];
/*      }
      else if (p_mean_risk[p_riskset[nr_ind]] < p_mean_risk[order[nr_riskset]]&higher_mean) {
	p_riskset[nr_ind] = order[nr_riskset];
      }
       else if (abs(p_mean_risk[p_riskset[nr_ind]] -p_mean_risk[0])< abs(p_mean_risk[order[nr_riskset]]-p_mean_risk[0])&!higher_mean) {
	p_riskset[nr_ind] = order[nr_riskset];
      }
*/
      nonempty = 1;
      label:
      continue;
    }
  }
  delete[] order;

  if (nonempty == 0) {
    return (0);
  }				//returns 0 if empty riskset
  return (1);			//returns 1 if non-empty riskset
}



//////////////////////////////////////////
void Matrix::update_adjust (float adjust_window = 0.2){
     
  if(adj==0)
    return;
  int chose=rg.IRandom(0,adj-1);
  p_adj_risk[chose] = FRandom (adj_risk[chose] - adjust_window, adj_risk[chose] + adjust_window);
  bool boolean=test();

  if(!boolean)
   p_adj_risk[chose]=adj_risk[chose];
  else
   adj_risk[chose]=p_adj_risk[chose];
}

/////////////////////////////////
void Matrix::delete_matrix(){
     
  delete[] matrix;

  delete[] adj_risk;
  delete[] p_adj_risk;
  delete[] risk_parameters;
  delete[] p_risk_parameters;
delete[]   max_env;
delete[]   min_env;
delete[]   trait_vektor;
delete[]   riskset;
delete[]   p_riskset;
delete[]   mean_risk;
delete[]   p_mean_risk;
delete[]   sum_mean_risk;
delete[]   mean_mean_risk;
delete[]   nr_active;
delete[]   p_nr_active;
delete[]   active_risk_parameters;
delete[]   p_active_risk_parameters;
delete[]   active;
delete[]   p_active;
delete[]   expected_mean;
delete[]   antal_in_riskset;
delete[]   riskset_prior;
delete[]   active_prior;
delete[]   active_dist_prior;     
     
     
}
//////////////////////////////////////////////////////
void Matrix::calc_mean (float das_mean[]){
  float *number_in_riskset;
  number_in_riskset=new float[nr_mean_risk];
  float adjust=0;
  for (int tal = 0; tal < nr_mean_risk; tal++) {
    das_mean[tal]=0;
    number_in_riskset[tal]=0;
  }
  for (int tal = 0; tal < ind; tal++) {
    if(adj>0){
      adjust=0;
      for(int adj_nr=0;adj_nr<adj;adj_nr++)
        adjust=adjust+matrix[tal][thetrait+snp+1+env]*p_adj_risk[adj_nr];
    }
    das_mean[riskset[tal]]+=p_mean_risk[p_riskset[tal]]+adjust;
    number_in_riskset[riskset[tal]]++;
  }
  for (int tal = 0; tal < nr_mean_risk; tal++)
    if(number_in_riskset[tal])
      das_mean[tal]=das_mean[tal]/number_in_riskset[tal];
    
  delete[] number_in_riskset;  
}
////////////////////////////////////////////////////////
int Matrix::check_env(){
  for(int env_nr=0;env_nr<env;env_nr++){
    for(int riskset_nr=1;riskset_nr<nr_mean_risk;riskset_nr++){
     if(active_risk_parameters[snp+env_nr][riskset_nr])
        if(abs(risk_parameters[snp+env_nr][riskset_nr])<min_env[env_nr]){
          cout<<"ENV"<<env_nr<<"riskset"<<riskset_nr<<"= "<<risk_parameters[snp+env_nr][riskset_nr]<<" minE "<<min_env[env_nr]<<" maxE "<<max_env[env_nr]<<endl;
          return(1);
        }
    }
    for(int riskset_nr=1;riskset_nr<p_nr_mean_risk;riskset_nr++){
     if(p_active_risk_parameters[snp+env_nr][riskset_nr])
         if(abs(p_risk_parameters[snp+env_nr][riskset_nr])<min_env[env_nr]){
          cout<<"pENV"<<env_nr<<"riskset"<<riskset_nr<<"= "<<p_risk_parameters[snp+env_nr][riskset_nr]<<" minE "<<min_env[env_nr]<<" maxE "<<max_env[env_nr]<<endl;
          return(1);
        }
    }
  }
  return(0);
}

/////////////////////////////////////////////////////
void Matrix::check_snp(){
  for(int nr_snp=0;nr_snp<snp;nr_snp++)
    for(int nr_ind=0;nr_ind<ind;nr_ind++)
      if(matrix[nr_ind][nr_snp+1]!=0&matrix[nr_ind][nr_snp+1]!=1&matrix[nr_ind][nr_snp+1]!=2&matrix[nr_ind][nr_snp+1]!=3)
	err_stop("SNP parameters must consist of 1,2,3 only!");
}

////////////////////////////////////////////////////
void Matrix::calc_maf(){
  for(int nr_snp=0;nr_snp<snp;nr_snp++)
    for(int nr_ind=0;nr_ind<ind;nr_ind++)
      miss_matrix[nr_ind][nr_snp]=0;

  for(int nr_snp=0;nr_snp<snp;nr_snp++){
    maf[nr_snp]=0;
    miss[nr_snp]=0;
    for(int nr_ind=0;nr_ind<ind;nr_ind++){
      maf[nr_snp]+=matrix[nr_ind][nr_snp+1]-1;
      if(matrix[nr_ind][nr_snp+1]==0){
        miss_matrix[miss[nr_snp]][nr_snp]=nr_ind;
        miss[nr_snp]++;
      }
    }
  maf[nr_snp]=maf[nr_snp]/((ind-miss[nr_snp])*2);
   // cout<<miss[nr_snp]<<" "<<maf[nr_snp]<<endl;
  }
  for(int nr_snp=0;nr_snp<snp;nr_snp++)
    for(int nr_miss=0;nr_miss<miss[nr_snp];nr_miss++)
      matrix[miss_matrix[nr_miss][nr_snp]][nr_snp+1]=(int)(FRandom(0,1)<maf[nr_snp])+(int)(FRandom(0,1)>maf[nr_snp]);
}
////////////////////////////////////////////////////
void Matrix::random_missing(int nr_ind=0,int nr_snp=0){
  if(nr_snp==0)
    nr_snp=rg.IRandom(0,(snp-1));  
  if(nr_ind==0)
    nr_ind=rg.IRandom(0,miss[nr_snp]-1);  
  matrix[nr_ind][nr_snp+1]=(int)(FRandom(0,1)<maf[nr_snp])+(int)(FRandom(0,1)>maf[nr_snp]);     
}

////////////////////////////////////////////////////
void Matrix::update_missing(int nr_snp=0){
  int *temp;
  int *p_temp;

  if(nr_snp==0)
    nr_snp=rg.IRandom(0,(snp-1));  
  if(miss[nr_snp]==0)
    return;
  temp =new int[miss[nr_snp]];
  p_temp =new int[miss[nr_snp]];
  for(int nr_miss=0;nr_miss<miss[nr_snp];nr_miss++){
    temp[nr_miss]=(int)matrix[miss_matrix[nr_miss][nr_snp]][nr_snp+1];
    p_temp[nr_miss]=(int)(FRandom(0,1)<maf[nr_snp])+(int)(FRandom(0,1)>maf[nr_snp]);
    matrix[miss_matrix[nr_miss][nr_snp]][nr_snp+1]=p_temp[nr_miss]; 
  }
  generate_p_riskset();
  if(!test()){
    for(int nr_miss=0;nr_miss<miss[nr_snp];nr_miss++)
      matrix[miss_matrix[nr_miss][nr_snp]][nr_snp+1]=temp[nr_miss];
  }
  delete[] temp;
  delete[] p_temp;
}
//////////////////////////////////////////////////
void Matrix::perm_matrix(){
    int new_seed = rg.IRandom (0, 1000000);
    int *perm;
    float *temp;
    perm = new int[ind];
    temp = new float[ind];
    perm_uniform (ind, &new_seed,perm);
    for(int nr_ind=0;nr_ind<ind;nr_ind++)
      temp[nr_ind]=matrix[nr_ind][0];
    for(int nr_ind=0;nr_ind<ind;nr_ind++)
      matrix[nr_ind][0]=temp[perm[nr_ind]-1];
    delete[] perm; 
}

////////////////////////////////////////////////
void Matrix::open_missingfile(float missing_cutoff){
  float (*temp)[6];
  temp=new float[snp*ind][6];
  char buffer[5000];
  ifstream exsamplefile (the_missingfile);
  if(!exsamplefile.is_open())
    err_stop("missingfile is not open");
  int tal = 0;
  int done=0;
  while (!exsamplefile.eof ()) {
    exsamplefile.getline (buffer, 5000);
    char *buffer_ar[5000] = { 0 };
    for (int tal2 = 0; tal2 < 5; tal2++) {
      if (tal2)
      	buffer_ar[tal2] = strtok (NULL, " ");
      else
	      buffer_ar[tal2] = strtok (buffer, " ");
      if(tal2==5){
        if(!buffer_ar[tal2])
         break;
	       cout<<"line "<<tal+1<<" in the input file"<<endl;
	       err_stop("in missingfile")
      }
      if(buffer_ar[tal2]==0&tal2<5){
        done=1;
	      break;
      }
      temp[tal][tal2] = atof (buffer_ar[tal2]);
    }
    if(done)
      break;
    tal++;
  }
  exsamplefile.close ();
  int n_missing=tal;
  for(int nr_mis=0;nr_mis<n_missing;nr_mis++){
    if(temp[nr_mis][2]>missing_cutoff){
      matrix[(int)temp[nr_mis][0]][(int)temp[nr_mis][1]+1]=1;
    }
    else if(temp[nr_mis][3]>missing_cutoff){
      matrix[(int)temp[nr_mis][0]][(int)temp[nr_mis][1]+1]=2;
    }
    else if(temp[nr_mis][4]>missing_cutoff){
      matrix[(int)temp[nr_mis][0]][(int)temp[nr_mis][1]+1]=3;
    }    
    else{
    prior_miss_matrix[(int)temp[nr_mis][0]][(int)temp[nr_mis][1]][0]=temp[nr_mis][2];
    prior_miss_matrix[(int)temp[nr_mis][0]][(int)temp[nr_mis][1]][1]=temp[nr_mis][3]+temp[nr_mis][2];
    prior_miss_matrix[(int)temp[nr_mis][0]][(int)temp[nr_mis][1]][2]=temp[nr_mis][4]+temp[nr_mis][3]+temp[nr_mis][2];
    }
  }
  delete[] temp;
}

////////////////////////////////////////////////////
void Matrix::update_missing_fP(int nr_snp=0){
  int *temp;
  int *p_temp;

  if(nr_snp==0)
    nr_snp=rg.IRandom(0,(snp-1));  
  if(miss[nr_snp]==0)
    return;
  temp =new int[miss[nr_snp]];
  p_temp =new int[miss[nr_snp]];
  for(int nr_miss=0;nr_miss<miss[nr_snp];nr_miss++){
    temp[nr_miss]=(int)matrix[miss_matrix[nr_miss][nr_snp]][nr_snp+1];
    float rand=FRandom(0,1);
    if(rand<prior_miss_matrix[miss_matrix[nr_miss][nr_snp]][nr_snp][0])
      p_temp[nr_miss]=1;
    else if(rand<prior_miss_matrix[miss_matrix[nr_miss][nr_snp]][nr_snp][1])
      p_temp[nr_miss]=2;
    else
      p_temp[nr_miss]=3;
                     
    matrix[miss_matrix[nr_miss][nr_snp]][nr_snp+1]=p_temp[nr_miss]; 
  }
  generate_p_riskset();
  if(!test()){
    for(int nr_miss=0;nr_miss<miss[nr_snp];nr_miss++)
      matrix[miss_matrix[nr_miss][nr_snp]][nr_snp+1]=temp[nr_miss];
  }
  delete[] temp;
  delete[] p_temp;
}

///////////////////////////////////////////////////////////////////////////////
void  Matrix::gibbs_move_active_risk_parameters ()
{
  return;
  if (nr_mean_risk == 1)
    return;
  for (int nr_ind = 0; nr_ind < ind; nr_ind++) {
    if (p_riskset[nr_ind] > nr_mean_risk) {
      cout << "fejl in move" << endl;
      cin.get ();
    }
  }
  int nr_riskset = rg.IRandom (1, nr_mean_risk - 1);	//nr_riskset
  if (nr_active[nr_riskset] >= snp + env - 1)
    return;
  int tal2 = rg.IRandom (0, nr_active[nr_riskset] - 1);	//active parameter

  int old_active = active[tal2][nr_riskset];
  if(old_active>=snp)
    return;

  p_active_risk_parameters[old_active][nr_riskset] = 0;
  int nr_inactive_snp=0;
  for(int nr_snp=0;nr_snp<snp;nr_snp++)
    if(p_active_risk_parameters[nr_snp][nr_riskset] == 0)
      nr_inactive_snp++;

  
  float *conditional_likelihood;
  conditional_likelihood = new float[nr_inactive_snp*6]; 
  int free_snp=0; 
  for(int nr_snp=0;nr_snp<snp;nr_snp++){
    if (p_active_risk_parameters[nr_snp][nr_riskset] != 0)
      continue;

    p_active[tal2][nr_riskset] = nr_snp;
    p_active_risk_parameters[nr_snp][nr_riskset] = 1;
    
    p_risk_parameters[nr_snp][nr_riskset]=1;
    generate_p_riskset ();
    calc_p_likelihood();
    conditional_likelihood[free_snp*6+0]=p_likelihood;

    p_risk_parameters[nr_snp][nr_riskset]=2;
    generate_p_riskset ();
    calc_p_likelihood();
    conditional_likelihood[free_snp*6+1]=p_likelihood;

    p_risk_parameters[nr_snp][nr_riskset]=3;
    generate_p_riskset ();
    calc_p_likelihood();
    conditional_likelihood[free_snp*6+2]=p_likelihood;

    p_risk_parameters[nr_snp][nr_riskset]=5;
    generate_p_riskset ();
    calc_p_likelihood();
    conditional_likelihood[free_snp*6+3]=p_likelihood;

    p_risk_parameters[nr_snp][nr_riskset]=6;
    generate_p_riskset ();
    calc_p_likelihood();
    conditional_likelihood[free_snp*6+4]=p_likelihood;

    p_risk_parameters[nr_snp][nr_riskset]=7;
    generate_p_riskset ();
    calc_p_likelihood();
    conditional_likelihood[free_snp*6+5]=p_likelihood;


    p_active_risk_parameters[nr_snp][nr_riskset] = 0;
    if(free_snp>nr_inactive_snp)
      err_stop("error in gibbs move ");
    free_snp++;
  }
  if(free_snp!=nr_inactive_snp){
    cout<<"free_snp="<<free_snp<<" inactive="<<nr_inactive_snp<<endl;
    err_stop("error in gibbs move2 ");
  }
  int choose=gibbs_pick(conditional_likelihood,6*nr_inactive_snp);
  int new_state=0;
  int free_num=0;
  for(int num=0;num<snp;num++){
    if (p_active_risk_parameters[num][nr_riskset] != 0)
      continue;

    if(free_num+1>choose*1.0/6){
      p_active_risk_parameters[num][nr_riskset] = 1;
      p_active[tal2][nr_riskset] = num;
      new_state=choose-free_num*6+1;
      if(new_state>3)
        new_state++;

      p_risk_parameters[num][nr_riskset] = new_state;
      break;
    }
    free_num++;
  }
  if(new_state==4)
    err_stop("error in gibbs move2.5")

  if(free_num>nr_inactive_snp)
    err_stop("error in gibbs move3")
      //                cout<<"snp="<<p_active[tal2][nr_riskset]<<" state="<<new_state<<" chose="<<choose<<endl;

  for(int tal=0;tal<nr_active[nr_riskset];tal++){
      if(p_active_risk_parameters[p_active[tal][nr_riskset]][nr_riskset] != 1)
	err_stop("error in gibbs move4");
  }
 
 
  int nr_inactive_snp2=0;
  for(int nr_snp=0;nr_snp<snp;nr_snp++)
    if(p_active_risk_parameters[nr_snp][nr_riskset] == 0)
      nr_inactive_snp2++;

  if(nr_inactive_snp!=nr_inactive_snp2+1){
    cout<<nr_inactive_snp<<" 2="<<nr_inactive_snp2<<endl;
    cout<<"after:";
    for(int tal=0;tal<nr_active[nr_riskset];tal++)
      cout<<p_active[tal][nr_riskset]<<". ";
      cout<<endl;
    cout<<"before:";
    for(int tal=0;tal<nr_active[nr_riskset];tal++)
      cout<<active[tal][nr_riskset]<<". ";
      cout<<endl;
    err_stop("error in gibbs move5")
    
  }
  generate_p_riskset ();
  calc_p_likelihood();
  accept();
  delete conditional_likelihood;

}
/*
//calculates the likelihood
void Matrix::calc_p_likelihood (){
     
     
  float adjust=0;
  double  temp_likelihood = 0;
  if(!binary){
    if(adj>0){
      for (int nr_ind = 0; nr_ind < ind; nr_ind++) {
        adjust=0;
        for(int adj_nr=0;adj_nr<adj;adj_nr++)
          adjust=adjust+matrix[nr_ind][thetrait+snp+1+env]*p_adj_risk[adj_nr];
        temp_likelihood += Kvadrat (trait_vektor[nr_ind] - p_mean_risk[p_riskset[nr_ind]]-adjust);
      }
    }
    else{
      for (int nr_ind = 0; nr_ind < ind; nr_ind++) {
        temp_likelihood += Kvadrat (trait_vektor[nr_ind] - p_mean_risk[p_riskset[nr_ind]]-adjust);
      }
    }  
    p_likelihood =-log (2 * Pi) * ind / 2 - ind * log (p_varians) -    temp_likelihood / (2 * Kvadrat (p_varians));
  }
  else{
    for (int nr_ind=0;nr_ind<ind;nr_ind++)
    temp_likelihood +=trait_vektor[nr_ind]*log(p_mean_risk[p_riskset[nr_ind]])+(1-trait_vektor[nr_ind])*log(1-p_mean_risk[p_riskset[nr_ind]]);
    p_likelihood=temp_likelihood;
  }  
  if (p_likelihood < -100000000) {
    for (int tal = 0; tal < p_nr_mean_risk; tal++)
      cout << "mean" << tal << " " << p_mean_risk[tal] << endl;
    cout << " fejl i calc likelihood, p_like=" << p_likelihood <<
      " varians " << p_varians << " temp=" << temp_likelihood << endl;
    cin.get ();
  }
}


*/
