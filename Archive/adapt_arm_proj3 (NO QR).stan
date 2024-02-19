
 data {
  int<lower=1> N_obs; //number of observations
  int<lower=1> N_site; //number of sites
  int<lower=1> N_trt_groups; //number of treatment groups

  array[N_obs] int site; //site ids
  array[N_obs] int trt; //trt group
  array[N_obs] int resp; //outcome 

  real int_prior_mu; // intercept prior mu
  real int_prior_sd; // intercept prior sd
  array[N_trt_groups-1] real trt_prior_mu; // trt group priors mu
  array[N_trt_groups-1] real trt_prior_sd; // trt group prior SD

}
transformed data { 
    matrix[N_obs, N_trt_groups-1] X_matrix = rep_matrix(0, N_obs, N_trt_groups-1);
    matrix[N_obs, N_site] Z_matrix = rep_matrix(0, N_obs, N_site);
    for (i in 1:N_obs){
      for (j in 2:N_trt_groups){
        if (j==trt[i]) X_matrix[i,j-1] = 1;
      }  
      for(k in 1:N_site)    
      if (k==site[i]) Z_matrix[i,k] = 1;
    }

  matrix[N_obs, N_trt_groups-1] Q_ast;
  matrix[N_trt_groups-1,N_trt_groups-1] R_ast;
  matrix[N_trt_groups-1,N_trt_groups-1] R_ast_inverse;

  Q_ast = qr_thin_Q(X_matrix) * sqrt(N_obs - 1);
  R_ast = qr_thin_R(X_matrix) / sqrt(N_obs - 1);
  R_ast_inverse = inverse(R_ast);



}
parameters {
  real b0;
  vector[N_trt_groups-1] beta_trt_; //beta coefficients for each treatment group
  vector[N_site] alpha_site_raw; //random effect coefficients for site
  cholesky_factor_corr[N_site] lkj_corr; //random effect correlation 
}

transformed parameters {   
 vector[N_site] b0_site = b0 + 0.2*lkj_corr * alpha_site_raw;//implies: random intercept for site is sampled from multivariate normal with mean 0 and 100 SD
 
 
}
model {
  b0~normal(int_prior_mu,int_prior_sd);
  beta_trt_[1]~normal(trt_prior_mu[1],trt_prior_sd[1]);
  beta_trt_[2]~normal(trt_prior_mu[2],trt_prior_sd[2]);
  beta_trt_[3]~normal(trt_prior_mu[3],trt_prior_sd[3]);
  lkj_corr ~ lkj_corr_cholesky(5); 
  alpha_site_raw~std_normal();
  
  resp ~ bernoulli_logit( X_matrix * beta_trt_ + Z_matrix*b0_site);
}

generated quantities {
    
    vector[N_trt_groups] beta_trt;
    beta_trt[1]=0;
      beta_trt[2:N_trt_groups]=beta_trt_;

    vector[N_trt_groups] pred_prob_trt;
        pred_prob_trt = inv_logit(mean(b0_site) + beta_trt);

array[N_obs] int ypred;

ypred = bernoulli_logit_rng(mean(b0_site) + X_matrix * beta_trt[2:4]);
	
	int pp_trt2;
		if(max(beta_trt) == beta_trt[2]) pp_trt2=1;
		else pp_trt2=0;

	int pp_trt3;
		if(max(beta_trt) == beta_trt[3]) pp_trt3=1;
		else pp_trt3=0;
		
	int pp_trt4;
		if(max(beta_trt) == beta_trt[4]) pp_trt4=1;
		else pp_trt4=0;

  int ov_fut; //overall futility rule
    if(max(beta_trt) == beta_trt[1]) ov_fut = 1;
    else ov_fut = 0;

//prior predictive check;
array[N_trt_groups-1] real thetarep;

  thetarep[N_trt_groups-1] = normal_rng(trt_prior_mu[N_trt_groups-1],trt_prior_sd[N_trt_groups-1]);


}
