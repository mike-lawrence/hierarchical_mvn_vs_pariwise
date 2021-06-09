//aria: compile=1
//aria: compile_debug=1
//aria: run_debug=0 #because the auto-generated debug data doesn't work
#include helper_functions.stan
data{

	// o: number of trials
	int<lower=1> o ;

	// obs: obs on each trial
	vector[o] obs ;

	// n: number of subj
	int<lower=1> n ;

	// kw: num of within coefficients
	int<lower=1> kw ;

	// w: unique entries in the within predictor matrix
	matrix[kw,kw] w ;

	// obs_index: index of each obs in flattened subject-by-condition value matrix
	int obs_index[o] ;

	// is_intercept: binary indicator of columns that reflect intercept parameters
	int<lower=0,upper=1> is_intercept[kw] ;

}
transformed data{

	// num_r: number of correlations implied by kw
	int num_r = (kw*(kw-1)) %/% 2 ;

	// obs_mean: mean of obs vector
	real obs_mean = mean(obs) ;

	// obs_sd: sd of obs vector
	real obs_sd = sd(obs) ;

	// obs_: observations scaled to have zero mean and unit variance
	vector[o] obs_ = (obs-obs_mean)/obs_sd ;

	//w_n: an array where each element is a row from w repeated for n rows
	matrix[n,kw] w_n[kw] ;
	for(i_kw in 1:kw){
		w_n[i_kw] = rep_matrix(w[i_kw],n) ;
	}

}
parameters{

	//for parameters below, trailing underscore denotes that they need to be un-scaled in generated quantities

	// noise_: observation-level measurement noise
	real<lower=0> noise_ ;

	// z_m_: mean (across subj) for each coefficient
	row_vector[kw] z_m_ ;

	// z_s_: sd (across subj) for each coefficient
	vector<lower=0>[kw] z_s_ ;

	// z_: latent values
	row_vector[kw] z_[n] ;

	// r_: population-level correlations (on cholesky factor scale) amongst within-subject predictors
	cholesky_factor_corr[kw] r_ ;

}
transformed parameters{
	real z_arr[n,kw] ;
	for(i_n in 1:n){
		for(i_kw in 1:kw){
			z_arr[i_n,i_kw] = z_[i_n,i_kw] ;
		}
	}
	matrix[n,kw] z_mat = to_matrix(z_arr) ;

}
model{

	////
	// Priors
	////

	// low-near-zero prior on measurement noise
	noise_ ~ weibull(2,1) ; // weibull(2,1) is peaked around .8

	// normal(0,1) priors on all means
	z_m_ ~ std_normal() ;

	// normal(0,1) priors on all sds
	z_s_ ~ weibull(2,1) ;

	// lkj prior on correlations
	r_ ~ lkj_corr_cholesky(1) ;


	////
	// Compute mid-hierarchy quantities
	////

	// compute z_ (non-centered parameterization)
	z_ ~ multi_normal_cholesky(
		z_m_
		, diag_pre_multiply(z_s_,r_)
	) ;

	// Loop to compute the loc for each subject in each condition implied by the design matrix
	matrix[n,kw] m ;
	for(i_kw in 1:kw){
		m[,i_kw] = rows_dot_product( z_mat , w_n[i_kw]	) ;
	}

	////
	// Compute likelihood
	////

	obs_ ~ normal( to_vector(m)[obs_index] , noise_ ) ;

}
generated quantities{

	// noise: *unscaled* measurement noise
	real noise = noise_ * obs_sd ;

	// z_s: *unscaled* sd (across subj) for each coefficient
	vector[kw] z_s = z_s_ * obs_sd ;

	// z_m: *unscaled* mean (across subj) for each coefficient
	row_vector[kw] z_m = z_m_ * obs_sd ;

	// z: *unscaled* coefficients for each subject
	matrix[n,kw] z = z_mat * obs_sd ;

	// add the mean to any intercept coefficients:
	for(i_kw in 1:kw){
		if(is_intercept[i_kw]==1){
			z_m[i_kw] = z_m[i_kw] + obs_mean ;
			z[,i_kw] = z[,i_kw] + obs_mean ;
		}
	}

	// r: the upper-tri of the correlation matrix, flattened to a vector for efficient storage
	vector[num_r] r = flatten_lower_tri( multiply_lower_tri_self_transpose(r_) ) ;

}
