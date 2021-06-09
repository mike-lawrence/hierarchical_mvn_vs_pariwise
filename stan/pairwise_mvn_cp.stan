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

	// which_r: matrix-like array to facilitate extracting the pairwise correlations
	//    from the lower-tri-vector parameter representation
	array[kw,kw] int which_r ;
	int k = 0 ;
	for(i_kw in 1:(kw-1)){
		for(j_kw in (i_kw+1):kw){
			k = k + 1;
			which_r[i_kw,j_kw] = k ;
		}
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
	// row_vector[kw] z_[n] ;
	matrix[n,kw] z_ ;

	// r: population-level correlations (on cholesky factor scale) amongst within-subject predictors
	vector<lower=-1,upper=1>[num_r] r ;

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

	r ~ uniform(-1,1) ;

	////
	// Mid-hierarchy structure
	////

	// looping over all pairs of columns in z_
	for(i_kw in 1:(kw-1)){
		for(j_kw in (i_kw+1):kw){
			array[n] row_vector[2] pairz_ ;
			for(i_n in 1:n){
				pairz_[i_n] = [ z_[i_n,i_kw] , z_[i_n,j_kw] ] ;
			}
			pairz_ ~ multi_normal(
				[ z_m_[i_kw] , z_m_[j_kw] ]
				, quad_form_diag(
					to_matrix(
						[ 1 , r[which_r[i_kw,j_kw]]
						    , r[which_r[i_kw,j_kw]] , 1 ]
						,2,2
					)
					, [ z_s_[i_kw] , z_s_[j_kw] ]
				)
			) ;
		}
	}

	// Loop to compute the loc for each subject in each condition implied by the design matrix
	matrix[n,kw] m ;
	for(i_kw in 1:kw){
		m[,i_kw] = rows_dot_product( z_ , w_n[i_kw]	) ;
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
	matrix[n,kw] z = z_ * obs_sd ;

	// add the mean to any intercept coefficients:
	for(i_kw in 1:kw){
		if(is_intercept[i_kw]==1){
			z_m[i_kw] = z_m[i_kw] + obs_mean ;
			z[,i_kw] = z[,i_kw] + obs_mean ;
		}
	}


}
