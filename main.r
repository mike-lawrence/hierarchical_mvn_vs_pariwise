# preamble (options, installs, imports & custom functions) ----

options(warn=1) # really should be default in R

# specify the packages used:
required_packages = c(
	'github.com/rmcelreath/rethinking' # for rlkjcorr & rmvrnom2
	, 'github.com/stan-dev/cmdstanr' # for Stan stuff
	, 'github.com/stan-dev/posterior' # for posterior diagnostics & summaries
	, 'tidyverse' # for all that is good and holy
)

# load the helper functions:
for(file in fs::dir_ls('r_helpers')){
	cat('Loading function: ',fs::path_ext_remove(fs::path_file(file)),'()\n',sep='')
	source(file)
}
install_if_missing(required_packages)

# load the tidyverse (instead of imports bc lazy)
library(tidyverse)
library(magrittr) #for other pipes

# simulate data ----
set.seed(1) # change this to make different data

# setting the data simulation parameters
lst(
	# parameters you can play with
	num_subj = 3e1 # number of subjects, must be an integer >1
	, num_vars = 3 # number of 2-level variables manipulated as crossed and within each subject, must be an integer >0
	#num_trial will determine whether the cp or ncp models will sample better
	, num_trials = 1e3 # number of trials per subject/condition combo, must be an integer >1
	# the rest of these you shouldn't touch
	, num_coef = 2^(num_vars)
	, coef_means = rep(0,num_coef) #rnorm(num_coef)
	, coef_sds = rep(1,num_coef) #rweibull(num_coef,2,1)
	, cor_mat = diag(nrow=num_coef,ncol=num_coef) #rethinking::rlkjcorr(1,num_coef,eta=1)
	, noise = rweibull(1,2,1)
) -> sim_pars
sim_pars$cor_mat
summary(sim_pars$cor_mat[lower.tri(sim_pars$cor_mat)])
hist(sim_pars$cor_mat[lower.tri(sim_pars$cor_mat)],br=100,xlim=c(-1,1))

y = rethinking::rmvnorm2(
	n = 1e5 # yes, many more than we actually need; will eventually subset to sim_pars$num_subj
	, Mu = sim_pars$coef_means
	, sigma = sim_pars$coef_sds
	, Rho = sim_pars$cor_mat
)

#create odd/even as slight deviates from y
odd = y
dimnames(odd)=list(NULL,paste0('odd_v',1:ncol(odd)))
even = y
dimnames(odd)=list(NULL,paste0('even_v',1:ncol(even)))
for(i_k in 1:sim_pars$num_coef){
	odd[,i_k] = y[,i_k] + rnorm(nrow(y),0,sim_pars$coef_sds[i_k]/3)
	even[,i_k] = y[,i_k] + rnorm(nrow(y),0,sim_pars$coef_sds[i_k]/3)
}
y = cbind(odd,even)
dimnames(y) = list(
	NULL
	, c(
		paste0('odd_v',1:ncol(odd))
		, paste0('even_v',1:ncol(even))
	)
)

#now estimate the true values given this large-N set
sim_pars$num_coef = 2*sim_pars$num_coef
sim_pars$coef_means = apply(y,2,mean)
sim_pars$coef_sds = apply(y,2,sd)
sim_pars$cor_mat = cor(y)
summary(sim_pars$cor_mat[lower.tri(sim_pars$cor_mat)])
hist(sim_pars$cor_mat[lower.tri(sim_pars$cor_mat)],br=100,xlim=c(-1,1))

#now we trim down to the requested N
(
	y[1:sim_pars$num_subj,]
	%>% as_tibble()
	%>% mutate(
		subj = 1:n()
	)
) -> subj_coef


# compute the contrast matrix
(
	1:sim_pars$num_vars
	%>% map(.f=function(x){
		factor(c('lo','hi'))
	})
	%>% (function(x){
		names(x) = paste0('v',1:sim_pars$num_vars)
		return(x)
	})
	%>% cross_df()
	%>% get_contrast_matrix(
		formula = as.formula(paste('~1+',paste0('v',1:sim_pars$num_vars,collapse='*')))
		, contrast_kind = halfsum_contrasts
	)
) -> contrast_matrix

# get condition means implied by subject coefficients and contrast matrix
(
	subj_coef
	%>% pivot_longer(
		cols = -subj
		, names_to = c('split','name')
		, names_sep = '_'
	)
	%>% pivot_wider()
	%>% group_by(subj,split)
	%>% summarise(
		(function(x){
			out = attr(contrast_matrix,'data')
			out$cond_mean = as.vector(contrast_matrix %*% t(x))
			return(out)
		})(cur_data())
		, .groups = 'drop'
	)
) -> subj_cond

# get noisy measurements in each condition for each subject
(
	subj_cond
	%>% expand_grid(trial = 1:sim_pars$num_trials)
	%>% mutate(
		obs = rnorm(n(),cond_mean,sim_pars$noise)
	)
	%>% select(-cond_mean)
) -> dat

# dat is now what one would have if a real dataset had been collected
print(dat)

# Compute inputs to model ----

# W: the full trial-by-trial contrast matrix
(
	dat
	%>% mutate(dat_row=1:n())
	%>% group_by(subj,split)
	%>% summarise(
		dat_row = dat_row
		, as_tibble(get_contrast_matrix(
			data = cur_data()
			, as.formula(paste0('~1+',paste0('v',1:sim_pars$num_vars,collapse='*')))
			, contrast_kind = halfsum_contrasts
		))
		, .groups = 'drop'
	)
	%>% pivot_longer(
		cols = c(-split,-subj,-dat_row)
	)
	%>% unite(col='name',split,name)
	%>% pivot_wider(values_fill=0)
	%>% arrange(dat_row)
) -> W

# quick glimpse; lots of rows
glimpse(W)

# get the unique entries in W
(
	W
	%>% select(-subj,-dat_row)
	%>% distinct()
) -> uW
glimpse(uW)
# far fewer rows!

# for each unique condition specified by uW, the stan model will
# work out values for that condition for each subject, and we'll need to index
# into the resulting subject-by-condition matrix. So we need to create our own
# subject-by-condition matrix and get the indices of the observed data into a
# the array produced when that matrix is flattened.
(
	uW
	%>% mutate(
		uW_row = 1:n()
	)
	# first repeat the matrix so there's a copy for each subject
	%>% slice(
		rep(
			row_number()
			, each=length(unique(dat$subj))
		)
	)
	# add subject labels
	%>% mutate(
		subj = rep(sort(unique(dat$subj)),times=nrow(uW))
	)
	%>% arrange(uW_row,subj)
	# add a row index
	%>% mutate(
		uW_by_subj_row = 1:n()
	)
	# join to the full contrast matrix W
	%>%	right_join(
		W
		, by = c(names(uW),'subj')
	)
	# ensure ordered as the original data
	%>% arrange(dat_row)
	# pull the row label
	%>%	pull(uW_by_subj_row)

) -> obs_index

# package for stan & sample ----

data_for_stan = lst( # lst permits later entries to refer to earlier entries

	####
	# Entries we need to specify ourselves
	####

	# W: within predictor matrix
	w = as.matrix(uW)

	# num_subj: number of subjects
	, n = length(unique(dat$subj))

	# outcome: outcome on each trial
	, obs = dat$obs

	# obs_index: index of each trial in flattened version of subject-by-condition value matrix
	, obs_index = obs_index

	####
	# Entries computable from the above
	####

	# num_obs: total number of observations
	, o = length(obs)

	# num_rows_W: num rows in within predictor matrix W
	, kw = nrow(w)

	, is_intercept = as.numeric(guess_if_intercept(w))
)

# double-check:
glimpse(data_for_stan)

# sample the posterior ----
mod = cmdstanr::cmdstan_model(
	'stan/mvn_cp.stan'
	# 'stan/pairwise_mvn_cp.stan'
	, include = 'stan'
)
iter_warmup = 1e3
iter_sampling = 1e3
fit = mod$sample(
	data = data_for_stan
	, chains = parallel::detectCores()/2
	, parallel_chains = parallel::detectCores()/2
	, iter_warmup = iter_warmup
	, iter_sampling = iter_sampling
	, seed = 1
	, refresh = (iter_warmup+iter_sampling)/10
)

#save in case we want to look at the fit later
fit$save_object(
	paste(
		str_remove(fit$metadata()$model_name,'_model')
		, sim_pars$num_subj
		, sim_pars$num_vars
		, sim_pars$num_trials
		, 'fit.rds'
		, sep='_'
		, collapse='_'
	)
)
# see ?CmdStanMCMC for methods available for `fit`

#run cmdstan diagnostics for divergences, E-BFMI, etc
fit$cmdstan_diagnose()

# gather summary for core parameters (inc. rhat & ess)
(
	c('lp__','noise','z_m','z_s','r')
	%>% fit$draws()
	%>% posterior::summarise_draws(
		'rhat'
		, 'ess_bulk'
		, 'ess_tail'
		, 'median'
		, 'quantile2'
	)
	#adjust names a bit so proper_level_order works
	%>% mutate(
		variable = str_remove(variable,']')
		, variable = str_replace(variable,fixed('['),fixed('.'))
		, variable = str_replace(variable,',',fixed('.'))
	)
) -> fit_summary

# check the range of rhat & ess
(
	fit_summary
	%>% select(rhat,contains('ess'))
	%>% summary()
)

# show the full tibble
# print(fit_summary,n=nrow(fit_summary))

# plot the group-level parameters
(
	fit_summary
	%>% filter(str_starts(variable,'z_m'))
	%>% mutate(
		variable = factor(variable,levels=proper_level_order(variable))
	)
	%>% ggplot()
	+ geom_hline(yintercept = 0)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q5
			, ymax = q95
			, colour = ess_tail
		)
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = median
			, fill = rhat
		)
		, shape=21
		, size = 5
	)
	+ coord_flip()
	+ scale_color_gradient(
		high = 'white'
		, low = scales::muted('red')
	)
	+ scale_fill_gradient(
		low = 'white'
		, high = scales::muted('red')
	)
	+ labs(y = 'posterior')
	+ theme(
		panel.background = element_rect(fill='grey50')
		, panel.grid = element_blank()
	)
)

(
	fit_summary
	%>% filter(str_starts(variable,'z_s'))
	%>% mutate(
		variable = factor(variable,levels=proper_level_order(variable))
	)
	%>% ggplot()
	+ geom_hline(yintercept=1)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q5
			, ymax = q95
			, colour = ess_tail
		)
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = median
			, fill = rhat
		)
		, shape=21
		, size = 5
	)
	+ coord_flip()
	+ scale_color_gradient(
		high = 'white'
		, low = scales::muted('red')
	)
	+ scale_fill_gradient(
		low = 'white'
		, high = scales::muted('red')
	)
	+ labs(y = 'posterior')
	+ theme(
		panel.background = element_rect(fill='grey50')
		, panel.grid = element_blank()
	)
)

(
	fit$summary('r')
	%>% mutate(
		variable = factor(variable,levels=proper_level_order(variable))
	)
	%>% ggplot()
	+ geom_hline(yintercept = 0)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q5
			, ymax = q95
			, colour = ess_tail
		)
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = median
			, fill = rhat
		)
		, shape = 21
		, size = 5
	)
	+ coord_flip()
	+ scale_color_gradient(
		high = 'white'
		, low = scales::muted('red')
	)
	+ scale_fill_gradient(
		low = 'white'
		, high = scales::muted('red')
	)
	+ labs(y = 'posterior')
	+ scale_y_continuous(limits=c(-1,1))
	+ theme(
		panel.background = element_rect(fill='grey50')
		, panel.grid = element_blank()
	)
)

# plot the subject-level deviations-from-true
(
	fit$summary('z')
	%>% left_join(
		(
			subj_coef
			%>% (function(x){
				ncxm1 = ncol(x) - 1
				names(x)[names(x)!='subj'] = as.character(1:ncxm1)
				return(x)
			})
			%>% pivot_longer(
				-subj
				# , names_prefix = 'v'
			)
			%>% mutate(
				variable = paste0('z[',subj,',',name,']')
			)
		)
	)
	%>% mutate(
		median = median - value
		, q5 = q5 - value
		, q95 = q95 - value
		, name = factor(name,levels=proper_level_order(name))
	)
	%>% group_by(name)
	%>% arrange(median)
	%>% mutate(x=1:n())
	%>% ggplot()
	+ facet_grid(
		name~.
		, scales = 'free'
	)
	+ geom_hline(yintercept = 0)
	+ geom_linerange(
		mapping = aes(
			x = x
			, ymin = q5
			, ymax = q95
			, colour = median
		)
	)
	+ scale_color_gradient2(guide=F)
	+ theme(
		panel.background = element_rect(fill='grey50')
		, panel.grid = element_blank()
	)
)

