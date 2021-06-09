# preamble (options, installs, imports & custom functions) ----

options(warn=1) # really should be default in R

# specify the packages used:
required_packages = c(
	'github.com/rmcelreath/rethinking' # for rlkjcorr & rmvrnom2
	, 'github.com/mike-lawrence/aria/aria' # for Stan stuff
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

# check and compile the stan code ----
aria::enable_rstudio_syntax_compile()
# code_path = 'stan/pairwise_mvn_cp.stan'
code_path = 'stan/mvn_cp.stan'
file.edit(code_path)
#go enable "check syntax on save" and save to trigger check & compile

# sample ----
fs::dir_create('sampled')
out_path = fs::path(
	'sampled'
	, paste(
		fs::path_file(fs::path_ext_remove(code_path))
		, sim_pars$num_subj
		, sim_pars$num_vars
		, sim_pars$num_trials
		, sep = '_'
	)
	, ext = 'qs'
)
aria::compose(
	data = data_for_stan
	, code_path = code_path
	, out_path = out_path
	# , exe_args_list = list(
	# 	sample = list(
	# 		num_warmup = 2e3
	# 		, num_samples = 2e3
	# 	)
	# 	, init = .1
	# )
)

#examine posterior ----
fit = aria::coda(out_path)

#toss warmup
fit = filter(fit,!warmup)

# Check treedepth, divergences, & rebfmi
(
	fit
	%>% filter(!warmup)
	%>% group_by(chain)
	%>% summarise(
		max_treedepth = max(treedepth__)
		, num_divergent = sum(divergent__)
		, rebfmi = var(energy__)/(sum(diff(energy__)^2)/n())
	)
	# %>% summarise(
	# 	across(-chain,max)
	# )
)

# gather summary for core parameters (inc. rhat & ess)
(
	fit
	%>% filter(!warmup)
	%>% select(-(warmup:energy__))
	%>% select(!ends_with(fixed('_')))
	%>% select(!contains(fixed('_.')))
	%>% select(!contains(fixed('z_arr')))
	%>% select(!contains(fixed('z_mat')))
	# %>% select(!contains(fixed('z.')))
	%>% pivot_longer(
		cols = c(-chain,-iteration)
		, names_to = 'variable'
	)
	%>% group_by(variable)
	%>% summarise(
		rhat = posterior::rhat(matrix(value,ncol=length(unique(chain))))
		, ess_bulk = posterior::ess_bulk(matrix(value,ncol=length(unique(chain))))
		, ess_tail = posterior::ess_tail(matrix(value,ncol=length(unique(chain))))
		, as_tibble(t(posterior::quantile2(value,c(.1,.25,.5,.75,.9))))
		, .groups = 'drop'
	)
) ->
	fit_summary

# check the range of rhat & ess
(
	fit_summary
	%>% select(rhat,contains('ess'))
	%>% summary()
)

#means:
(
	fit_summary
	%>% filter(str_starts(variable,fixed('z_m.')))
	%>% left_join(
		(
			tibble(true=sim_pars$coef_means)
			%>% mutate(variable = paste('z_m',1:n(),sep='.'))
		)
	)
	%>% mutate(
		variable = factor(variable,levels=proper_level_order(variable))
	)
	%>% ggplot()
	+ geom_hline(yintercept = 0)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q10
			, ymax = q90
			, colour = ess_tail
		)
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q25
			, ymax = q75
			, colour = ess_bulk
		)
		, size = 2
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = q50
			, fill = rhat
		)
		, shape = 21
		, size = 2
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
	+ labs(
		y = 'Mean'
		, colour = 'ESS'
		, fill = 'Rhat'
	)
	+ theme(
		panel.background = element_rect(fill='grey50')
		, panel.grid = element_blank()
	)
)

# standard deviations:
(
	fit_summary
	%>% filter(str_starts(variable,'z_s'))
	%>% left_join(
		(
			tibble(true=sim_pars$coef_sds)
			%>% mutate(variable = paste('z_s',1:n(),sep='.'))
		)
	)
	%>% mutate(
		variable = factor(variable,levels=proper_level_order(variable))
	)
	%>% ggplot()
	+ geom_hline(yintercept = 1)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q10
			, ymax = q90
			, colour = ess_tail
		)
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q25
			, ymax = q75
			, colour = ess_bulk
		)
		, size = 2
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = q50
			, fill = rhat
		)
		, shape = 21
		, size = 2
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
	+ labs(
		y = 'SD'
		, colour = 'ESS'
		, fill = 'Rhat'
	)
	+ theme(
		panel.background = element_rect(fill='grey50')
		, panel.grid = element_blank()
	)
)

# correlations:
(
	fit_summary
	%>% filter(str_starts(variable,fixed('r.')))
	%>% left_join(
		(
			tibble(true=sim_pars$cor_mat[lower.tri(sim_pars$cor_mat)])
			%>% mutate(variable = paste('r',1:n(),sep='.'))
		)
	)
	%>% mutate(
		variable = factor(variable,levels=proper_level_order(variable))
	)
	%>% ggplot()
	+ geom_hline(yintercept = .9)
	+ geom_hline(yintercept = 0)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q10
			, ymax = q90
			, colour = ess_tail
		)
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q25
			, ymax = q75
			, colour = ess_bulk
		)
		, size = 2
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = q50
			, fill = rhat
		)
		, shape = 21
		, size = 2
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
	+ labs(
		y = 'Correlation'
		, colour = 'ESS'
		, fill = 'Rhat'
	)
	+ scale_y_continuous(limits=c(-1,1))
	+ theme(
		panel.background = element_rect(fill='grey50')
		, panel.grid = element_blank()
	)
)
