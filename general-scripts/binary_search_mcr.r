#Notes - 
# Throughout this code we refer to e_orig as e0, and to e_switch as e1.
# The general inputs for our binary search are
	# 1) a function to minimize linear combinations of e0 & e1 (`get_h_min`) for the model class of interest
	# 2) loop-invariant parts of binary search (denoted by `sst`), passed to `get_h_min`


library(MASS)
library(Rcgmin)
library(dfoptim)
library(optimr)


# General function for cross-validation
	# fit_fun takes `y`, `X` as named arguments.
	# report_fun takes named inputs `y`,`X`,and `mod`, where `mod` is the output of fit_fun. report_fun should return a scalar.
CV_mod <- function(y,X, n_folds=min(15,length(y)), fit_fun, report_fun, fit_in_sample=FALSE){
	n_dat <- length(y)
	cv_breaks <- cut(1:n_dat,breaks=n_folds)
	for(i in 1:n_folds){
		cv_ind <- cv_breaks == levels(cv_breaks)[i]
		mod_i <- fit_fun(y=y[!cv_ind], X=X[!cv_ind,])
		CV_report_i <- report_fun(mod=mod_i, X = X[cv_ind,], y = y[cv_ind])
		if(i==1){
			in_sample_report <- 
			CV_report <- matrix(NA,n_folds, length(CV_report_i))
			colnames(in_sample_report) <-
			colnames(CV_report) <- names(CV_report_i)
		}
		if(fit_in_sample) in_sample_report[i,] <- report_fun(mod=mod_i, X = X[!cv_ind,], y = y[!cv_ind])
		CV_report[i,] <- CV_report_i
	}
	cv_out <- colMeans(CV_report)
	if(!fit_in_sample) return(cv_out)
	if(fit_in_sample) return(c(
		'in_sample'=colMeans(in_sample_report),
		'cv'  = cv_out
	))
	
}


#' Helper function for loop invariant calculations in binary search, for general model classes
#' Returns an augmented ("full") dataset, containing the original dataset, in addition to terms used to approximate or calculate e_switch.
#' @param p1 indeces for the variables to be switched
#' @param times setting `times=2` corresponds to e_divide to approximate e_switch. Increasing `times` further increases the number of terms used in the approximation of e_switch. If `times =n,` all permutations are returned.
get_full_sample <- function(y,X,p1,n=length(y),times=2){
	p <- dim(X)[2]
	p2 <- setdiff(1:p,p1)

	if(times<=1) stop('If times=1, then no permuted values are created.')
	if(times%%1 != 0 | times>n) stop('times must be an integer <= n.')

	# Steps:
	# 	1 - create several groups of indeces (n_groups of them)
			# note, times is also the size of each group
	# 	2 - in each group, compute all possible combinations
	#	3 - tag elements by whether they correspond to permutations or whether they are contained in the original sample `is_perm`.

	combn_inds <- c() # to hold indeces for permuting data.
	n_groups <- floor(n/times)
	if(n_groups != n/times) warning(paste(n%%times, 'observation(s) dropped before creating permuted data'))
	starts <- 1 + times*(0:(n_groups-1)) #beginning of each group
	for(j in 1:n_groups){
		# store all combinations within the group
		group_inds <- which(1:n >= starts[j] & 1:n < starts[j]+times)
		combn_inds <- rbind(combn_inds,expand.grid(group_inds,group_inds))
	}
	if(any(duplicated(combn_inds))) stop('dup error')
	# combn_inds
	colnames(combn_inds) <- c('p1','p2')
	

	full_n <- dim(combn_inds)[1]
	if(n>50 & times==n){
		if(full_n!=n^2) stop('n error') #workcheck
		warning('full sample is of size n^2 = ',full_n,'. Consider setting "times" smaller') 
	}

	pre_mix <- sample(n)
	full_X <- array(NA, dim=c( full_n, dim(X)[2]))
	full_X[,p1] <- (X[pre_mix,])[combn_inds[,1],p1]
	full_X[,p2] <- (X[pre_mix,])[combn_inds[,2],p2]

	return(list(
		is_perm = combn_inds[,1] != combn_inds[,2],
		full_X = full_X,
		full_y = (y[pre_mix])[combn_inds[,2]]
	))
}



#  __  __    _____   _____
# |  \/  |  / ____| |  __ \
# | \  / | | |      | |__) |
# | |\/| | | |      |  _  /
# | |  | | | |____  | | \ \
# |_|  |_|  \_____| |_|  \_\





#' @param s If searching for MCR+, set `s` to 1. If searching for MCR-, set `s` to -1. 
#' @param eps loss threshold on the absolute scale
getMCR_internal <- function(s, eps, get_h_min, tol_mcr = 2^(-20), force_lower_0=FALSE, maxiter_BS = 400, verbose=TRUE,...){
	#We search to set gam as low as possible while still meeting cond_move_down

	# Verbose output
	vcat <- function(...){ if(verbose){ cat(...) }}

	#Check feasibility
	e0_min_model <- get_h_min(s=1,gam=0,...)
	if(e0_min_model$e0_value > eps){
		stop('Best model does not meet eps constraint')
	}
	validMR <- e0_min_model$e1_value / e0_min_model$e0_value #baseline to compare against

	#### Conditions
	cond_stop <- function(model_gam){
		((model_gam$h_value == 0 & model_gam$e0_value <= eps) | 
		 (model_gam$h_value >= 0 & model_gam$e0_value == eps))
	}
	cond_move_down <- function(model_gam){
		 (model_gam$h_value  > 0 & model_gam$e0_value <  eps)
	}
	cond_move_up <- function(model_gam){
		 !( cond_move_down(model_gam) | cond_stop(model_gam) )
	}


	#### initialize
	gam_lower <- -1
	gam_upper <- 1
	if(s==1) gam_upper <- 0
	if(s==-1) gam_lower <- 0  #if force_lower_0, only do convex searches within binary search

	#### Expand binary search region
	max_expand <- 100
	mod_upper <- get_h_min(s=s, gam=gam_upper, ...)
	for(i in 1:max_expand){
		if(!cond_move_up(mod_upper)) break
		gam_upper <- 2*abs(gam_upper+1)
		mod_upper <- get_h_min(s=s, gam=gam_upper, ...)
	}
	mod_lower <- get_h_min(s=s, gam=gam_lower, ...)
	if(!force_lower_0 | s==1){
		for(i in 1:max_expand){
			if(!cond_move_down(mod_lower)) break
			gam_lower <- -2*abs(gam_lower-1)
			mod_lower <- get_h_min(s=s, gam=gam_lower, ...)
		}
	}

	gam_opt <- NA
	do_BS_search <- TRUE
	if(cond_move_down(mod_lower)){
		warning('lower limit reached, result may be conservative')
		gam_opt <- gam_lower
		model_gam_opt <- mod_lower
		do_BS_search <- FALSE
	}
	if(cond_move_up(mod_upper)){
		stop('upper limit reached - size of Rashomon set may be too small? Stopping search to avoid errors.')
	}
	if(gam_lower>gam_upper) stop('search endpoint error')

	#### quick optimality check
	if(cond_stop(mod_lower)){gam_opt <- gam_lower; model_gam_opt <- mod_lower}
	if(cond_stop(mod_upper)){gam_opt <- gam_upper; model_gam_opt <- mod_upper}

	#### Run binary search
	searchpath<-data.frame(
		'e0'=mod_upper$e0_value,
		'e1'=mod_upper$e1_value,
		'h'=mod_upper$h_value,
		'gam'=gam_upper)
	if(do_BS_search){
		vcat('\n\nStarting binary search\n')
		for(i in 1:maxiter_BS){
			gam_mid <- (gam_lower+gam_upper)/2	
			time_i <- system.time({
				mod_mid <- get_h_min(s=s, gam=gam_mid, mod1 =mod_lower, mod2=mod_upper, ...)
			})['elapsed']

			r5 <- function(x) round(x,digits=5)
			vcat('\n sign:',s,'; iter:', i,'; gam:', r5(gam_mid),'; converg:', r5(gam_upper-gam_lower),'; time:',time_i)

			searchpath<-rbind(searchpath,c(
				'e0'=mod_mid$e0_value,
				'e1'=mod_mid$e1_value,
				'h'=mod_mid$h_value,
				'gam'=gam_mid
			))

			if(cond_stop(mod_mid)){
				gam_upper <- gam_lower <- gam_mid
				mod_upper <- mod_lower <- mod_mid
				break
			}
			if(cond_move_down(mod_mid)){
				gam_upper <- gam_mid
				mod_upper <- mod_mid
				vcat('.. down')
			}
			if(cond_move_up(mod_mid)){
				gam_lower <- gam_mid
				mod_lower <- mod_mid
				vcat('.. up')
			}

			if(abs(gam_upper-gam_lower) < abs(tol_mcr*eps)) break
		}
		gam_opt <- gam_upper
		model_gam_opt <- mod_upper
	}


	#### Summarize Results
	vcat('\n\nFound gamma at ',gam_opt,'\n\n')

	if(s==1)  out <- c((model_gam_opt$h_value/eps - 1) * gam_opt^-1)
	if(s==-1) out <- c((model_gam_opt$h_value/eps) - gam_opt)
	if(model_gam_opt$e0_value <= eps){
		validMR <- model_gam_opt$e1_value / model_gam_opt$e0_value
	}else{
		stop('Error in e0 calculation; final value of search should satisfy cond_move_down. Rashomon set may be too small?')
	}
	if(out<0) warning('Negative error')


	#Approximation error compared to a nonlinear approach
	approximation_error_linear <- abs(out - validMR)

	#Approximation error due to binary search resolution, and multiple solutions for get_h_min due to flat parts of the search space.
	approximation_error_search <- NA
	if(mod_lower$h_value > 0 & mod_upper$h_value > 0 & do_BS_search){
		if(sum( c(mod_lower$e0_value, mod_upper$e0_value) > eps) != 1) stop('move condition error') #we should have one model with e0>eps, and one with e0<= eps

		# find the slope and intercept of the line L
		# connecting the coordinates of 
		# mod_upper & mod_lower in e0,e1 space,
		# and take it's intersection with the line e0=eps.
		slope <- (mod_lower$e1_value - mod_upper$e1_value)/(mod_lower$e0_value - mod_upper$e0_value)
		int <- mod_lower$e1_value - slope * mod_lower$e0_value
		worst_case_out <- (int + slope * eps)/eps
		approximation_error_search <- out - worst_case_out
	}


	return(list('mcr'=out,
			'approximation_error_linear'=approximation_error_linear,
			'approximation_error_search'=approximation_error_search,
			'gam_opt'=gam_opt,
			s=s, 
			full_model_output=model_gam_opt,
			searchpath=searchpath
		))

}

get_empirical_MCR <- function(...){
	minus <- getMCR_internal(s=-1,...)
	plus  <- getMCR_internal(s= 1,...)
	return(list(
		range=c(minus$mcr,plus$mcr),
		minus=minus,
		plus=plus
	))
}

# Plot search path from binary search
#' @param emp_mcr output from `get_empirical_MCR`
plot_empirical_MCR_searchpath <- function(emp_mcr, eps=NA, ...){
	par(mfrow=c(1,2))
	e0all <- c(emp_mcr$minus$searchpath[,'e0'],
			emp_mcr$plus$searchpath[,'e0'])
	e1all <- c(emp_mcr$minus$searchpath[,'e1'],
			emp_mcr$plus$searchpath[,'e1'])
	hall <- c(emp_mcr$minus$searchpath[,'h'],
			emp_mcr$plus$searchpath[,'h'])
	e0all <- e0all[hall>=0]
	e1all <- e1all[hall>=0]
	plot(y=(e1all/e0all),x=e0all, 
		xlab='original loss',ylab='MR', ...)
	abline(h=1,lty=3)
	if(!is.na(eps)) abline(v=eps)
	legend('topleft',c('MR=1','bounds','eps'),lty=c(3,2,1))

	plot(x=e0all,y=e1all,
		xlab='original loss',ylab='permuted loss', asp=1, ...)
	if(!is.na(eps)) abline(v=eps)

	#Show boundaries from search lemmas
	abline(b=1,a=0,lty=3)
	abline(b=-emp_mcr$minus$gam_opt,a=emp_mcr$minus$full_model_output$h_value,lty=2)
	abline(b=-1/emp_mcr$plus$gam_opt,a=emp_mcr$plus$full_model_output$h_value/emp_mcr$plus$gam_opt,lty=2)

}

# adhoc method to get MCR+ and MCR- from stock optimization methods
get_bound_tightness <- function(s, eps, get_e0, get_e1, get_norm, reg_threshold, start, threshold_tol = 10^-10, K_start = 50, short_nmk_lim = 20,long_nmk_lim = 5000, maxit_SA = 400){

	invalid_constr <- TRUE
	max_iter <- 20
	iter <- 1
	K <- K_start
	while(invalid_constr & iter < max_iter){
		obj_fun <- function(theta){
			e0_t <- get_e0(theta)
			e1_t <- get_e1(theta)
			norm_t <- get_norm(theta)
			penalty <- 0
			if(e0_t > eps) 				penalty <- penalty + K *(e0_t - eps)
			if(norm_t > reg_threshold) 	penalty <- penalty + K *(norm_t - reg_threshold)

			if(s== -1) return( e1_t/e0_t + penalty)
			if(s==  1) return(-e1_t/e0_t + penalty)
		}

		warm_start <- nmk(par=start, fn=obj_fun,control=list(maxfeval=100))$par

		### we always want SANN to give the same answer, regardless of when it's run in a script
			old_seed <- rnorm(1)*1000000
			set.seed(0) 
				soln1_start_point <- optim(par=warm_start,
					fn=function(proposal){
						nmk(par=proposal, fn=obj_fun,control=list(maxfeval=short_nmk_lim))$value
					},
					method='SANN',
					control=list(maxit=maxit_SA, parscale= abs(warm_start))
				)
			set.seed(old_seed)
		### 

		nmk_soln <- nmk(par=soln1_start_point$par, fn=obj_fun,control=list(maxfeval=long_nmk_lim))
		theta_soln <- nmk_soln$par
		# str(nmk_soln)
		
		iter <- iter+1
		invalid_constr <- (get_norm(theta_soln) > reg_threshold + threshold_tol) |
						  (get_e0(theta_soln)   > eps + threshold_tol)
		K <- max(K,2)*10
	}

	e0_sol <- get_e0(theta_soln)
	e1_sol <- get_e1(theta_soln)
	norm_sol <- get_norm(theta_soln)
	return(list(
			par = theta_soln,
			norm = norm_sol,
			e0 = e0_sol,
			e1 = e1_sol,
			MR = e1_sol / e0_sol,
			full_soln =nmk_soln
		))

}
