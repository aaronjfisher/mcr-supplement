

library(MASS)
library(Rcgmin)
library(dfoptim)
library(optimr)


###### Weighted loss functions, objective, gradient
hinge_pred <- function(pred, y, m=1){
	(m-pred*y)*(m-pred*y>0)
}
hinge_w <- function(w, X, ...){
	hinge_pred(pred=w[1]+ X %*% w[-1],...)
}
get_hinge_grad <- function(X,y,w, reg_threshold, reg_matrix, K, case.weights){
	pos_ind <- c(hinge_w(w=w, X=X, y=y, m=1) > 0)
	#1/n sum (m-yi xi w), deriv is
	#1/n sum -yi xi where xi is a vector
	penalty_grad <- rep(0,length(w))
	norm_w <- t(w[-1])%*%reg_matrix %*% w[-1]
	if(norm_w > reg_threshold) penalty_grad[-1] <- K * (t(reg_matrix) + reg_matrix) %*% w[-1]

	n_y_pi_cw <- -y * pos_ind * case.weights 
	err_pos <- c(sum(n_y_pi_cw), crossprod(n_y_pi_cw, X)) #first entry is for the intercept

	err_pos + penalty_grad
}
get_hinge_obj <- function(X,y,w, reg_threshold, reg_matrix, K, case.weights){
	un_penalized <- sum(
			case.weights * 
			c(hinge_w(w=w, X=X, y=y, m=1))
		)
	norm_w <- w[-1] %*% reg_matrix %*% w[-1]
	penalty <- 0
	if(norm_w>reg_threshold) penalty <- K* (norm_w-reg_threshold)
	un_penalized + penalty
}


##### Workcheck
	if(FALSE){ 
		library(numDeriv)

		p <- 15
		N <- 50
		X <- matrix(rnorm(N*p),N,p)
		B_gen <- c(1,(1:p)/p)
		logity <- c(cbind(1,X) %*% B_gen)
		y <- rbinom(N,1,prob = 1/(1+exp(-logity)))
		K <- 100
		w_eval <- c(0,(1:p)/p)
		reg_matrix <- matrix(log(1+1:p^2), p, p )
		crossprod(w_eval)

		for(test_cw in 1:2){
		for(test_reg in 1:2){
		for(test_w in 1:2){
			
			if(test_cw==1) cw <- rnorm(N)
			if(test_cw==2) cw <- runif(N)
			
			if(test_w==1) w_eval <- c(0,(1:p)/p)
			if(test_w==2) w_eval <- B_gen
			
			if(test_reg==1) regt <- Inf
			if(test_reg==2) regt <- crossprod(w_eval)*.75

			fn_K <- function(w) get_hinge_obj(X=X,y=y,w=w, reg_threshold=regt, reg_matrix=reg_matrix, K=K, case.weights=cw)
			gr_K <- function(w) get_hinge_grad(X=X,y=y,w=w, reg_threshold=regt, reg_matrix=reg_matrix, K=K, case.weights=cw)
			fn_K(w_eval)
			# gr_K(w_eval)

			diffs <- abs(grad(fn_K,w_eval) - gr_K(w_eval))
			str(max(abs(diffs)) / max(abs(gr_K(w_eval))))
		}}}
	}
#####




# Optimize a linear classifier with possibly negative observation weights
optimize_hinge_general <- function(y, X, reg_matrix=diag(p), n=length(y), case.weights = rep(1/n, n), reg_threshold, K=1000, p = dim(X)[2], start = rep(0,p+1), constr_tol = 10^-4, method='grad-based', sparse_warning=TRUE, maxit_SA = 1000, extra_step_NM=TRUE,
	short_cg_lim =15, long_cg_lim =500){


	if(reg_threshold<=0) stop('invalid threshold')

	##### Set contant variables to have coefficient = 0, only analyze remaining variables in body of this function
		is_intercept <- apply(X,2,var)==0
		if(any(is_intercept) & sparse_warning) warning(sum(is_intercept), ' constant variables detected. Assigning these to have coefficient 0.')
		X_analyze <- X[,!is_intercept]
		reg_matrix_analyze <- reg_matrix[!is_intercept,!is_intercept]
		start_analyze <- start[c(TRUE,!is_intercept)]
		X_names <- colnames(X)

		rm(list=c('X','p','start'))
	#####

	invalid_constr <- TRUE 
	max_iter <- 20
	iter <- 1
	while(invalid_constr & iter < max_iter){
	
		fn_K <- function(w){
			get_hinge_obj(X=X_analyze,y=y,w=w, reg_threshold=reg_threshold, reg_matrix=reg_matrix_analyze, K=K, case.weights=case.weights)
		}
		gr_K <- function(w){
			get_hinge_grad(X=X_analyze,y=y,w=w, reg_threshold=reg_threshold, reg_matrix=reg_matrix_analyze, K=K, case.weights=case.weights)
		}
		
		# grad_based_step <- function(proposal, maxit_cg){Rcgmin(par=proposal, fn=fn_K, gr=gr_K,control=list(maxit=maxit_cg))}
		# grad_based_step <- function(proposal, maxit_cg){optim(par=proposal, fn=fn_K, gr=gr_K, method='BFGS',control=list(maxit=maxit_cg))}
		# grad_based_step <- function(proposal, maxit_cg){optim(par=proposal, fn=fn_K, gr=gr_K, method='CG',control=list(maxit=maxit_cg))}
		# grad_based_step <- function(proposal, maxit_cg){optimr(par=proposal, fn=fn_K, gr=gr_K, method='CG',control=list(maxit=maxit_cg))}
		grad_based_step <- function(proposal, maxit_cg){optimr(par=proposal, fn=fn_K, gr=gr_K, method='BFGS',control=list(maxit=maxit_cg))}
			#in ad hoc tests, this solver appeared to work best.
		

		soln <-
		warm_start <- grad_based_step(start_analyze, short_cg_lim)

		if(!method %in% c('grad-based','grad-based-plus-SA')){
			warning('non-standard solving method')
			soln <-
			soln1 <- optimr(par=warm_start$par, fn=fn_K, method=method, ...)
		}

		if(method=='grad-based-plus-SA'){			

			parscale_SA <- pmax(
				max(abs(warm_start$par))*10^-8,
				abs(warm_start$par)
			)
			if(all(parscale_SA==parscale_SA[1])) parscale_SA <- rep(1,length(warm_start$par))
			### we would like SANN to give the same answer regardless of where it is run during a script -- change seed to ensure this.
			old_seed <- rnorm(1)*1000000
			set.seed(0) 
				soln1_start_point <- optim(par=warm_start$par,
					fn=function(proposal){
						fn_K(grad_based_step(proposal, short_cg_lim)$par)
					}, method='SANN',
					control=list(
						maxit=maxit_SA,
						parscale= parscale_SA))
			set.seed(old_seed)
			### 

			soln <-
			soln1 <- grad_based_step(proposal=soln1_start_point$par, maxit_cg=long_cg_lim)
		}
		if(method=='grad-based'){
			soln <-
			soln1 <- grad_based_step(proposal=start_analyze, maxit_cg=long_cg_lim)
		}

		if(extra_step_NM) soln <- nmk(par=soln1$par, fn=fn_K)

		out <- soln$par
		
		invalid_constr <- t(out[-1])%*%reg_matrix_analyze %*% out[-1] > reg_threshold * (1+constr_tol)
		K <- max(c(K,1000))*1000 #increase penalty if constraint is not met, and optimization must be re-run
		iter <- iter+1
	}
	if(invalid_constr) stop('Invalid constraint value')

	### Assign sparse categories to have 0 coefficient
	coef_vec <- c(NA,rep(0,length(is_intercept)))
	names(coef_vec) <- c('intercept', X_names)
	coef_vec[c(TRUE,!is_intercept)] <- out
	### 

	coef_vec
}






get_suff_stats_linear_hinge <- function(reg_matrix, reg_threshold, start, ...){
	# compute loop invariant quantities for MCR binary search.
	fs <- get_full_sample(...)
	c(fs, list(reg_matrix=reg_matrix, start = start, reg_threshold=reg_threshold))
}
get_e0_hinge <- function(w,sst){
	mean(hinge_w(w,
		X=sst$full_X[!sst$is_perm,],
		y=sst$full_y[!sst$is_perm]))
}
get_e1_hinge <- function(w,sst){
	mean(hinge_w(w,
		X=sst$full_X[sst$is_perm,],
		y=sst$full_y[sst$is_perm]))
}

get_h_min_linear_hinge <- function(s=-1,gam,suff_stats,upper_method, mod1=NA,mod2=NA){
	# minimize linear combinations of e_orig & e_switch

	sst <- suff_stats
	if(!s %in% c(-1, 1)) stop('s input error')
	
	weights_gam <- rep(NA,length(sst$is_perm))
	if(s==-1){ #MCR-
		weights_gam[ sst$is_perm] <- 1/sum(sst$is_perm) #avoids using n here, since sample size may != n if e_divide approach is used.
		weights_gam[!sst$is_perm] <- gam/sum(!sst$is_perm)
	}else{ #MCR+
		weights_gam[ sst$is_perm] <- gam/sum(sst$is_perm)
		weights_gam[!sst$is_perm] <- 1/sum(!sst$is_perm)
		if(gam!=0 & upper_method != 'grad-based-plus-SA') warning('Optimization is approximate, multiple starts are recommended')
	}
	# sum(weights_gam[ sst$is_perm])
	# sum(weights_gam[!sst$is_perm])
	
	if(s== 1) opt_method <- upper_method #for MCR+
	if(s==-1) opt_method <- 'grad-based'
	
	start_w <- sst$start
	if(!all(is.na(mod1))){
		start_w <- (mod1$w_gam + mod2$w_gam)/2
		#mod1 & mod2 are outputs from two previous optimizations, used as a starting point.
	}

	w_gam <- optimize_hinge_general(y=sst$full_y, X=sst$full_X, case.weights=weights_gam, start = start_w, reg_threshold= sst$reg_threshold, method=opt_method, extra_step_NM = TRUE)

	pred_gam <- cbind('intercept'=1,as.matrix(sst$full_X))%*%w_gam
	e_orig_gam <- mean(hinge_pred(
		pred_gam[!sst$is_perm], sst$full_y[!sst$is_perm]))
	e_switch_gam <-  mean(hinge_pred(
		pred_gam[ sst$is_perm], sst$full_y[ sst$is_perm]))
	h_value <- crossprod(hinge_pred(pred_gam, sst$full_y), weights_gam)

	###### redundant workchecks	
		if(s==-1 & h_value - (e_switch_gam+e_orig_gam*gam) > 10^-10) stop('h_value calculation error')
		if(s== 1 & h_value - (e_switch_gam*gam+e_orig_gam) > 10^-10) stop('h_value calculation error')
	######  

	return(list(
		h_value = h_value,
		e0_value  =e_orig_gam,
		e1_value  =e_switch_gam,
		Reg = t(w_gam[-1])%*% sst$reg_matrix %*%w_gam[-1],
		w_gam = w_gam
	))

}


