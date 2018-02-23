rm(list = ls())

library('RColorBrewer')
# library('gmum.r')
library('ggplot2')
library('SparseM')
library('Matrix')
library('pbapply')


script_folder <- '../general-scripts/'

source(paste0(script_folder,'binary_search_mcr.r'))
source(paste0(script_folder,'hinge_mcr.r'))

if(!dir.exists('plots')) dir.create('plots')

set.seed(0)






dat_type <- 'circles'


######### polynomial features of different degrees
to_poly1 <- function(x){
	c(
		'x1'=x[1],
		'x2'=x[2]
	)
}
to_poly2 <- function(x){
	c(
		to_poly1(x),
		'x1sq'=x[1]^2,
		'x2sq'=x[2]^2,
		'x1x2'=x[1]*x[2]
	)
}
to_poly3 <- function(x){
	c(
		to_poly2(x),
		x1cu=x[1]^3,
		x2cu=x[2]^3,
		x1sqx2=(x[1]^2)*x[2],
		x1x2sq=x[1]*(x[2]^2)
		)
}




########## generate raw coordinate values
gen_Xy <- function(n, dat_type){ 
	p_base <- 2

	y <- rep(c(-1,1),times = n/2)[sample(1:n)]
	if(any(y==0)) stop('not correct y specification') #redundant note on how y is coded

	X <- matrix(NA,n,p_base)
	X_base <-  matrix(rnorm(n*p_base),n,p_base)/3

	if(dat_type=='linear'){
		separation <- matrix(NA,n,p_base)
		separation[,1] <- 0.4
		separation[,2] <- 0.2
		X <- X_base + separation * y
	}

	if(dat_type=='circles'){
		circ1 <- array(NA,dim=c(n,2))
		dimnames(circ1) <- list('obs'=1:n,'position'=c('x1','x2'))
		circ_position <- runif(n,-pi,pi)
		circ1[,'x1'] <- cos(circ_position)
		circ1[,'x2'] <- sin(circ_position)
		circ2 <- circ1 * 0
		# plot(rbind(circ1,circ2))

		X[y== 1,] <- (X_base + circ1)[y== 1,] #1 goes on the outside
		X[y==-1,] <- (X_base + circ2)[y==-1,] #-1 goes on the inside
	}

	if(dat_type=='halfmoons'){
		halfmoon1 <- array(NA,dim=c(n,2))
		dimnames(halfmoon1) <- list('obs'=1:n,'position'=c('x1','x2'))
		halfmoon2 <- halfmoon1
		start_point <- -pi/4
		moon_position <- runif(n,start_point,start_point+pi)
		halfmoon1[,'x1'] <- cos(moon_position)
		halfmoon1[,'x2'] <- sin(moon_position)
		halfmoon2[,'x1'] <- cos(moon_position+pi) + cos(start_point) + cos(start_point+pi/2)/2
		halfmoon2[,'x2'] <- sin(moon_position+pi) + sin(start_point) + sin(start_point+pi/2)/2
		# plot(rbind(halfmoon1,halfmoon2))

		X[y== 1,] <- (X_base + halfmoon1)[y== 1,]
		X[y==-1,] <- (X_base + halfmoon2)[y==-1,]
	}
	return(list(X=X, y=y))
}




######### Simulate data
n_train <- 300
n_test <- 300
p1_base <- 1 #variable of interest, to switch
to_poly_deg <- to_poly3

Xy_base_test <- gen_Xy(n_test,dat_type)
Xy_base_train <- gen_Xy(n_train,dat_type)
Xy_base_train$X[1,]
round(to_poly_deg(Xy_base_train$X[1,]),4) #example of covariates in full feature space

X_train <- t(apply(Xy_base_train$X, 1, to_poly_deg))
X_test  <- t(apply(Xy_base_test$X, 1, to_poly_deg))
y_train <- Xy_base_train$y
y_test  <- Xy_base_test$y
(p_poly <- dim(X_train)[2])


plot_train_dat <- function(col=c('blue','black')[(y_train==1)+1],...){
	plot(X_train[,c('x1','x2')],col=col,pch=c(1,4)[(y_train==1)+1],...)
}
# plot_train_dat()
plot_test_dat <- function(col=c('blue','black')[(y_test==1)+1],...){
	plot(X_test[,c('x1','x2')],col=col,pch=c(1,4)[(y_test==1)+1],...)
}
# plot_test_dat()




#  _                    _
# | |                  (_)
# | |_   _ __    __ _   _   _ __
# | __| | '__|  / _` | | | | '_ \
# | |_  | |    | (_| | | | | | | |
#  \__| |_|     \__,_| |_| |_| |_|



####################   Cross-validate

# Start by fitting a basic model
system.time({w_unconstrained <- optimize_hinge_general(y=y_train, X=X_train, reg_threshold=Inf, extra_step_NM=TRUE ,method='grad-based', short_cg_lim =15, long_cg_lim =500)})
mean(hinge_w(w=w_unconstrained,X=X_train,y=y_train))

J <- 60
w_start_cv <- rep(0,p_poly+1) #Starting point for optimizations. Don't use unrestricted solution (w_unconstrained) because this uses out-of-fold data.
(reg_unlim <- c(crossprod(w_unconstrained[-1])))
reg_seq <- reg_unlim * (2^c(0,seq(-3,0.3,length=J-1)))
reg_seq <- reg_seq[order(reg_seq)]
cv_report <- t(pbsapply(reg_seq,function(rr){
	CV_mod(y=y_train,X=X_train,n_folds=10,
	fit_fun = function(y,X){optimize_hinge_general(
		y=y, X=X, reg_threshold=rr, start=w_start_cv,
		method='grad-based',
		sparse_warning=FALSE, 
		extra_step_NM = TRUE,
		long_cg_lim =500,
		maxit_SA=100
		)},
	report_fun = function(y,X,mod){
		c('err'=mean(hinge_w(y=y,X=X,w=mod)),
		  'reg'=crossprod(mod[-1]))
	}, fit_in_sample=TRUE)
}))

(reg_threshold <- reg_seq[which(cv_report[,'cv.err'] == min(cv_report[,'cv.err']))][1])

pdf(paste0('plots/',Sys.Date(),'_10fold_CV_plot.pdf'))
	plot(reg_seq,cv_report[,'cv.err'],type='o',col='orange',lwd=2, ylab='Loss', xlab='regularization threshold', main='Cross-validation on training data', ylim=range(cv_report[,c('cv.err','in_sample.err')]))
	lines(reg_seq,cv_report[,'in_sample.err'],type='o',col='darkgreen',lwd=2)
	abline(h=min(cv_report[,'cv.err']), lty=3)
	abline(v=reg_threshold,lty=2)
	legend('topright', pch=1, lwd=2, 
		   c('training loss','CV loss'),
		col=c('darkgreen','orange'))
dev.off()

# Select a reference model (w_tr)
w_tr <- optimize_hinge_general(y=y_train, X=X_train, reg_threshold=reg_threshold, extra_step_NM=TRUE)
crossprod(w_tr[-1])/reg_threshold # workcheck
crossprod(w_unconstrained[-1])/reg_threshold # workcheck

err_est_w_tr <- min(cv_report[,'cv.err'])
####################





#################### AR - Train models that only use x1 or x2
#tag-cv-AR-sim
noX1 <- which(!grepl('x1',names(to_poly_deg(Xy_base_train$X[1,]))))
noX2 <- which(!grepl('x2',names(to_poly_deg(Xy_base_train$X[1,]))))
w_poly_drop2 <- 
w_poly_drop1 <- rep(0, 1+p_poly)
names(w_poly_drop1) <-
names(w_poly_drop2) <- c('intercept',names(to_poly_deg(Xy_base_train$X[1,])))

w_poly_drop1[c(1,1+noX1)] <- optimize_hinge_general(y=y_train, X=as.matrix(X_train[,noX1]), reg_threshold=reg_threshold)
w_poly_drop2[c(1,1+noX2)] <- optimize_hinge_general(y=y_train, X=as.matrix(X_train[,noX2]), reg_threshold=reg_threshold)
####################





# Remove information from training data, with the exception of:
	# the regularization threshold,
	# the reference model w_tr,
	# AR models
	# the estimated cv error of w_tr.
rm(w_unconstrained, Xy_base_train, X_train, y_train, cv_report)


#  _                  _
# | |                | |
# | |_    ___   ___  | |_
# | __|  / _ \ / __| | __|
# | |_  |  __/ \__ \ | |_
#  \__|  \___| |___/  \__|


#### Calculate loop-invariant statistics for MCR binary search.
# here, `sst` stands for "sufficient statistics"
sstpoly <- 
sstpoly_pre <- get_suff_stats_linear_hinge(y=y_test,X=Xy_base_test$X,p1=p1_base,times=2, reg_matrix=diag(p_poly), start = w_start_cv, reg_threshold=reg_threshold)
sstpoly$full_X <- t(apply(sstpoly_pre$full_X, 1, to_poly_deg))
str(sstpoly)



#### Set the value of eps_abs -- epsilon on the absolute scale rather than the relative scale
(eps_poly <- 
	mean(hinge_w(y=y_test,X=X_test,w=w_tr)) + # E(loss(Ref-model))
	err_est_w_tr * 0.10) #approx 10% of the expected loss of training model (Ref-model), estimated from training data.
		#tag-cv-est-eps-10-percent
		# tag-choice-of-eps

#### Compute MR & MCR
getMRpolyw <- function(w){
	get_e1_hinge(w,sstpoly)/get_e0_hinge(w,sstpoly)
}
round(getMRpolyw(w_tr),2) #MR of reference model; #tag f-ref-mr

system.time({
	hinge_poly_mcr <- get_empirical_MCR(
			eps=eps_poly,
			get_h_min=function(...){get_h_min_linear_hinge(upper_method='grad-based-plus-SA',...)},
			suff_stats=sstpoly)
})
hinge_poly_mcr$range
str(lapply(hinge_poly_mcr[-1], function(z) {z$approximation_error_linear})) # tag-approx-sim
str(lapply(hinge_poly_mcr[-1], function(z) {z$approximation_error_search}))
# str(hinge_poly_mcr) 





#        _       _
#       | |     | |
#  _ __ | | ___ | |_
# | '_ \| |/ _ \| __|
# | |_) | | (_) | |_
# | .__/|_|\___/ \__|
# | |
# |_|


# Plot the classification boundary associated with a model
plot.poly.w <- function(w, hull.omit=FALSE, ...){
	nn <- 100
	sX <- seq(min(Xy_base_test$X),max(Xy_base_test$X),length=nn)
	z_contour <- matrix(NA, nn,nn)
	for(i in 1:nn){
	for(j in 1:nn){
		z_contour[i,j] <- w %*% c(1,to_poly_deg(sX[c(i,j)]))
	}}
	contour(x=sX,y=sX,z_contour,levels=0, drawlabels=FALSE, add=TRUE,...)
}


w_best <- optimize_hinge_general(y=y_test, X=X_test, reg_threshold=reg_threshold, extra_step_NM=TRUE)

# Function to generate a set of models, among which a certain percentage must be well-performing. 
#based on a matrix of reference points.
#' @param eps_R loss threshold on absolute scale.
generate_R_set <- function(sst, eps_R, w_best, w_ref_mat=NULL, nreps=1000, prop_good=0.8, shrink_factor = 0.2){

	p <- length(w_best)-1
	is_good <- rbinom(n=nreps,1,p=prop_good) #which models should be well performing (in the Rashomon set)
	e0_sim <- e1_sim <- rep(NA,nreps)
	w_mat <- matrix(NA,nreps, p+1)

	if( mean(hinge_w(w_best, X=sst$full_X[!sst$is_perm,], y=sst$full_y[!sst$is_perm])) > eps_R) stop('invalid R set')

	all_references_points <- rbind(w_best,w_ref_mat)
	picks <- rep(1:dim(all_references_points)[1], length=nreps)
	var_ref <- apply(all_references_points,2,var)
	if(any(is.na(var_ref))) var_ref <- rep(1, p+1)

	for(j in 1:nreps){
		w_ref_j <- all_references_points[picks[j],]
		
		w_sim_j <- w_ref_j + mvrnorm(1,
			mu=rep(0,p+1),
			Sigma=diag(var_ref))
		e0j <-   mean(hinge_w(w_sim_j, X=sst$full_X[!sst$is_perm,], y=sst$full_y[!sst$is_perm]))
		iter <- 1
		check_if_too_big <- function(w){
			c((w[-1] %*% sst$reg_matrix %*% w[-1]) > sst$reg_threshold )
		}
		while(( e0j > eps_R | check_if_too_big(w_sim_j)) & iter<1000 & is_good[j]==1){
			if(e0j > eps_R) w_sim_j <- shrink_factor * w_best + (1-shrink_factor) * w_sim_j
			if(check_if_too_big(w_sim_j)) w_sim_j <- 0 + (1-shrink_factor) * w_sim_j
			iter <- iter+1
			e0j <-   mean(hinge_w(w_sim_j, X=sst$full_X[!sst$is_perm,], y=sst$full_y[!sst$is_perm]))
		}

		e0_sim[j] <-e0j
		e1_sim[j] <- mean(hinge_w(w_sim_j, X=sst$full_X[sst$is_perm,], y=sst$full_y[sst$is_perm]))
		w_mat[j,] <- w_sim_j
	}
	if(any(e0_sim[is_good==1] > eps_R)) warning('Less than the requested proportion of models may be well performing.')
	return(list(e0=e0_sim, e1=e1_sim, w_mat = w_mat))
}


mypal <- brewer.pal(9,'Paired')

pdf(paste0('plots/',Sys.Date(),'_MCR+AR_1panel_example_poly_hinge_',dat_type,'.pdf'),6,6)
{
	plot_test_dat(col=mypal[c(7,9)][(y_test==1)+1], main=paste0('AR + MCR Example: Polynomial Classifier with\nHinge Loss (',dat_type,')'), xlab='X1', ylab='X2')
	plot.poly.w(w_tr,lwd=3,lty=3)
	plot.poly.w(hinge_poly_mcr$minus$full_model$w,lty=2,lwd=2,col=mypal[4])
	plot.poly.w(hinge_poly_mcr$plus$full_model$w,lty=1,lwd=2,col=mypal[4])
	plot.poly.w(w_poly_drop1,lty=2,lwd=2,col=mypal[2])
	plot.poly.w(w_poly_drop2,lty=1,lwd=2,col=mypal[2])
	legend('bottomright',c('AR-','AR+','MCR-','MCR+','ref'), col=c(mypal[c(2,2,4,4)],'black'),lwd=c(2,2,2,2,3), lty=c(2,1,2,1,3),bg='white')

}
dev.off()

pdf(paste0('plots/',Sys.Date(),'_3panel_example_poly_hinge_',dat_type,'.pdf'),height=5,width=10.5,pointsize=14)
{

	fref <- expression(paste(italic('f')['ref']))
	fpe <-   expression(paste(hat(italic('f'))['+,'][epsilon],))
	fme <-   expression(paste(hat(italic('f'))['-,'][epsilon],))
	f1 <-    expression(paste(hat(italic('f'))[1]))
	f2 <-    expression(paste(hat(italic('f'))[2]))
	x1exp <-    expression(paste(italic('X'))[1])
	x2exp <-    expression(paste(italic('X'))[2])

	layout(matrix(c(1,2,3),nrow=1), widths=c(2,2,1.25))
	par(oma=c(0,0,2,0), mar = c( 5.1,4.1,4.1,2.1))

	plot_test_dat(col=mypal[c(7,9)][(y_test==1)+1], main='', xlab='', ylab='', cex.lab=1.1)
	mtext(x1exp,1,2.8);mtext(x2exp,2,2)
	mtext(text='1) Algorithm reliance (AR)',cex=1.1, line=1.62,adj=0)
	plot.poly.w(w_tr,lwd=3,lty=3)
	plot.poly.w(w_poly_drop1,lty=2,lwd=2,col=mypal[2])
	plot.poly.w(w_poly_drop2,lty=1,lwd=2,col=mypal[2])
	legend('bottomright',c(fref,f1,f2), col=c(1,mypal[c(2,2)]),lwd=c(3,2,2), lty=c(3,1,2),bg='white', cex=1.24, pt.cex=1)




	mtext(text='Example: AR, MCR & MR for polynomial classifiers',outer=TRUE,line=0, font=2, cex=1.15)




	#### R set
	set_size <- 15
	w_ref_mat_poly <- rbind(
		hinge_poly_mcr$minus$full_model$w,
		hinge_poly_mcr$plus$full_model$w)
	set.seed(0) # fix seed now that analysis is done, to keep the same plot contents
	poly_R_set <- generate_R_set(sstpoly, eps_R=eps_poly, w_ref_mat=w_ref_mat_poly, w_best=w_best, nreps=set_size, prop_good=1, shrink_factor=0.2)

	plot_test_dat(col=mypal[c(7,9)][(y_test==1)+1], main='', xlab='', ylab='', cex.lab=1.1)
	mtext(x1exp,1,2.8);mtext(x2exp,2,2)
	mtext(text='2) Model class reliance (MCR)',cex=1.1, line=1.62,adj=0)
	# plot.poly.w(w_tr,lty=1,lwd=2,col='black')
	for(j in 1:set_size){
		plot.poly.w(poly_R_set$w_mat[j,],lty=1,lwd=1,col='lightgray')	
	}
	plot.poly.w(hinge_poly_mcr$minus$full_model$w,lty=2,lwd=2,col=mypal[4])
	plot.poly.w(hinge_poly_mcr$plus$full_model$w,lty=1,lwd=2,col=mypal[4])
	legend('bottomright',c('good',fpe,fme), col=c('darkgray',mypal[c(4,4)]),lwd=c(1,2,2), lty=c(1,1,2),bg='white', cex=1.24, pt.cex=1)
	####



	goodMR <- apply(poly_R_set$w_mat,1,getMRpolyw)
	drop2MR <- getMRpolyw(w_poly_drop2)
	drop1MR <- getMRpolyw(w_poly_drop1)
	highMR <- getMRpolyw(hinge_poly_mcr$plus$full_model$w)
	lowMR <- getMRpolyw(hinge_poly_mcr$minus$full_model$w)
	
	plot(c(), ylim=c(1,0.5*ceiling(max(highMR,drop2MR)*2)),xlim=c(-.07,.9),axes=F,ylab='',xlab='', cex.lab=1.1)
	mtext('3) Model\nreliance (MR)', cex=1.1, line=0.0,adj=0)
	mtext('MR',2,2)
	axis(2,lwd=2)
	segments(x0=rep(0,set_size),x1=rep(0.2,set_size),y0=c(goodMR),y1=c(goodMR),col='darkgray')
	segments(x0=0,x1=0.2,y0=highMR,y1=highMR,lty=1,col=mypal[4],lwd=2)
	segments(x0=0,x1=0.2,y0=lowMR,y1=lowMR,lty=2,col=mypal[4],lwd=2)
	segments(x0=0,x1=0.2,y0=getMRpolyw(w_poly_drop2),y1=getMRpolyw(w_poly_drop2),lty=1,col=mypal[2],lwd=2)
	segments(x0=0,x1=0.2,y0=getMRpolyw(w_poly_drop1),y1=getMRpolyw(w_poly_drop1),lty=2,col=mypal[2],lwd=2)
	segments(x0=0,x1=0.2,y0=getMRpolyw(w_tr),y1=getMRpolyw(w_tr),lty=3,col='black',lwd=3)

	points(x=rep(0.1,set_size),y=c(goodMR),col='darkgray')
	points(x=0.1,y=highMR,pch=25,col=mypal[4],lwd=2)
	points(x=0.1,y=lowMR,pch=24,col=mypal[4],lwd=2)

	points(x=0.1,y=getMRpolyw(w_poly_drop2),col=mypal[2],pch=0,lwd=2)
	points(x=0.1,y=getMRpolyw(w_poly_drop1),col=mypal[2],pch=5,lwd=2)
	points(x=0.1,y=getMRpolyw(w_tr),pch=19,col='black',lwd=2)

	legend('bottomright',c(fref,'good',fpe,fme,f1,f2),
		lwd=c(2,1,2,2,2,2),
		pch=c(19,1,25,24,0,5),bg='white'
		,col=c('black','darkgray',mypal[c(4,4,2,2)]),
		lty=c(3,1,1,2,1,2),
		cex=1.24, pt.cex=1)

}
dev.off()





