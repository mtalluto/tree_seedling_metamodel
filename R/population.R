#' Population model data
#' @description Create data list for laplaces demon for survival model
#' @param dat Data frame for the survival model
#' @export
pop_ld_dat <- function(dat)
{
	Y <- dat[,2]
	Nt <- dat[,3]
	time <- dat[,4]
	X <- as.matrix(dat[,c(5:8)])
	J <- ncol(X)

	list(
		X_p = X,
		Y_p = Y,
		Nt_p = Nt,
		t_p = time,
		J = J,
		N = length(Y),
		PGF = function(Data) rnorm(3 + Data$J), # one for intercept, one for log_sigma, one for b_p
		parm.names = LaplacesDemon::as.parm.names(list(alpha_p=0, beta_p=rep(0,J), sigma_p=0, b_p=0)),
		mon.names = c('LP')
	)
}


#' Population model from seedling metamodel
#' Created by Paige E. Copenhaver-Parry  28 January 2017
#' Modified 6 Dec 2017
#' 
#' @param parm Parameter list
#' @param Data Data list
#' @param nested Boolean, is the model nested within another likelihood?
#' 
#' @details If \code{nested} is \code{TRUE}, only the log likelihood and 
#' log posterior are returned, along with the model predictions. Otherwise
#' a list conforming to LaplacesDemon model functions is returned
#' @export
population_lp <- function(parm, Data, nested=FALSE)
{

	## likelihood
	ll <- 0 
	# unpacking params
	alpha_p <- parm[param_index(Data, 'alpha_p')]
	beta_p <- parm[param_index(Data, 'beta_p')]
	#constrain b-P to be positive
	parm[param_index(Data, 'b_p')] <- b_p <- interval(parm[param_index(Data, 'b_p')], 0, Inf)
	# we track log sigma, to avoid having to truncate sigma (it is strictly positive)
	parm[param_index(Data, 'sigma_p')] <- sigma_p <- interval(parm[param_index(Data, 'sigma_p')], 0, Inf)
	# unpack data
	y_p <- Data$Y_p
	x_p <- Data$X_p
	Nt_p <- Data$Nt_p
	t_p <- Data$t_p

	# model
	a <- alpha_p + x_p %*% beta_p
	rhat <- a - b_p * Nt_p # b_p is the density dependence parameter

	ll <- ll + sum(dnorm(y_p, rhat, sigma_p, log = TRUE))

	### priors on parameters
	LP <- ll + dgamma(sigma_p, 2, 1, log=TRUE) + 
		dgamma(b_p, 2, 1, log=TRUE) + 
	 	sum(dcauchy(beta_p, 0, 2.5, log=TRUE)) + 
	 	dcauchy(alpha_p, 0, 2.5, log=TRUE)

	if(nested)
	{
		return(list(ll=ll, LP = LP, r = rhat))
	} else {
		return(list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rnorm(length(rhat), rhat, sigma_p), parm=parm))
	}

}