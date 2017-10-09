#' Population model data
#' @description Create data list for laplaces demon for survival model
#' @param dat Data frame for the survival model
#' @export
pop_ld_dat <- function(dat)
{
	Y <- dat[,2]
	X <- as.matrix(dat[,-(1:3)])
	J <- ncol(X)

	list(
		X_p = X,
		Y_p = Y,
		J = J,
		N = length(Y),
		PGF = function(Data) rnorm(2 + Data$J), # one for intercept, one for log_sigma
		parm.names = LaplacesDemon::as.parm.names(list(alpha_p=0, beta_p=rep(0,J), log_sigma_p=0)),
		mon.names = c('LP')
	)
}


#' Population model from seedling metamodel
#' Created by Paige E. Copenhaver-Parry  28 January 2017
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

	# unpacking data and params
	alpha_p <- parm[param_index(Data, 'alpha_p')]
	beta_p <- parm[param_index(Data, 'beta_p')]
	# we track log sigma, to avoid having to truncate sigma (it is strictly positive)
	sigma_p <- exp(parm[param_index(Data, 'log_sigma_p')])
	y_p <- Data$Y_p
	x_p <- Data$X_p

	# model
	rhat <- alpha_p + x_p %*% beta_p

	ll <- ll + sum(dnorm(y_p, rhat, sigma_p, log = TRUE))

	### priors on parameters
	LP <- ll + dgamma(sigma_p, 2, 0.1, log=TRUE) + 
		sum(LaplacesDemon::dhalfcauchy(beta_p, 2.5, log=TRUE)) + 
		LaplacesDemon::dhalfcauchy(alpha_p, 2.5, log=TRUE)

	if(nested)
	{
		return(list(ll=ll, LP = LP, r = rhat))
	} else {
		return(list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rnorm(length(rhat), rhat, sigma_p), parm=parm))
	}

}