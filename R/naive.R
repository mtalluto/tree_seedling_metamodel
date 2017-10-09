#' Naive model data
#' @description Create data list for laplaces demon
#' @param dat Data frame with species data for the SDM
#' @export
naive_ld_dat <- function(dat)
{
	Y <- dat[,1]
	X <- as.matrix(dat[,-c(1:3)])
	J <- ncol(X)

	list(
		X_n = X,
		Y_n = Y,
		J_n = J,
		N = length(Y),
		PGF = function(Data) rnorm(1+Data$J),
		parm.names = LaplacesDemon::as.parm.names(list(alpha_n = 0, beta_n=rep(0,J))),
		mon.names = c('LP')
	)
}


#' Naive distribution model
#' Naive distribution model from seedling metamodel
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
naive_lp <- function(parm, Data, nested=FALSE)
{
	ll <- 0

	# unpack parameters and data
	alpha_n <- parm[param_index(Data, 'alpha_n')]
	beta_n <- parm[param_index(Data, 'beta_n')]
	x_n <- Data$X_n
	y_n <- Data$Y_n

	probs <- plogis(alpha_n + x_n %*% beta_n)
	## bernoulli likelihood: p where y == 1; 1-p where y == 0
	ll <- ll + sum(dbinom(y_n, 1, probs, log=TRUE))
	# guard against numerical problems with plogis
	if(is.infinite(ll)) ll <- -1e7
	LP <- ll + sum(LaplacesDemon::dhalfcauchy(beta_n, 2.5, log=TRUE)) + LaplacesDemon::dhalfcauchy(alpha_n, 5, log=TRUE)

	if(nested)
	{
		return(list(ll=ll, LP = LP, probs = probs))
	} else {
		return(list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rbinom(length(probs), 1, probs), parm=parm))
	}
}

