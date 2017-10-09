#' Survival model data
#' @description Create data list for laplaces demon for survival model
#' @param dat Data frame for the survival model
#' @export
survival_ld_dat <- function(dat)
{
	Y <- dat[,2]
	X <- as.matrix(dat[,c('sizet', 'dd0', 'tdiff')])
	J <- ncol(X)

	list(
		X_s = X,
		Y_s = Y,
		J = J,
		N = length(Y),
		PGF = function(Data) rnorm(1 + Data$J),
		parm.names = LaplacesDemon::as.parm.names(list(alpha_s = 0, beta_s=rep(0,J))),
		mon.names = c('LP')
	)
}

#' Survival model
#' @description Survival model from seedling metamodel
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
survival_lp <- function(parm, Data, nested=FALSE)
{
	ll <- 0

	# unpack parameters and data
	alpha_s <- parm[param_index(Data, 'alpha_s')]
	beta_s <- parm[param_index(Data, 'beta_s')]
	x_s <- Data$X_s
	y_s <- Data$Y_s

	probs <- plogis(alpha_s + x_s %*% beta_s)
	## bernoulli likelihood: p where y == 1; 1-p where y == 0
	ll <- ll + sum(dbinom(y_s, 1, probs, log=TRUE))
	LP <- ll + sum(LaplacesDemon::dhalfcauchy(beta_s, 2.5, log=TRUE)) + LaplacesDemon::dhalfcauchy(alpha_s, 5, log=TRUE)

	if(nested)
	{
		return(list(ll=ll, LP = LP, probs = probs))
	} else {
		return(list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rbinom(length(probs), 1, probs), parm=parm))
	}

}



