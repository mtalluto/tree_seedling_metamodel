#' Naive Distribution Model
#'
#' @description Naive distribution model from seedling metamodel, impemented in rstan.
#' Created by Paige E. Copenhaver-Parry  28 January 2017
#'
#' @param iter number of hmc iterations
#' @param chains number of hmc chains
#' @param repeatable boolean, should we use a fixed seed?
#' @param print_stan_messages boolean, show the (copius) stan messages?
#'
#' @return list of stan models, one per species
#' @import treeSeedlingMetamodelData
#' @export
naive_model_stan <- function(iter = 1000, chains = 1, repeatable = FALSE, print_stan_messages = FALSE)
{
	if(repeatable)
		set.seed(20170328)

	stanfile = system.file('stan', 'naiveSDM.stan', package = 'treeSeedlingMetamodel')

	dat <- seedlings$fia_adults
	# Y is a dataframe, one column per species; X is a list, one dataframe per species
	Y <- dat[,names(naive_covariates)] # species names are located in the names of the naive_covariates list
	X <- prep_naive_covars(dat, naive_covariates)

	# run stan for each species
	lapply(names(naive_covariates), function(modname) {
		standat <- list(
			N = nrow(X[[modname]]),
			K = ncol(X[[modname]]),
			x = X[[modname]],
			y = Y[,modname])
		msg <- capture.output(mod <- stan(file = stanfile, data=standat, iter=iter, chains=chains), type='output')
		if(print_stan_messages)
			print(msg)
		mod
	})
}

#' convenience function to automate data prep for the naive model
#' @param dat Data frame; expected to be the seedlings$fia_adults data
#' @param covars Covariate list; expected to be the internal dataset naive_covariates
#' @keywords internal
prep_naive_covars <- function(dat, covars)
{
	X <- as.matrix(dat[, c('dd0', 'tdiff', 'gsp', 'winp')])
	# organize data and scale
	X <- scale(X)  ## minor change from original; scale to unit variance, instead of by the range
	# add columns for interactions
	Xint <- do.call(cbind, sapply(1:ncol(X), function(i) X[,i] * X[,i:ncol(X), drop=FALSE]))
	colnames(Xint) <- paste(colnames(Xint), rep(colnames(X), ncol(X):1), sep='X')
	Xall <- cbind(X, Xint)

	# pull out the needed covars for each species
	xx <- sapply(names(covars), function(modname) {
		ret <- Xall[,covars[[modname]]]
	}, simplify=FALSE, USE.NAMES = TRUE)
	attr(xx, "scaled:center") <- attr(X, "scaled:center")
	attr(xx, "scaled:scale") <- attr(X, "scaled:scale")
	xx
}


#' Naive Distribution Model
#'
#' @description Naive distribution model from seedling metamodel, impemented in LaplacesDemon.
#' Original model by Paige E. Copenhaver-Parry  28 January 2017
#'
#' @param method Estimation method; either \code{'laplace'} for Laplace approximation or \code{'metropolis'}
#'      for Metropolis-Hastings (currently not implemented)
#' @param iter number of mcmc iterations (only for \code{method='metropolis'})
#' @param chains number of mcmc chains (only for \code{method='metropolis'})
#' @param repeatable boolean, should we use a fixed seed?
#'
#' @return
##' @import treeSeedlingMetamodelData
#' @export
naive_model_LD <- function(method=c('laplace', 'metropolis'), iter=1000, chains=1, repeatable = FALSE)
{
	method <- match.arg(method)
	if(repeatable)
		set.seed(20170328)

	dat <- seedlings$fia_adults
	Y <- dat[,names(naive_covariates)] # species names are located in the names of the naive_covariates list
	X = prep_naive_covars(dat, naive_covariates)

	# specify data for LD
	ldData <- lapply(names(naive_covariates), function(modname) {
		J <- ncol(X[[modname]]) + 1 # add one for intercept
		parm.names <- LaplacesDemon::as.parm.names(list(beta=rep(0,J)))
	 	list(
			X = as.matrix(cbind(1, X[[modname]])),
			Y = Y[,modname],
			J = J,
			PGF = function(Data) rnorm(Data$J),
			parm.names = parm.names,
			pos.beta = grep("beta", parm.names),
			mon.names = c('LP')
		)
	})
	names(ldData) <- names(naive_covariates)


	## LD Model
	Model <- function(parm, Data)
	{
		beta <- parm[Data$pos.beta]
		x <- Data$X
		y <- Data$Y
		probs <- plogis(x %*% beta)
		## bernoulli likelihood: p where y == 1; 1-p where y == 0
		ll <- sum(dbinom(y, 1, probs, log=TRUE))

		# prior on betas
		# using a cauchy 0,2.5
		LP <- ll + sum(dcauchy(beta, 0, 2.5, log=TRUE))

		list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rbinom(length(probs), 1, probs), parm=parm)
	}

	lapply(names(naive_covariates), function(modname)
		LaplacesDemon::LaplaceApproximation(Model, LaplacesDemon::GIV(Model, ldData[[modname]], PGF=TRUE),
			ldData[[modname]], Iterations = 1000))

}

