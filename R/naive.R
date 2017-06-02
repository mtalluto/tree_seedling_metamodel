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

	X <- as.matrix(dat[, c('dd0', 'tdiff', 'gsp', 'winp')])
	Y <- dat[,names(naive_covariates)]

	# organize data and scale
	X <- scale(X)  ## minor change from original; scale to unit variance, instead of by the range

	# add columns for interactions
	Xint <- do.call(cbind, sapply(1:ncol(X), function(i) X[,i] * X[,i:ncol(X), drop=FALSE]))
	colnames(Xint) <- paste(colnames(Xint), rep(colnames(X), ncol(X):1), sep='X')
	X <- cbind(X, Xint)

	# select variables (in internal dataset naive_covariates) and run stan
	lapply(names(naive_covariates), function(modname) {
		X[,naive_covariates[[modname]]]
		standat <- list(
			N = nrow(X),
			K = length(naive_covariates[[modname]]),
			x = X[,naive_covariates[[modname]]],
			y = Y[,modname])
		if(print_stan_messages) {
			mod <- stan(file = stanfile, data=standat, iter=iter, chains=chains)
		} else {
			msg <- capture.output(mod <- stan(file = stanfile, data=standat, iter=iter, chains=chains), type='output')
		}
		mod
	})
}
