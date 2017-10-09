#' @description Create data list for laplaces demon for survival model
#' @param dat Data frame for the survival model
survival_ld_dat <- function(dat)
{
	Y <- dat[,1]
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






#' Seedling survival model
#'
# #' @description Seedling-only survival model, impemented in LaplacesDemon.
# #'
# #' @param method Estimation method; either \code{'laplace'} for Laplace approximation or \code{'metropolis'}
# #'      for Metropolis-Hastings (currently not implemented)
# #' @param iter number of mcmc iterations (only for \code{method='metropolis'})
# #' @param chains number of mcmc chains (only for \code{method='metropolis'})
# #' @param repeatable boolean, should we use a fixed seed?
# #'
# #' @return
# ##' @import treeSeedlingMetamodelData
# #' @export
# survival_model_LD <- function(method=c('laplace', 'metropolis'), iter=1000, chains=1, repeatable = FALSE)
# {
# 	# method <- match.arg(method)
# 	# if(repeatable)
# 	# 	set.seed(20170328)

# 	# dat <- seedlings$survival
# 	# spNames <- unique(dat$species)

# 	# ldData <- sapply(spNames, function(sp)
# 	# {
# 	# 	Data <- dat[dat$species == sp,]
# 	# 	X <- data.frame(sizet = Data$x, dd0 = Data$dd0, tdiff = Data$tdiff)
# 	# 	X <- scale(X)
# 	# 	X <- as.matrix(cbind(1, X))
# 	# 	J <- ncol(X)
# 	# 	parm.names <- LaplacesDemon::as.parm.names(list(beta_surv=rep(0,J)))
# 	#  	list(
# 	# 		X = X,
# 	# 		Y = Data$surv,
# 	# 		J = J,
# 	# 		PGF = function(Data) rnorm(Data$J),
# 	# 		parm.names = parm.names,
# 	# 		pos.beta_surv = grep("beta_surv", parm.names),
# 	# 		mon.names = c('LP')
# 	# 	)
# 	# }, simplify = FALSE, USE.NAMES = TRUE)


# 	# ## LD Model
# 	# Model <- function(parm, Data)
# 	# {
# 	# 	beta_surv <- parm[Data$pos.beta_surv]
# 	# 	x <- Data$X
# 	# 	y <- Data$Y
# 	# 	probs <- plogis(x %*% beta_surv)
# 	# 	## bernoulli likelihood: p where y == 1; 1-p where y == 0
# 	# 	ll <- sum(dbinom(y, 1, probs, log=TRUE))

# 	# 	# prior on beta_surv
# 	# 	LP <- ll + sum(dcauchy(beta_surv, 0, 2.5, log=TRUE))

# 	# 	list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rbinom(length(probs), 1, probs), parm=parm)
# 	# }

# 	# sapply(ldData, function(D)
# 	# 	LaplacesDemon::LaplaceApproximation(Model, LaplacesDemon::GIV(Model, D, PGF=TRUE),
# 	# 		D, Iterations = 1000), USE.NAMES = TRUE, simplify=FALSE)

# }



#' Seedling survival integrated model
#'
#' @description Integrated model combining naive and seedling survival, impemented in LaplacesDemon.
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
survival_model_integrated_LD <- function(method=c('laplace', 'metropolis'), iter=1000, chains=1, repeatable = FALSE)
{
	method <- match.arg(method)
	if(repeatable)
		set.seed(20170328)

	naive_dat <- seedlings$fia_adults
	spNames <- names(naive_covariates)
	naiveX <- prep_naive_covars(naive_dat, naive_covariates)
	naiveY <- naive_dat[,names(naive_covariates)] # species names are located in the names of the naive_covariates list

	survival_dat <- seedlings$survival
	ldData <- sapply(spNames, function(sp)
	{
		X_naive <- as.matrix(cbind(1, naiveX[[sp]]))
		J_naive <- ncol(X_naive)

		D_surv <- survival_dat[survival_dat$species == sp,]
		X_surv <- data.frame(sizet = D_surv$x, dd0 = D_surv$dd0, tdiff = D_surv$tdiff)

		# scale the two matching variables to the same scale as the naive dataset
		X_surv$sizet <- scale(X_surv$sizet)
		X_surv$dd0 <- scale(X_surv$dd0, center = attr(naiveX, "scaled:center")['dd0'], scale = attr(naiveX, "scaled:scale")['dd0'])
		X_surv$tdiff <- scale(X_surv$tdiff, center = attr(naiveX, "scaled:center")['tdiff'], scale = attr(naiveX, "scaled:scale")['tdiff'])
		X_surv <- as.matrix(cbind(1, X_surv))
		J_surv <- ncol(X_surv)

		# create an integration x dataset with the same cols as survival, but with interactions we know about to be used if needed
		X_int <- as.data.frame(X_surv[,which(! colnames(X_surv) %in% 'sizet')])
		if('dd0Xdd0' %in% colnames(X_naive))
			X_int$dd0Xdd0 <- X_int$dd0 * X_int$dd0
		if('tdiffXtdiff' %in% colnames(X_naive))
			X_int$tdiffXtdiff <- X_int$tdiff * X_int$tdiff
		if('tdiffXdd0' %in% colnames(X_naive))
			X_int$tdiffXdd0 <- X_int$tdiff * X_int$dd0
		X_int <- as.matrix(X_int)
		int.indices <- which(colnames(X_naive) %in% colnames(X_int))
		int.names <- paste0("beta_naive[", int.indices, "]")

		parm.names <- LaplacesDemon::as.parm.names(list(beta_surv=rep(0,J_surv), beta_naive=rep(0,J_naive), threshold=0))
	 	list(
	 		X_naive = X_naive,
	 		Y_naive = naiveY[,sp],
	 		J_naive = J_naive,
	 		X_surv = X_surv,
			Y_surv = D_surv$surv,
			J_surv = J_surv,
			X_int = X_int,
			PGF = function(Data) rnorm(Data$J_naive + Data$J_surv + 1),
			parm.names = parm.names,
			pos.beta_surv = grep("beta_surv", parm.names),
			pos.threshold = grep("threshold", parm.names),
			integrated_betas =which(parm.names %in% int.names),
			pos.beta_naive = grep("beta_naive", parm.names),
			mon.names = c('LP')
		)
	}, simplify = FALSE, USE.NAMES = TRUE)

	test <- sapply(ldData, function(D)
	 	LaplacesDemon::LaplaceApproximation(integrated_survival, LaplacesDemon::GIV(integrated_survival, D, PGF=TRUE),
	 										D, Iterations = 1000), USE.NAMES = TRUE, simplify=FALSE)

}


integrated_survival <- function(pars, dat)
{
	ll <- 0
	## SUBMODELS
	# first submodel - here we set up the log posterior
	# survival model
	beta_surv <- pars[dat$pos.beta_surv]
	x_surv <- dat$X_surv
	y_surv <- dat$Y_surv
	probs_surv <- plogis(x_surv %*% beta_surv)
	# bernoulli likelihood: p where y == 1; 1-p where y == 0
	ll <- sum(dbinom(y_surv, 1, probs_surv, log=TRUE))


	# naive model - compute the contribute to the LP of the naive SDM
	beta_naive <- pars[dat$pos.beta_naive]
	x_naive <- dat$X_naive
	y_naive <- dat$Y_naive
	probs_naive <- plogis(x_naive %*% beta_naive)
	## bernoulli likelihood: p where y == 1; 1-p where y == 0
	ll <- ll + sum(dbinom(y_naive, 1, probs_naive, log=TRUE))


	# integrated part of model
	# compute a pseudo y dataset from the survival data using a threshold model
	# the threshold is a parameter that will be estimated from data
	threshold <- plogis(pars[dat$pos.threshold])
	y_int_surv <- as.integer(probs_surv >= threshold)
	x_int_surv <- dat$X_surv[,dat$integrated_vars]
	beta_integrated <- pars[dat$integrated_betas] ## these should overlap with some of the params in beta_naive
	probs_integrated <- plogis(x_int_surv %*% beta_integrated) # note we only use the params we actually care about for integration; other params are set to 0
	ll <- ll + sum(dbinom(y_int_surv, 1, probs_integrated, log = TRUE))
	## NO PRIOR on beta_integrated, because they are the same parameters as in beta_naive

	## END OF LIKELIHOOD SECTION
	LP <- ll


	## PRIORS
	# prior on beta_naive
	LP <- LP + sum(dcauchy(beta_naive, 0, 2.5, log=TRUE))
	# prior on beta_surv
	LP <- LP + sum(dcauchy(beta_surv, 0, 2.5, log=TRUE))
	LP <- LP + dcauchy(qlogis(threshold), 0, 5, log=TRUE) # have to back transform the threshold


	list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rbinom(length(probs_naive), 1, probs_naive), parm=pars)
}


