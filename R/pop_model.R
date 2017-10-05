# needed for log sum exp
# devtools::install_github("rmcelreath/rethinking")


pop_model <- function(parm, Data, nested=FALSE)
{

	# this gets massively changed in my version
	# first; the mixture model is not necessary; if the species is absent, we cna't know anything about it's growth rate
	#	 so we just drop any rows where a species is absent; it is missing data
	# second; the old model was a bit funky with log lambdas and such
	#	here, we pre-compute the annual rate of increase r; this is our response variable, which we compute in a simple bayesian regression
	#	so a = alpha + x %*% beta; with a gaussian likelihood of a (in other words, a is a prediction of r1)
	#	time is not required in this model; it is pre-computed
	#	this vastly simplifies this model; for the integratoin we then predict rhat and decide whether r > 0; the integration sample
	#	size will be different for each species because of the differing number of obs (where the species was present)



	# unpacking data and params
	y <- Data$Y_pop
	x <- Data$X_pop
	time <- Data$time_pop
	beta <- parm[Data$pos.beta_pop]
	alpha <- parm[Data$pos.alpha_pop]
	theta <- parm[Data$pos.theta_pop]   ### mixing parameter
	sigma <- parm[Data$pos.sigma_pop]

	# model
	a <- alpha + x %*% beta

	## likelihood
	ll <- 0

	# the likihood for observations where y == 0; this is a slow computation because we have to use log_sum_exp to compute the likelihood for each individual point
	inds <- which(y == 0)
	ll <- ll + sum(sapply(inds, function(i) {
		rethinking::log_sum_exp(c(
			dbinom(1, 1, theta, log = TRUE),
			dbinom(0, 1, theta, log = TRUE) + dnorm(0, exp(a[i], sigma, log=TRUE)
		))
	}))

	# likelihood for all nonzero parts
	ll <- ll sum(dbinom(0, 1, theta, log=TRUE) + dnorm(y[y != 0], exp(a[y != 0], sigma, log = TRUE)))
	lp <- ll

	### priors on parameters
	lp <- lp + dgamma(sigma, 2, 0.1, log=TRUE) + LaplacesDemon::dhalfcauchy(beta, 2.5, log=TRUE) + 
		LaplacesDemon::dhalfcauchy(alpha, 2.5, log=TRUE) + dbeta(theta, 1, 1, log=TRUE)

	if(nested)
	{
		return(list(ll=ll, LP = lp, r = log(a) / time))
	} else {
		return(list(LP=lp, Dev=-2*ll, Monitor=lp, yhat=rnorm(length(a), a, sigma), parm=parm))
	}

}