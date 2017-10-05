# needed for log sum exp
# devtools::install_github("rmcelreath/rethinking")


# this model is a bit bizarre, because it is modelling the ratio of population sizes;
pop_model <- function(parm, Data, nested=FALSE)
{
	# unpacking data and params
	y <- Data$Y_pop
	x <- Data$X_pop
	time <- Data$time_pop
	beta <- parm[Data$pos.beta_pop]
	sigma <- parm[Data$pos.sigma_pop]

	# model
	a <- x %*% beta

	## likelihood
	ll <- 0

	# the likihood for observations where y == 0; this is a slow computation because we have to use log_sum_exp to compute the likelihood for each individual point
	inds <- which(y == 0)
	ll <- ll + sum(sapply(inds, function(i) {
		rethinking::log_sum_exp(c(
			dbinom(1, 1, theta_pop, log = TRUE),
			dbinom(0, 1, theta, log = TRUE) + dnorm(0, exp(a[i], sigma, log=TRUE)
		))
	}))

	# likelihood for all nonzero parts
	ll <- ll sum(dbinom(0, 1, theta, log=TRUE) + dnorm(y[y != 0], exp(a[y != 0], sigma, log = TRUE)))
	lp <- ll

	### priors on parameters
	lp <- lp + dgamma(sigma_pop, 2, 0.1, log=TRUE) + LaplacesDemon::dhalfcauchy(beta, 2.5, log=TRUE)

	if(nested)
	{
		return(list(ll=ll, LP = LP, r = log(a) / time);  stop("need to fix calculation of r")
	} else {
		return(list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rnorm(length(a), a, sigma), parm=parm))
	}

}