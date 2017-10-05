integrated_model <- function(parm, Data, nested=FALSE)
{
	ll <- 0
	LP <- 0

	# naive model -- CHECK THIS AND SURVIVAL FOR HALF CAUCHY NOT FULL CAUCHY
	nm <- naive_model(parm, Data, TRUE)
	ll <- ll + nm$ll
	LP <- LP + nm$LP
	probs <- nm$probs

	#survival model
	sm <- survival_model(parm, Data, TRUE)
	ll <- ll + sm$ll
	LP <- LP + sm$LP
	probs_s <- sm$probs_s


	# population model
	pm <- population_model(parm, Data, TRUE)
	ll <- ll + pm$ll
	LP <- LP + pm$LP
	pop_r <- pm$r

	# integrate survival model with naive
	# compute a pseudo y dataset from the survival data using a threshold model
	# threshold will be set by choosing a value that maximizes TSS
	threshold <- 0.2
	stop("fix threshold")
	y_is <- as.integer(probs_s >= threshold)
	x_is <- Data$X_is
	beta_is <- parm[Data$pos.integrated_surv]
	probs_is <- plogis(x_is %*% beta_is)
	newll <- sum(dbinom(y_is, 1, probs_is, log = TRUE))
	ll <- ll + newll
	LP <- LP + newll # no additional parameters here to include


	# integrate population model
	# in this case the pseudo-y's come from the assumption that pop is present if r is > 0
	y_ip <- (pop_r > 0) * 1.0
	x_ip <- Data$x_ip
	beta_ip <- parm[Data$pos.integrated_pop]
	probs_ip <- plogis(x_ip %*% beta_ip)
	newll <- sum(dbinom(y_ip, 1, probs_ip, log=TRUE))
	ll <- ll + newll
	LP <- LP + newll

	# no priors on integrated betas, they've already been incorporated in the naive model
	list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rbinom(length(probs), 1, probs), parm=parm)
}