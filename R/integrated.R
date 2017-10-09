#' Utility function to find the index position in a LD parameter object of a named parameter
param_index <- function(dat, parm) 
{
	grep(parm, dat$parm.names)
}

#' Compute the true skill statistic
#' 
#' @param thres The threshold to apply if probabilities (pr) are specified
#' @param y Data points; observed 1s and 0s
#' @param pr Predicted probabilities
tss <- function(thresh, y, pr)
{
	yhat <- as.integer(pr > thresh)
	
	a <- sum(yhat & y)
	b <- sum(yhat & !y)
	c <- sum(!yhat & y)
	d <- sum(!yhat & !y)
	sens <- a / (a+c)
	spec <- d / (b+d)
	sens + spec - 1
}


#' Integrated data
#' Create data list for laplaces demon for integrated model
#' @param ndat Data frame for the naive model
#' @param sdat Data frame for the survival model
#' @param pdat Data frame for population model
#' @export
integrated_ld_dat <- function(ndat, sdat, pdat)
{
	naive_data <- naive_ld_dat(ndat)
	surv_data <- survival_ld_dat(sdat)
	pop_data <- pop_ld_dat(pdat)
	parm.names <- c(naive_data$parm.names, surv_data$parm.names, pop_data$parm.names)
	naive_data$parm.names <- surv_data$parm.names <- pop_data$parm.names <- parm.names 

	# map integrated parameters onto existing parameters
	indices_ns <- which(colnames(naive_data$X_n) %in% colnames(sdat))
	names_ns <- paste0("beta_n[", indices_ns, "]")
	indices_np <- which(colnames(naive_data$X_n) %in% colnames(pdat))
	names_np <- paste0("beta_n[", indices_np, "]")

	list(
		nData = naive_data,
		sData = surv_data,
		pData = pop_data,
		X_ns = as.matrix(sdat[,colnames(naive_data$X_n)[indices_ns]]), # taking care to preserve ordering of columns
		X_np = as.matrix(pdat[,colnames(naive_data$X_n)[indices_np]]), # taking care to preserve ordering of columns
		J = naive_data$J + surv_data$J + pop_data$J,
		N = naive_data$N + surv_data$N + pop_data$N,
		PGF = function(Data) rnorm(Data$J + 4), # add one for each intercept and one for population sigma
		parm.names = parm.names,
		mon.names = c('LP'),
		params_ns = names_ns,
		params_np = names_np
	)
}



#' Integrated model
#' Log probability integrating naive and survival and population models
#' 
#' @param parm Parameter list
#' @param Data Data list
#' 
#' @export
integrated_lp <- function(parm, Data)
{
	ll <- 0
	LP <- 0
	alpha_n <- parm[param_index(Data, 'alpha_n')]

	# naive model
	nm <- naive_lp(parm, Data$nData, TRUE)
	ll <- ll + nm$ll
	LP <- LP + nm$LP
	probs_n <- nm$probs

	#survival model
	sm <- survival_lp(parm, Data$sData, TRUE)
	ll <- ll + sm$ll
	LP <- LP + sm$LP
	probs_s <- sm$probs

	# population model
	pm <- population_lp(parm, Data$pData, TRUE)
	ll <- ll + pm$ll
	LP <- LP + pm$LP
	pop_r <- pm$r


	# compute survival threshold using TSS
	threshold <- optimize(tss, interval = c(0, 1), maximum = TRUE, y = Data$sData$Y_s, pr = probs_s)
	threshold <- threshold$maximum

	# data and parameters for integrated model
	y_ns <- as.integer(probs_s >= threshold)
	x_ns <- Data$X_ns
	beta_ns <- parm[which(Data$parm.names %in% Data$params_ns)]

	# integrate population model
	# in this case the pseudo-y's come from the assumption that pop is present if r is > 0
	y_np <- (pop_r > 0) * 1.0
	x_np <- Data$X_np
	beta_np <- parm[which(Data$parm.names %in% Data$params_np)]

	# integrated model likelihood
	probs_ns <- plogis(alpha_n + x_ns %*% beta_ns) # always use naive intercept
	ll_ns <- sum(dbinom(y_ns, 1, probs_ns, log = TRUE))
	probs_np <- plogis(alpha_n + x_np %*% beta_np)
	ll_np <- sum(dbinom(y_np, 1, probs_np, log = TRUE))
	ll <- ll + ll_ns + ll_np
	LP <- LP + ll_ns + ll_np # no new parameters, so no additional priors

	list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rbinom(length(probs), 1, probs), parm=parm)
}



