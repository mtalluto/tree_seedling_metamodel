library("treeSeedlingMetamodelData")
library(LaplacesDemon)
data(seedlings, envir = environment())
naive_covariates <- as.matrix(read.csv("data-raw/mods_to_test.csv", row.names=1))
naive_covariates <- apply(naive_covariates, 1, function(x) names(which(x == 1)))


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

prep_surv_covars <- function(dat, naiveDat, sp)
{

	D_surv <- dat[dat$species == sp,]
	X_surv <- data.frame(sizet = D_surv$x, dd0 = D_surv$dd0, tdiff = D_surv$tdiff)
	X_surv$sizet <- scale(X_surv$sizet)
	X_surv$dd0 <- scale(X_surv$dd0, center = attr(naiveDat, "scaled:center")['dd0'], scale = attr(naiveDat, "scaled:scale")['dd0'])
	X_surv$tdiff <- scale(X_surv$tdiff, center = attr(naiveDat, "scaled:center")['tdiff'], scale = attr(naiveDat, "scaled:scale")['tdiff'])
	as.matrix(cbind(1, X_surv))
}

prep_interactions <- function(X_s, X_n)
{
	X_int <- as.data.frame(X_s[,which(! colnames(X_s) %in% 'sizet')])
	if('dd0Xdd0' %in% colnames(X_n))
		X_int$dd0Xdd0 <- X_int$dd0 * X_int$dd0
	if('tdiffXtdiff' %in% colnames(X_n))
		X_int$tdiffXtdiff <- X_int$tdiff * X_int$tdiff
	if('tdiffXdd0' %in% colnames(X_n))
		X_int$tdiffXdd0 <- X_int$tdiff * X_int$dd0
	as.matrix(X_int)	
}


prep_naive_dat <- function(X,Y)
{
	X <- as.matrix(cbind(1, X))
	J <- ncol(X)

	parm.names <- LaplacesDemon::as.parm.names(list(beta_n=rep(0,J)))
	list(
		X = X,
		Y = Y,
		J = J,
		PGF = function(Data) rnorm(Data$J),
		parm.names = parm.names,
		pos.beta = grep("beta_n", parm.names),
		mon.names = c('LP')
		)
}

prep_int_dat <- function(X, Y, X_s, Y_s)
{
	
	dat <- prep_naive_dat(X,Y)	
	J_s <- ncol(X_s)

	parm.names <- c(dat$parm.names, LaplacesDemon::as.parm.names(list(beta_s=rep(0,J_s))))
	# parm.names <- LaplacesDemon::as.parm.names(list(beta_n=rep(0,J), beta_s=rep(0,J_s), threshold = 0))
	# create an integration x dataset with the same cols as survival, but with interactions we know about to be used if needed
	X_i <- prep_interactions(X_s, dat$X)
	int.indices <- c(1, which(colnames(dat$X) %in% colnames(X_i))) #always add the intercept
	int.names <- paste0("beta_n[", int.indices, "]")
	
	# add integrated items to the data structure
	dat$X_s <- X_s
	dat$Y_s <- Y_s
	dat$J_s <- J_s
	dat$X_i <- X_i
	dat$PGF <- function(Data) rnorm(Data$J + Data$J_s)
	dat$parm.names <- parm.names
	dat$pos.beta = grep("beta_n", parm.names)
	dat$pos.beta_surv = grep("beta_s", parm.names)
	# dat$pos.threshold = grep("threshold", parm.names)
	dat$pos.integrated = which(parm.names %in% int.names)
	dat$mon.names = c('LP')
	dat
}


naive_model <- function(parm, Data, nested=FALSE)
{
	ll <- 0

	beta <- parm[Data$pos.beta]
	x <- Data$X
	y <- Data$Y
	probs <- plogis(x %*% beta)
	## bernoulli likelihood: p where y == 1; 1-p where y == 0
	ll <- ll + sum(dbinom(y, 1, probs, log=TRUE))

	LP <- ll + sum(dcauchy(beta, 0, 2.5, log=TRUE))
	stop("half cauchy!")

	if(nested)
	{
		return(list(ll=ll, LP = LP, probs = probs))
	} else {
		return(list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rbinom(length(probs), 1, probs), parm=parm))
	}
}


survival_model <- function(parm, Data, nested=FALSE)
{
	ll <- 0

	beta_surv <- parm[Data$pos.beta_surv]
	x_s <- Data$X_s
	y_s <- Data$Y_s
	probs_s <- plogis(x_s %*% beta_surv)
	## bernoulli likelihood: p where y == 1; 1-p where y == 0
	ll <- ll + sum(dbinom(y_s, 1, probs_s, log=TRUE))
	LP <- ll + sum(dcauchy(beta_surv, 0, 2.5, log=TRUE))

	if(nested)
	{
		return(list(ll=ll, LP = LP, probs_s = probs_s))
	} else {
		return(list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rbinom(length(probs_s), 1, probs_s), parm=parm))
	}

}


int_model <- function(parm, Data)
{
	ll <- 0
	LP <- 0

	# naive model
	nm <- naive_model(parm, Data, TRUE)
	ll <- ll + nm$ll
	LP <- LP + nm$LP
	probs <- nm$probs

	#survival model
	sm <- survival_model(parm, Data, TRUE)
	ll <- ll + sm$ll
	LP <- LP + sm$LP
	probs_s <- sm$probs_s


	# # integrated part of model
	# # compute a pseudo y dataset from the survival data using a threshold model
	# # the threshold is a parameter that will be estimated from data
	threshold <- 0.2
	# # threshold <- plogis(parm[Data$pos.threshold])
	y_i <- as.integer(probs_s >= threshold)
	x_i <- Data$X_i
	beta_i <- parm[Data$pos.integrated]
	probs_i <- plogis(x_i %*% beta_i)
	ll <- ll + sum(dbinom(y_i, 1, probs_i, log = TRUE))
	LP <- LP + ll

	# # prior on betas
	# # using a cauchy 0,2.5
	# LP <- LP + dcauchy(qlogis(threshold), 0, 5, log=TRUE) # have to back transform the threshold

	list(LP=LP, Dev=-2*ll, Monitor=LP, yhat=rbinom(length(probs), 1, probs), parm=parm)
}




spNames <- names(naive_covariates)
dat_n <- seedlings$fia_adults
dat_s <- seedlings$survival
X_all <- prep_naive_covars(dat_n, naive_covariates)
Y_all <- dat_n[,spNames] 

#### START LOOP OVER SPECIES
sp <- spNames[1]
# for(sp in spNames)
# {
	X <- X_all[[sp]]
	Y <- Y_all[,sp]
	X_s <- prep_surv_covars(dat_s, X_all, sp)
	Y_s <- dat_s$surv[dat_s$species == sp]

	ldDat <- prep_int_dat(X, Y, X_s, Y_s)
	naiveDat <- prep_naive_dat(X,Y)

	# ldMod <- LaplacesDemon::LaplaceApproximation(int_model, LaplacesDemon::GIV(int_model, ldDat, PGF=TRUE),
	# 		ldDat, Iterations = 10)

	# ldMod <- LaplacesDemon::LaplacesDemon(int_model, ldDat, LaplacesDemon::GIV(int_model, ldDat, PGF=TRUE), 
	# 	Algorithm = "AMWG", specs=list(B = NULL, n=0, Periodicity = 50), Iterations = 10000)

## final algo used for abla
	# Initial.Values <- as.initial.values(ldMod)
	# ldMod <- LaplacesDemon(int_model, Data=ldDat, Initial.Values,
 #     Covar=ldMod$Covar, Iterations=140000, Status=1168, Thinning=140,
 #     Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100,
 #     n=0, w=1))
	### note that because A is inf, these aren;t true final samples - need to set a = 0 for that
	# saveRDS(ldMod, paste0('res/int_', sp, '.rds'))

# 	naiveMod <- LaplacesDemon::LaplacesDemon(naive_model, naiveDat, LaplacesDemon::GIV(naive_model, naiveDat, PGF=TRUE), 
# 		Algorithm = "AMWG", specs=list(B = NULL, n=0, Periodicity = 50), Iterations = 10000)

# saveRDS(naiveMod, 'res/naive_initial_abla.rds')

naiveMod <- readRDS('res/naive_initial_abla.rds')
Initial.Values <- as.initial.values(naiveMod)
naiveMod <- LaplacesDemon(naive_model, Data=naiveDat, Initial.Values,
     Covar=naiveMod$Covar, Iterations=100000, Status=1168, Thinning=100,
     Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100,
     n=0, w=1))

saveRDS(naiveMod, "res/naive_abla.rds")
# }







### quick and dirty response curve
new_dd0 <- seq(-1.4, 4.5, 0.1)
new_tdiff <- seq(-3, 3.5, 0.1)
# newdat <- data.frame(i=1, dd0=new_dd0, tdiff=0, gsp=0, winp=0, dd0Xdd0 = new_dd0^2, winpXdd0=0, gspXtdiff=0, winpXgsp=0, winpXwinp=0)
newdat <- data.frame(i=1, dd0=0, tdiff=new_tdiff, gsp=0, winp=0, dd0Xdd0 = 0, winpXdd0=0, gspXtdiff=0, winpXgsp=0, winpXwinp=0)
newdat = as.matrix(newdat)
preds_naive <- newdat %*% t(naiveMod$Posterior2)
preds_integrated <- newdat %*% t(intMod$Posterior2[,1:10])

plot(new_tdiff, plogis(rowMeans(preds_naive)), xlab='tdiff', ylab = 'presence', lwd=2, col='red', type='l')
lines(new_tdiff, plogis(rowMeans(preds_integrated)), col='blue')