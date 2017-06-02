.onLoad <- function(libname, pkgname){
	library(rstan)
	library("treeSeedlingMetamodelData")
	data(seedlings, envir = environment())
}
