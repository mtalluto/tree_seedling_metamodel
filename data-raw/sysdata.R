naive_covariates <- as.matrix(read.csv("data-raw/mods_to_test.csv", row.names=1))
naive_covariates <- apply(naive_covariates, 1, function(x) names(which(x == 1)))
devtools::use_data(naive_covariates, overwrite = TRUE, internal = TRUE)
