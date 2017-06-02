//Naive Distribution Model
//Logistic Regression
//Seedling Metamodel
//Created by: Paige E. Copenhaver-Parry
//Created on: 1/28/2017
//Last modified: 3/28/2017

data{
	int<lower=0> N; //number of observations
	int<lower=0> K; //number of predictors
	matrix[N,K] x; //predictor matrix
	int<lower=0,upper=1> y[N]; //binary response 

}

parameters{
	vector[K] beta; //regression parameters 
	real<lower=0> tau; //variance on betas 
	real<lower=0> sigma; //standard deviation on y regression
	real alpha; //intercept
}



model{
	//priors
	tau~cauchy(0,2.5);
	for (k in 1:K){
	beta[k]~normal(0,tau);
	}
	
	//likelihood
	y~bernoulli_logit(alpha + x*beta);
}

generated quantities{
	vector[N] odds;
	vector[N] prob;
	for(n in 1:N){
		odds[n]=exp(alpha + x[n]*beta);
		prob[n]=odds[n]/(odds[n]+1);
	}
}
