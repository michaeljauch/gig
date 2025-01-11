####################################
# Benchmarking Posterior Inference # 
####################################

library(rbenchmark)
library(rstan)
expose_stan_functions("~/Documents/gig/rig_rng.stan", show_compiler_warnings = TRUE)
library(GIGrvg) # For Hormann & Leydold method 
library(boodist) # For Devroye method
library(coda) # For effectiveSize function

# Simulate the data 
n <- 100
mu_true <- 1
sigma2_true <- 1

# Define the prior hyperparameters 
theta0 <- 0 
tau02 <- 1e2
p0 <- 3/4
a0 <- 1
b0 <- 1

# Plot the GIG prior for sigma^2
# xgrid <- seq(0, 5, length.out=1000)
# plot(x=xgrid, y=dgig(xgrid, lambda=p0, chi=b0, psi=a0), type="l")

benchmark("HormannLeydold" = {
  
  # Simulate data 
  y <- rnorm(n, mu_true, sd=sqrt(sigma2_true))
  ybar <- mean(y)
  
  # Gibbs sampler
  num_iter <- 5000
  draws <- matrix(rep(NA, 2*num_iter), ncol=2)
  sigma2 <- 1 # Initial value for the Markov Chain 
  for(iter in 1:num_iter){
    
    # Sample mu 
    taun2 <- 1/(n/sigma2 + 1/tau02)
    mu <- rnorm(n=1, mean=taun2*(n*ybar/sigma2 + theta0/tau02), sd=sqrt(taun2))
    
    # Sample sigma2 (lambda=p, psi=a, chi=b)
    bn <- b0 + sum((y-mu)^2)
    sigma2 <- rgig(n=1, lambda=p0-n/2, psi=a0, chi=bn)
    
    # Store the results
    draws[iter,] <- c(mu, sigma2) 
    
  }
},
"Devroye" = {
  
  # Simulate data 
  y <- rnorm(n, mu_true, sd=sqrt(sigma2_true))
  ybar <- mean(y)
  
  # Gibbs sampler
  num_iter <- 5000
  draws <- matrix(rep(NA, 2*num_iter), ncol=2)
  sigma2 <- 1 # Initial value for the Markov Chain
  for(iter in 1:num_iter){
    
    # Sample mu 
    taun2 <- 1/(n/sigma2 + 1/tau02)
    mu <- rnorm(n=1, mean=taun2*(n*ybar/sigma2 + theta0/tau02), sd=sqrt(taun2))
    
    # Sample sigma2
    bn <- b0 + sum((y-mu)^2)
    sigma2 <- GeneralizedInverseGaussian$new(theta = sqrt(a0*(b0 + sum((y-mu)^2))), 
                                             eta = sqrt((b0 + sum((y-mu)^2))), 
                                             lambda=p0-n/2)$r(1)
    
    # Store the results
    draws[iter,] <- c(mu, sigma2) 
    
  }
},
"ProposedMethod" = {
  
  # Simulate data 
  y <- rnorm(n, mu_true, sd=sqrt(sigma2_true))
  ybar <- mean(y)
  
  # Gibbs sampler
  num_iter <- 5000
  draws <- matrix(rep(NA, 2*num_iter), ncol=2)
  sigma2 <- 1 # Initial value for the Markov Chain
  for(iter in 1:num_iter){
    
    # Sample mu 
    taun2 <- 1/(n/sigma2 + 1/tau02)
    mu <- rnorm(n=1, mean=taun2*(n*ybar/sigma2 + theta0/tau02), sd=sqrt(taun2))
    
    # Sample omega and sigma2 simultaneously: 
    bn <- b0 + sum((y-mu)^2)
    sigma2 <- gibbs_update_rng(p=p0-n/2, a=a0, b=bn, x_prev=sigma2)
    
    # Store the results
    draws[iter,] <- c(mu, sigma2) 
    
  }
}, 
replications = 1000,
columns = c("test", "replications", "elapsed",
            "relative", "user.self", "sys.self"))


####################################
# Evaluating Effective Sample Size #  
####################################

num_repetitions <- 100 # Number of times we run each Gibbs sampler
num_iter <- 50000 # Length of each Markov chain 
warmup <- 5000 # Number of warmup iterations 

# Create a matrix where we store the ESS values 
# The first column stores the ESS value for mu with the DA Gibbs sampler. 
# The second column stores the ESS value for sigma^2 with the DA Gibbs sampler. 
# The third column stores the ESS value for mu with the original Gibbs sampler. 
# The fourth column stores the ESS value for sigma^2 with the original Gibbs 
# sampler. 
# The ESS values should not differ depending on which package (GIGrvg or 
# boodist) we use, so we only compare with GIGrvg. 
ESS_mat <- matrix(rep(NA, num_repetitions*4), nrow=num_repetitions)

for(rep in 1:num_repetitions){
  
  # Simulate data 
  y <- rnorm(n, mu_true, sd=sqrt(sigma2_true))
  ybar <- mean(y)
  
  
  ### Data-augmented Gibbs sampler 
  draws <- matrix(rep(NA, 2*num_iter), ncol=2)
  sigma2 <- 1 # Initial value for the Markov Chain 
  for(iter in 1:num_iter){
    
    # Sample mu 
    taun2 <- 1/(n/sigma2 + 1/tau02)
    mu <- rnorm(n=1, mean=taun2*(n*ybar/sigma2 + theta0/tau02), sd=sqrt(taun2))
    
    
    # Sample omega and sigma2 simultaneously: 
    bn <- b0 + sum((y-mu)^2)
    sigma2 <- gibbs_update_rng(p=p0-n/2, a=a0, b=bn, x_prev=sigma2)
    
    # Store the results
    draws[iter,] <- c(mu, sigma2) 
  }
  
  ESS_mat[rep, 1:2] <- effectiveSize(draws[(warmup+1):num_iter,])
  
  ### Hormann-Leydold 
  
  draws <- matrix(rep(NA, 2*num_iter), ncol=2)
  sigma2 <- 1 # Initial value for the Markov Chain 
  for(iter in 1:num_iter){
    
    # Sample mu 
    taun2 <- 1/(n/sigma2 + 1/tau02)
    mu <- rnorm(n=1, mean=taun2*(n*ybar/sigma2 + theta0/tau02), sd=sqrt(taun2))
    
    # Sample sigma2 (lambda=p, psi=a, chi=b)
    bn <- b0 + sum((y-mu)^2)
    sigma2 <- rgig(n=1, lambda=p0-n/2, psi=a0, chi=bn)
    
    # Store the results
    draws[iter,] <- c(mu, sigma2) 
    
  }
  
  ESS_mat[rep, 3:4] <- effectiveSize(draws[(warmup+1):num_iter,])
  
}

# Median ESS/iter for mu (DA Gibbs sampler)
median(ESS_mat[,1]/(num_iter-warmup))

# Median ESS/iter for sigma^2 (DA Gibbs sampler)
median(ESS_mat[,2]/(num_iter-warmup))

# Median ESS/iter for mu (original Gibbs sampler)
median(ESS_mat[,3]/(num_iter-warmup))

# Median ESS/iter for sigma^2 (original Gibbs sampler)
median(ESS_mat[,4]/(num_iter-warmup))
