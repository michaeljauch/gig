####################################
# Benchmarking Posterior Inference # 
####################################

set.seed(1)

expose_stan_functions("~/Desktop/rig_rng.stan", show_compiler_warnings = TRUE)
library(GIGrvg) # For Hormann & Leydold method 
library(boodist) # For Devroye method

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
replications = 250,
columns = c("test", "replications", "elapsed",
            "relative", "user.self", "sys.self"))


####################################
# Evaluating Effective Sample Size #  
####################################

# Simulate data 
y <- rnorm(n, mu_true, sd=sqrt(sigma2_true))
ybar <- mean(y)

num_iter <- 100000

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

effectiveSize(draws[10001:100000,])

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

effectiveSize(draws[10001:100000,])


