###############################################
# Benchmark exact simulation method for p=3/2 #  
###############################################

set.seed(1)

expose_stan_functions("~/Desktop/rig_rng.stan", show_compiler_warnings = TRUE)
library(GIGrvg) # For Hormann & Leydold method 
library(boodist) # For Devroye method

# Instead of producing a vector of GIG values, the stan f'n rgig32_rng produces 
# one GIG at a time. There are two reasons for doing it this way: 
# - It'd be easy to implement the vectorization badly and then have the entire 
# difference between the algorithms come down to that. 
# - In the paper, we discuss using the algorithms within a Gibbs sampler. This 
# reflects that "varying parameter case" (language from Hormann and Leydold). 

p <- 3/2; a <- 1; b <- 1; n <- 1

benchmark("HormannLeydold" = {
  x <- rep(0, n)
  for(i in 1:n){
    x[i] <- rgig(n=1, lambda=p, chi=b, psi=a)
  }
},
"Devroye" = {
  x <- rep(0, n)
  for(i in 1:n){
    x[i] <- GeneralizedInverseGaussian$new(theta = sqrt(a*b), 
                                           eta = sqrt(b/a), 
                                           lambda=p)$r(1)
  }
},
"ProposedMethod" = {
  x <- rep(0, n)
  for(i in 1:n){
    x[i] <- rgig32_rng(a=a, b=b)
  }
}, 
replications = 1000,
columns = c("test", "replications", "elapsed",
            "relative", "user.self", "sys.self"))