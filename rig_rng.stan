functions {
  // This function produces a single draw from the inverse Gaussian distribution 
  // via the algorithm of Michael et al. (1978). This implementation is a slight 
  // modification of the algorith Wikipedia page for the inverse 
  // Gaussian distribution. The changes were made in an effort to reduce the 
  // number of arithmetic operations. 
  real rig_single_rng(real mu, real lambda){
     real x; 
     real y;
     //nu = normal_rng(0,1);
     y = normal_rng(0,1)^2;
     //x = mu + mu^2*y/(2*lambda) - 
     //mu/(2*lambda)*sqrt(4*mu*lambda*y + mu^2*y^2);
     x = mu*(1 + mu*y/(2*lambda) - 
     1.0/(2*lambda)*sqrt(4*mu*lambda*y + mu^2*y^2));
     if(uniform_rng(0,1) <= mu/(mu+x)){
       return(x); 
     } else{
       return(mu^2/x);
     }
   }
   
   // This function uses the identity in Proposition 4 to simulate from the 
   // GIG in the case when p=3/2. The mixture weight w simplifies in that case. 
   real rgig32_rng(real a, real b){
     //if(uniform_rng(0,1) < 1/(1 + 1/sqrt(a*b))){
       if(uniform_rng(0,1) > 1.0/(1 + sqrt(a*b))){
       return(rig_single_rng(sqrt(b/a), b) + exponential_rng(.5*a));
       } else{
         return(1.0/rig_single_rng(sqrt(b/a), b) + exponential_rng(.5*a));
  }
   }
   
   // This function performs the update of sigma2 in the data-augmented 
   // Gibbs sampler described in the example at the end of Section 2. 
   real gibbs_update_rng(real p, real a, real b, real x_prev){
     //real x; 
     //real y;
     real b2y; 
     if(p < -.5){
       //y = gamma_rng(-(p+.5),x_prev);
       //x = rig_single_rng(sqrt(b/(a+2*y)), b);
       return(rig_single_rng(sqrt(b/(a+2*gamma_rng(-(p+.5),x_prev))), b));
     } else{
       //y = gamma_rng(p+.5,1/x_prev);
       b2y = b+2*gamma_rng(p+.5,1/x_prev);
       //x = rig_single_rng(sqrt((b+2*y)/a), b+2*y);
       return(rig_single_rng(sqrt(b2y/a), b2y));
     }
   }
   
}
