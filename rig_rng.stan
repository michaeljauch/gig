functions {
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
   
   real rgig32_rng(real a, real b){
     //if(uniform_rng(0,1) < 1/(1 + 1/sqrt(a*b))){
       if(uniform_rng(0,1) > 1.0/(1 + sqrt(a*b))){
       return(rig_single_rng(sqrt(b/a), b) + exponential_rng(.5*a));
       } else{
         return(1.0/rig_single_rng(sqrt(b/a), b) + exponential_rng(.5*a));
  }
   }
  
   // vector rig_multiple_rng(int n, real mu, real lambda){
   //   vector[n] output; 
   //   real nu; 
   //   real x; 
   //   real y;
   //   real z; 
   //   for(i in 1:n){
   //     output[i] = rig_single_rng(mu, lambda);
   //   }
   //   return(output);
   // }
   // 
   // vector gig_gibbs_rng(int n, real p, real a, real b){
   //   real x; 
   //   real y; 
   //   vector[n] x_vec; 
   //   // initial values
   //   x = 1; 
   //   if(p < -.5){
   //     for(i in 1:n){
   //     y = gamma_rng(-(p+.5),x);
   //     x = rig_single_rng(sqrt(b/(a+2*y)), b);
   //     x_vec[i] = x;
   //     }
   //     return(x_vec);
   //   } else{
   //     for(i in 1:n){
   //       y = gamma_rng(p+.5,1/x);
   //       x = rig_single_rng(sqrt((b+2*y)/a), b+2*y);
   //       x_vec[i] = x;
   //     }
   //     return(x_vec);
   //   }
   // }
   
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
