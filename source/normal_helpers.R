library(Rcpp)
library(RcppArmadillo)
sourceCpp("source/normal_power.cpp")

fast_sim_numerical_wrapper <-
  function(mu_0, alpha, beta, n_1, k, r, a_0, delta, alpha_0 = 0.05,
           gridsize = 100, tail_error = 0.0001){
    
    settings <-
      expand.grid(mu_0, alpha, beta, n_1, k, r, a_0, delta)
    colnames(settings) <- c("mu_0", "alpha", "beta", "n_1", "k", "r", "a_0", 
                            "delta")
    
    out <- fast_sim_numerical(
      mu_0_list = settings$mu_0, 
      delta_list = settings$delta,
      n_1_list = settings$n_1,
      n_0_list = floor(settings$n_1 / settings$k),
      a_0_list = settings$a_0,
      r_list = settings$r,
      alpha_list  = settings$alpha,
      beta_list = settings$beta,
      alpha_0 = alpha_0, 
      gridsize = gridsize,
      tail_error = tail_error
    )
    
    results <- cbind(settings, power = out[["power"]], type_1_error = out[["type_1_error"]])
    
    list(
      results = results
    )
  }

fast_sim_numerical_wrapper_flex <-
  function(mu_0, alpha, beta, n_1, k, r, a_0, delta, alpha_0 = 0.05,
           gridsize = 100, tail_error = 0.0001){
    
    settings <-
      data.frame(mu_0, alpha, beta, n_1, k, r, a_0, delta)
    colnames(settings) <- c("mu_0", "alpha", "beta", "n_1", "k", "r", "a_0", 
                            "delta")
    
    out <- fast_sim_numerical(
      mu_0_list = settings$mu_0, 
      delta_list = settings$delta,
      n_1_list = settings$n_1,
      n_0_list = floor(settings$n_1 / settings$k),
      a_0_list = settings$a_0,
      r_list = settings$r,
      alpha_list  = settings$alpha,
      beta_list = settings$beta,
      alpha_0 = alpha_0, 
      gridsize = gridsize,
      tail_error = tail_error
    )
    
    results <- cbind(settings, power = out[["power"]], type_1_error = out[["type_1_error"]])
    
    list(
      results = results
    )
  }
