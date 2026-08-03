library(Rcpp)
library(RcppArmadillo)
sourceCpp("source/binary_power.cpp")
binary_analysis <- function(R_0, delta, n_1, k, r, a_0,  
                            Delta = 0, gridsize = 100, alpha_0 = 0.05){
  
  settings <- 
    expand.grid(
      R_0 = R_0, delta = delta, n_1 = n_1, k = k, r = r, a_0 = a_0
    )
  
  out <-
    binary_power_vectorized(
      R_0_list = settings$R_0, 
      delta_list  = settings$delta, 
      n_1_list = settings$n_1,
      n_0_list  = with(settings, floor(n_1 /k)),
      r_list = settings$r,
      a_0_list = settings$a_0, 
      alpha_0 = alpha_0, 
      gridsize = gridsize, 
      Delta = Delta
    )
  
  results <- cbind(settings, power = out[["power"]], type_1_error = out[["type_1_error"]])
  
  results
  
}

