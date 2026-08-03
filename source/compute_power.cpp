// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


// This is a working version of skewnorm that will include different features
// Features checklist:
// Censoring
// Zero inflation


//Notes:
// * the plugins thing gives you a suite of nice functions from arma, such
// as normcdf. These are more flexible in that they accept arma objects.
// * Rcpp eigen has its own suite of functions and objects. Might be slightly
// faster for some things.
// * Rcpp operations should be automatically vectorized thanks to ``sugar". 
// * Standard functions requiring Rcpp NumericVector-s or NumericMatrix-s can be
// accessed as Rcpp::pnorm(), for example
// * Non-vectorized functions can be used supplying doubles to functions such as
// R::pnorm(). 
// * See https://statisticsglobe.com/convert-vector-from-rcpparmadillo-to-rcpp-r for
// converting between arma and rcpp objects
// * For now, RcppDist and armadillo extensions are probably not needed. 
// RcppArmadilloExtensions/sample.h just gives you the RcppArmadillo::sample function.
// Easy to look up. If using NumericVectors, Rcpp::sample should work fine.
// * RcppDist.h looks like it pulls in RcppArmadillo.h automatically, so you may
// not need the second one if you have the first one. (Yes! Confirmed)
// * RcppDist.h has the truncated normal built in.
// std gives standard functions

// DP


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//Note: notation for arma objects (e.g., arma::mat) and Rcpp objects (e.g., 
// NumericVector) are fairly different.

double cmin(double a, double b){
  
  if (a <= b){
    
    return(a);
  } else{
    return(b);
  }
  
}

double diff_prob(double alpha_n, double beta_n, double mu_0, double y_1, double mu_0n, 
                 int n_1, int n_0, int n_0h, double dx, double limit){
  
  double nu = 2 * alpha_n;
  double A_11 = alpha_n * (n_0h + n_0) / beta_n;
  double A_22 = alpha_n * n_1 / beta_n;
  double post_mu_0_sigma = pow(beta_n / alpha_n / (n_0 +  n_0h), 0.5);
  double F_t_mean = y_1;
  double F_t_sigma = pow((nu + A_11 * pow(mu_0 - mu_0n, 2.0)) / A_22 / (nu + 1), 0.5);
  double out = 0;
  double fun;
  
  for (double x = -limit; x < limit; x = x + dx){
    
    fun = 
      R::pt((x - F_t_mean) / F_t_sigma, nu + 1, 1, 0) * 
      R::dt((x - mu_0n) / post_mu_0_sigma, nu, 0) / post_mu_0_sigma;
    
    out = out + dx * fun; 
  }
  
  return(1-out);
}




double power_fixed(double mu_0, double mu_1, double lambda, 
                   double n_1, double n_0, double n_h, double y_h, 
                   double alpha, double beta){
  
  double y_0_upper = R::qnorm5(1 - 0.5 * pow(trunc,1.0/3.0), mu_0, 
                                   pow(lambda * n_0, -1.0), 1, 0);
  double y_0_lower = 2 * y_0 - y_0_upper;
  
  double y_1_upper = R::qnorm5(1 - 0.5 * pow(trunc,1.0/3.0), mu_1, 
                               pow(lambda * n_0, -1.0), 1, 0);
  double y_1_lower = 2 * y_1 - y_1_upper;
  
  double S2_upper = R::qchisq(1 - pow(trunc, 1.0/3.0), n_1 + n_0 - 2.0, 1, 0) / lambda;
  
  
  
  
}














