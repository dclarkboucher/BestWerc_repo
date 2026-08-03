// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

// [[Rcpp::export]]
double beta_diff_prob(double Delta, double alpha0, double beta0, double alpha1, double beta1,
                 double gridsize, double tail_error){


  double p0_lower = R::qbeta(tail_error/2.0,
                             alpha0, beta0, 1,0);
  double p0_upper = R::qbeta(1-tail_error/2.0,
                             alpha0, beta0, 1,0);
  double dx0 = (p0_upper - p0_lower)/gridsize;
  double prob = 0;
  for (double p0 = p0_lower; p0 <= p0_upper; p0 = p0 + dx0){

    prob = prob + dx0 * R::dbeta(p0, alpha0, beta0, 0) * 
      R::pbeta(p0 + Delta,alpha1,beta1,0,0); // Right tail
  }

  return(prob);

}

// [[Rcpp::export]]
double binary_power_exact(double alpha, double beta, 
                          int n_1, int n_0, double R_0, double R_1, 
                          double Delta = 0, double alpha_0 = 0.05, double gridsize = 100
){
  
  double power = 0;
  double a_0;
  double b_0;
  double prob;
  
  for (int x_0 = 0; x_0 <= n_0; x_0 ++){
    
    a_0 = alpha + x_0;
    b_0 = beta + n_0 - x_0;
    
    for (int x_1 = 0; x_1 <= n_1; x_1++){

      prob = beta_diff_prob(Delta, a_0, b_0, x_1 + 1, n_1 - x_1 + 1,
                            gridsize, 1e-6);
      if (prob >= 1-alpha_0){
              power = power + R::dbinom(x_0, n_0, R_0, 0) * R::dbinom(x_1, n_1, R_1, 0);
      }
    }
    
  }
  
  return(power);
  
}

// Calculates power by simulation instead of integration 
// double binary_power_random(double alpha, double beta, 
//                            int n_1, int n_0, double R_0, double R_1, 
//                            double Delta = 0, double alpha_0 = 0.05, double gridsize = 100, int reps = 1000
// ){
//   
//   double power = 0;
//   double a_0;
//   double b_0;
//   double a_1;
//   double b_1;
//   double x_0;
//   double x_1;
//   double prob;
//   
//   for (int rep = 0; rep < reps; rep++){
//     
//     set_seed(rep);
//     x_0 = R::rbinom(n_0, R_0);
//     x_1 = R::rbinom(n_1, R_1);
//   
//     a_0 = alpha + x_0;
//     b_0 = beta + n_0 - x_0;
//     a_1 = 1 + x_1;
//     b_1 = 1 + n_1 - x_1;
//     
//     prob = beta_diff_prob(Delta, a_0, b_0, a_1, b_1, gridsize, 1e-7);
//     if (prob >= 1-alpha_0) power = power + 1;
//     
//     
//   }
// 
//   return(power / double(reps));
//   
//   
// }
  
// [[Rcpp::export]]
List binary_power_vectorized(
    arma::vec R_0_list,
    arma::vec delta_list,
    arma::vec n_1_list,
    arma::vec n_0_list,
    arma::vec a_0_list,
    arma::vec r_list,
    double alpha_0,
    double gridsize = 100,
    double Delta = 0
){
  
  int size = R_0_list.n_elem;
  double R_0;
  double delta;
  int n_0;
  int n_1;
  double n_0h;
  double a_0;
  double r;
  double alpha;
  double beta;
  double r_0h;
  
  // Storage
  arma::vec power_list(size);
  arma::vec t1e_list(size);
  
  // Loop
  for (int i = 0; i < size; i++){
  
    R_0 = R_0_list(i);
    delta = delta_list(i);
    n_0 = n_0_list(i);
    n_1 = n_1_list(i);
    a_0 = a_0_list(i);
    r = r_list(i);
    n_0h = a_0 / (1 - a_0) * double(n_0);
    r_0h = r * delta + R_0;
    alpha = 1 + n_0h * r_0h;
    beta = 1 + n_0h * (1 - r_0h);
    
    power_list(i) =
      binary_power_exact(
        alpha, beta, n_1, n_0, R_0, R_0 + delta,Delta,
        alpha_0, gridsize
      );
    t1e_list(i) =
      binary_power_exact(
        alpha, beta, n_1, n_0, R_0, R_0,Delta,
        alpha_0, gridsize
      );
  }
  
  List out = List::create(Named("power") = power_list, 
                          Named("type_1_error") = t1e_list);
  return(out);
}


// [[Rcpp::export]]
double find_R1(double power, double R0, int n_1, int n_0, double Delta = 0, double alpha_0 = 0.05,
               double gridsize = 100, double err = 1e-5){
  
  double max = 1 - 1e-8;
  double min = R0;
  double R1 = 0.5 * (max + min);
  double pow;
  
  // Check maximum first.
  pow = 
    binary_power_exact(
      1, 1, n_1, n_0,  R0, max, Delta,
      alpha_0, gridsize
    );
  if (pow < power){
    return(1.1);
  }
  
  // Binary search
  while(max - min > err){
    
    pow = 
      binary_power_exact(
        1, 1, n_1, n_0, R0, R1,Delta,  
        alpha_0, gridsize
      );
    
    if (pow >= power){
      max = R1;
    } else {
      min = R1;
    }
    
    R1 = 0.5 * (min+max);
      
  }
  
  return(R1);

}

  
