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

double diff_prob_fast(double alpha_n, double beta_n, double y_1, double Delta, double mu_0n, 
                 int n_1, int n_0, double n_0h){
  
  y_1 = y_1 - Delta;
  double nu = 2 * alpha_n;
  double var = beta_n / alpha_n * (1/double(n_0+ n_0h) + 1/double(n_1));
  return(R::pt((y_1 - mu_0n) / pow(var, 0.5), nu, 1,0));
}

double find_S2_fast(double alpha_0, 
               double alpha, double beta, double y_0h, double Delta, 
               double y_1, double y_0, double mu_0n, int n_1, int n_0, double n_0h){
  
  double alpha_n = alpha + (n_0 + n_1 - 1) * 0.5;
  double out = 
    alpha_n / (1/double(n_1) + 1/double(n_0 + n_0h)) * 
    pow(y_1 - mu_0n - Delta, 2.0) / pow(R::qt(1-alpha_0, 2 * alpha_n, 1, 0),2.0) - 
    n_0 * n_0h / double(n_0 + n_0h) * pow(y_0 - y_0h, 2.0) - beta;
  return(2 * out);
  
}

// [[Rcpp::export]]
double random_numerical_power(
    double delta, double mu_1, double mu_0, double r, double sigma,
    int n_0, int n_1, double n_0h, double alpha, double beta, double alpha_0, 
    double Delta, double gridsize, double tail_error = 0.0001){
  
  double individual_error = 1 - pow(1 - tail_error, 1.0/3.0);
  double y_0h = r * delta + mu_0;
  double y0_lower = R::qnorm5(individual_error/2,mu_0,sigma/pow(n_0,0.5),1, 0) ;
  double y0_upper = R::qnorm5(1 - individual_error/2, mu_0, sigma/pow(n_0,0.5), 1, 0);
  double y1_lower = R::qnorm5(individual_error/2,mu_1, sigma/pow(n_1,0.5), 1, 0);
  double y1_upper = R::qnorm5(1-individual_error/2,mu_1, sigma/pow(n_1,0.5), 1, 0);
  double y1_lower_use;
  double mu_0n;
  double power = 0;
  double S2;
  double dx0 = (y0_upper - y0_lower) / gridsize;
  double dx1 = (y1_upper - y1_lower) / gridsize;
  
  for (double y_0 = y0_lower; y_0 <= y0_upper; y_0 = y_0 + dx0){
    
    mu_0n = (n_0h * y_0h + n_0 * y_0) / (n_0h + n_0);
    y1_lower_use = max(y1_lower, mu_0n + Delta);
    
    for (double y_1 = y1_lower_use; y_1 <= y1_upper; y_1 = y_1 + dx1){
      
      // Find S2 such that below S2, result is significant
      S2 = find_S2_fast(alpha_0, alpha, beta, y_0h, Delta, 
                        y_1, y_0, mu_0n, n_1, n_0, n_0h);
      power = power + dx0 * dx1 * 
        R::dnorm4(y_0, mu_0, sigma * pow(n_0, -0.5), 0) *
        R::dnorm4(y_1, mu_1, sigma * pow(n_1, -0.5), 0) * 
        R::pchisq(S2 / pow(sigma, 2.0), n_0 + n_1 - 2, 1, 0);
      
    }
    
  }
  
  return(power);
}

// [[Rcpp::export]]
List fast_sim_numerical(
    arma::vec mu_0_list,
    arma::vec delta_list,
    arma::vec n_1_list,
    arma::vec n_0_list,
    arma::vec a_0_list,
    arma::vec r_list,
    arma::vec alpha_list,
    arma::vec beta_list,
    double alpha_0,
    double gridsize = 100,
    double tail_error = 0.0001
){
  
  int size = mu_0_list.n_elem;
  double mu_0;
  double delta;
  int n_0;
  int n_1;
  double n_0h;
  double a_0;
  double r;
  double alpha;
  double beta;
  double sigma;
  
  // Storage
  arma::vec power_list(size);
  arma::vec t1e_list(size);

  // Simulations
  for (int i = 0; i < size; i++){
    
    set_seed(633788);
    mu_0 = mu_0_list(i);
    delta = delta_list(i);
    n_0 = n_0_list(i);
    n_1 = n_1_list(i);
    a_0 = a_0_list(i);
    n_0h = a_0 / (1 - a_0) * double(n_0);
    r = r_list(i);
    alpha = alpha_list(i);
    beta = beta_list(i);
    sigma = pow(alpha/beta, -0.5);
    
    power_list(i) = random_numerical_power(
      delta, mu_0 + delta, mu_0, r, sigma,
      n_0, n_1, n_0h, alpha, beta, alpha_0, 
      0.0, gridsize, tail_error);
    t1e_list(i) = random_numerical_power(
      delta, mu_0, mu_0, r, sigma,
      n_0, n_1, n_0h, alpha, beta, alpha_0, 
      0.0, gridsize, tail_error);

  }
  
  List out = List::create(Named("power") = power_list, 
                          Named("type_1_error") = t1e_list);
  return(out);
}




