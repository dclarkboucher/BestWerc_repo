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

double divg(double x, double shape, double scale){
  
  return(R::dgamma(1/x, shape, 1/scale, 0) * pow(x, -2.0));
  
}

double pivg(double x, double shape, double scale){
  
  return(R::pgamma(1/x, shape, 1/scale, 0, 0));
}

double qivg(double p, double shape, double scale){
  
  return(1 / R::qgamma(1-p, shape, 1/scale, 1, 0));
  
}


// [[Rcpp::export]] 
double gamma_diff_prob(double Delta, double alpha0, double beta0, double alpha1, double beta1,
                       double gridsize, double tail_error){
  
  double gamma0_lower = qivg(tail_error/2, alpha0, beta0);
  double gamma0_upper = qivg(1-tail_error/2, alpha0, beta0);
  double dx0 = (gamma0_upper - gamma0_lower)/gridsize;
  double prob = 0;
  for (double g0 = gamma0_lower; g0 <= gamma0_upper; g0 = g0 + dx0){
    prob = prob + dx0 * divg(g0, alpha0, beta0) * (1 - pivg(g0 + Delta,alpha1,beta1)); 
  }
  
  // double dx0 = 1.0 / (gridsize + 1.0);
  // arma::vec percentiles(gridsize + 1);
  // percentiles(0) = 0;
  // for (int i = 1; i < gridsize +1; i++){
  //   percentiles(i) = qivg(dx0 * i);
  // }
  // 
  // 
  // 
  // for(double g0_p = dx0;  g0_p < 1.0; g0_p = g0_p + dx0){
  //   double g0 = qivg(g0_p, alpha0, beta0);
  //   double g0_last = qivg(g0_p)
  //   
  // }
  
  
  return(prob);
  
}

// [[Rcpp::export]]
double survival_power_random(double alpha, double beta, double n_0, double n_1, 
                             double gamma_0, double gamma_1, double gamma_c, 
                             double Delta = 0, double alpha_0 = 0.05, int reps = 1000,
                             double gridsize = 100, bool fix_cens = false
){
  
  double power = 0;
  double a_0;
  double b_0;
  double a_1;
  double b_1;
  double prob;
  double gamma_x0 = 1.0/(1.0/gamma_0 + 1.0/gamma_c); 
  double gamma_x1 = 1.0/(1.0/gamma_1 + 1.0/gamma_c); 
  double prob0 = gamma_x0 / gamma_0; // Event prob
  double prob1 = gamma_x1 / gamma_1;
  double x_0;
  double x_1;
  double u_0 = prob0 * n_0;
  double u_1 = prob1 * n_1;
  for (int r = 0; r < reps; r++){
    set_seed(r);
    x_0 = R::rgamma(n_0, gamma_x0);
    x_1 = R::rgamma(n_1, gamma_x1);
    if (!fix_cens){
      u_0 = R::rbinom(n_0, prob0);
      u_1 = R::rbinom(n_1, prob1);
    }
    a_0 = alpha + u_0;
    b_0 = beta + x_0;
    a_1 = u_1;
    b_1 = x_1;
    prob = gamma_diff_prob(
      Delta, a_0, b_0, a_1, b_1, gridsize, 1e-6
    );
    if (prob > 1 - alpha_0){
      power = power + 1;
      
    }
  }
  return(power / double(reps));
}

// [[Rcpp::export]]
double survival_power_numerical(double alpha, double beta, double n_0, double n_1, 
                                double gamma_0, double gamma_1, double gamma_c, 
                                double Delta = 0, double alpha_0 = 0.05, 
                                double gridsize = 100, bool fix_cens = false
){
  double power = 0;
  double a_0;
  double b_0;
  double a_1;
  double b_1;
  double prob;
  
  double tail_error = 1 - pow(1.0 - 0.001, 0.25);
  if (fix_cens) tail_error = 1 - pow(1.0 - 0.001, 0.5);
  
  double gamma_x0 = 1.0/(1.0/gamma_0 + 1.0/gamma_c); 
  double gamma_x1 = 1.0/(1.0/gamma_1 + 1.0/gamma_c); 
  double prob0 = gamma_x0 / gamma_0; // Event prob
  double prob1 = gamma_x1 / gamma_1;
  
  // Integration limits
  double x_0_lower = R::qgamma(tail_error/2, n_0, gamma_x0, 1, 0);
  double x_0_upper = R::qgamma(tail_error/2, n_0, gamma_x0, 0, 0);
  double x_1_lower = R::qgamma(tail_error/2, n_1, gamma_x1, 1, 0);
  double x_1_upper = R::qgamma(tail_error/2, n_1, gamma_x1, 0, 0);
  double dx0 = (x_0_upper - x_0_lower) / gridsize;
  double dx1 = (x_1_upper - x_1_lower) / gridsize;
  
  if (fix_cens){
    
    double u_0 = prob0 * n_0;
    double u_1 = prob1 * n_1;
    
    // Integration limits
    for (double x_0 = x_0_lower; x_0 < x_0_upper; x_0 = x_0 + dx0){
      for (double x_1 = x_1_lower; x_1 < x_1_upper; x_1 = x_1 + dx1){
        // Calculate the decision rule
        a_0 = alpha + u_0;
        b_0 = beta + x_0;
        a_1 = u_1;
        b_1 = x_1;
        prob = gamma_diff_prob(
          Delta, a_0, b_0, a_1, b_1, gridsize, 1e-6
        );
        if (prob > 1 - alpha_0){
          power = power + dx0 * dx1 * R::dgamma(x_0, n_0, gamma_x0, 0) *
            R::dgamma(x_1, n_1, gamma_x1, 0);
        }
      }
    }
    
  } else {
    double u_0_lower = R::qbinom(tail_error/2, n_0, prob0, 1, 0);
    double u_1_lower = R::qbinom(tail_error/2, n_1, prob1, 1, 0);
    double u_0_upper = R::qbinom(tail_error/2, n_0, prob0, 0, 0);
    double u_1_upper = R::qbinom(tail_error/2, n_1, prob1, 0, 0);
    for (int u_0 = u_0_lower; u_0 <= u_0_upper; u_0 ++ ){
      for (int u_1 = u_1_lower; u_1 <= u_1_upper; u_1 ++ ){
        for (double x_0 = x_0_lower; x_0 < x_0_upper; x_0 = x_0 + dx0){
          for (double x_1 = x_1_lower; x_1 < x_1_upper; x_1 = x_1 + dx1){
            // Calculate the decision rule
            a_0 = alpha + u_0;
            b_0 = beta + x_0;
            a_1 = u_1;
            b_1 = x_1;
            prob = gamma_diff_prob(
              Delta, a_0, b_0, a_1, b_1, gridsize, 1e-6
            );
            if (prob > 1 - alpha_0){
              power = power + dx0 * dx1 * R::dgamma(x_0, n_0, gamma_x0, 0) *
                R::dgamma(x_1, n_1, gamma_x1, 0) * R::dbinom(u_0, n_0, prob0, 0) *
                R::dbinom(u_1, n_1, prob1, 0);
            }
          }
        }
      }
    }
    
  }
  
  return(power);
  
}

// [[Rcpp::export]]
List survival_power_vectorized(
    arma::vec gamma_0_list,
    arma::vec delta_list,
    arma::vec gamma_c_list,
    arma::vec n_1_list,
    arma::vec n_0_list,
    arma::vec a_0_list,
    arma::vec r_list,
    double alpha_0,
    int reps,
    double gridsize = 100,
    double Delta = 0,
    bool numerical = false,
    bool fix_cens = false
){
  
  int size = gamma_0_list.n_elem;
  double gamma_0;
  double delta;
  double gamma_c;
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
    
    gamma_0 = gamma_0_list(i);
    delta = delta_list(i);
    gamma_c = gamma_c_list(i);
    n_0 = n_0_list(i);
    n_1 = n_1_list(i);
    a_0 = a_0_list(i);
    r = r_list(i);
    n_0h = a_0 / (1 - a_0) * double(n_0);
    r_0h = r * delta + gamma_0;
    alpha = n_0h;
    beta = n_0h * r_0h;
    
    if (numerical){
      
      power_list(i) = survival_power_numerical(alpha, beta, n_0, n_1, 
                 gamma_0, gamma_0+delta, gamma_c, 
                 Delta, alpha_0, gridsize, fix_cens);
      t1e_list(i) = survival_power_numerical(alpha, beta, n_0, n_1, 
               gamma_0, gamma_0, gamma_c, 
               Delta, alpha_0, gridsize, fix_cens);
      
    } else{
      power_list(i) =
        survival_power_random(
          alpha, beta, n_0, n_1, gamma_0, gamma_0 + delta, gamma_c, Delta,
          alpha_0, reps, gridsize, fix_cens
        );
      t1e_list(i) =
        survival_power_random(
          alpha, beta, n_0, n_1, gamma_0, gamma_0, gamma_c, Delta,
          alpha_0, reps, gridsize, fix_cens
        );
      
    }
    
  }
  
  List out = List::create(Named("power") = power_list, 
                          Named("type_1_error") = t1e_list);
  return(out);
}


// [[Rcpp::export]]
double find_g1(double power, double gamma_0, double gamma_c, int n_1, int n_0, 
               double Delta = 0, double alpha_0 = 0.05, int reps = 1000, double gridsize = 100){
  
  double err = 0.0001;
  double max = gamma_0 * 10.0; // This is bad but it suffices for now. 
  double min = gamma_0 / 10.0;
  double gamma_1 = 0.5 * (max + min);
  double pow;
  
  // Check maximum first.
  pow = 
    survival_power_random(
      0, 0, n_0, n_1, gamma_0, max, gamma_c, Delta, 
      alpha_0, reps, gridsize
    );
  
  if (pow < power){
    stop("Error: maximum value not high enough");
  }
  
  // Binary search
  while(max - min > err){
    
    pow = 
      survival_power_random(
        0, 0, n_0, n_1, gamma_0, gamma_1, gamma_c, Delta, 
        alpha_0, reps, gridsize
      );
    
    if (pow >= power){
      max = gamma_1;
    } else {
      min = gamma_1;
    }
    
    gamma_1 = 0.5 * (min + max);
    
  }
  
  return(gamma_1);
  
}


// [[Rcpp::export]]
double survival_est_mean(
    double alpha, double beta, int n_0, 
    double gamma_0, double gamma_c
){
  
  double mean = 0;
  // double tail_error = 0.0000001;
  double gamma_x0 = 1.0/(1.0/gamma_0 + 1.0/gamma_c); // Follow-up mean
  double prob0 = gamma_x0 / gamma_0; // Event prob
  double numerator_mean = gamma_x0 * n_0 + beta;
  
    
  // Integration limits
  // double u_0_lower = R::qbinom(tail_error/2, n_0, prob0, 1, 0);
  // double u_0_upper = R::qbinom(tail_error/2, n_0, prob0, 0, 0);
  double u_0_lower = 1; // has to be 1 to prevent errors. 
  double u_0_upper = n_0;
  
  for (int u_0 = u_0_lower; u_0 <= u_0_upper; u_0 = u_0 + 1){
    
    mean = mean + R::dbinom(u_0, n_0, prob0, 0) / (alpha + u_0); 
    
  }
  return(mean * numerator_mean);
  
} 


// [[Rcpp::export]]
double survival_est_var(
    double alpha, double beta, int n_0, 
    double gamma_0, double gamma_c
){
  
  // double tail_error = 0.000000001;
  double gamma_x0 = 1.0/(1.0/gamma_0 + 1.0/gamma_c); // Follow-up mean
  double prob0 = gamma_x0 / gamma_0; // Event prob
  
  // First compute mean
  double mean = survival_est_mean(alpha, beta, n_0, gamma_0, gamma_c);
  
  // Now get mean of the square
  double mom2_part1 = pow(beta, 2.0) + 2 * beta * n_0 * gamma_x0  + pow(n_0 * gamma_x0, 2.0) + n_0 * pow(gamma_x0, 2.0);
  double mom2_part2 = 0;
  // double u_0_lower = R::qbinom(tail_error/2, n_0, prob0, 1, 0);
  // double u_0_upper = R::qbinom(tail_error/2, n_0, prob0, 0, 0);
  double u_0_lower = 1;
  double u_0_upper = n_0;
  for (int u_0 = u_0_lower; u_0 <= u_0_upper; u_0 = u_0 + 1){
    mom2_part2 = mom2_part2 + R::dbinom(u_0, n_0, prob0, 0) / pow(alpha + u_0, 2.0); 
  }
  return(mom2_part1 * mom2_part2 - pow(mean, 2));
  
}

// [[Rcpp::export]]
List survival_bias_mse_vectorized(
    arma::vec gamma_0_list,
    arma::vec delta_list,
    arma::vec gamma_c_list,
    arma::vec n_0_list,
    arma::vec a_0_list,
    arma::vec r_list
){
  
  int size = gamma_0_list.n_elem;
  double gamma_0;
  double delta;
  double gamma_c;
  int n_0;
  double n_0h;
  double a_0;
  double r;
  double alpha;
  double beta;
  double r_0h;
  
  // Storage
  arma::vec mean_list(size);
  arma::vec var_list(size);
  
  // Loop
  for (int i = 0; i < size; i++){
    
    gamma_0 = gamma_0_list(i);
    delta = delta_list(i);
    gamma_c = gamma_c_list(i);
    n_0 = n_0_list(i);
    a_0 = a_0_list(i);
    r = r_list(i);
    n_0h = a_0 / (1 - a_0) * double(n_0);
    r_0h = r * delta + gamma_0;
    alpha = n_0h;
    beta = n_0h * r_0h;
    mean_list(i) = 
      survival_est_mean(
        alpha, beta, n_0, gamma_0, gamma_c
      );
    var_list(i) = 
      survival_est_var(
        alpha, beta, n_0, gamma_0, gamma_c
      );
  }
  List out = List::create(Named("mean") = mean_list, 
                          Named("var") = var_list);
  return(out);
}


