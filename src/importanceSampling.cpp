// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "approxutils.h"

// [[Rcpp::depends(RcppEigen)]]

// Function to draw a sample of responses based on the linear predictor
Eigen::MatrixXd impsamp_draw_response_vector(Eigen::VectorXd expXB,
                                             Eigen::VectorXd X,
                                             int n_choice_set){
  
  // Get the number of alternatives per choice set
  int choice_set_size = X.rows() / n_choice_set;
  
  // Allocate space for the responses
  Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(X.rows(), 1);
  
  // Initialize variables for loop
  int ind           = 0;
  double maxval     = 0.0;
  double total_prob = 0.0;
  double cumu_prob  = 0.0;
  double rdraw      = 0.0;
  
  // Determine the most likely outcome for this sampled Y
  for (int ic = 0; ic < n_choice_set; ic++){
    
    total_prob = 0.0;
    for (int ialt = 0; ialt < choice_set_size; ialt++){
      total_prob += expXB( ic*choice_set_size + ialt );
    }
    
    // Reset
    ind = 0;
    double cumu_prob = 0.0;
    
    // Not exported via Rcpp, so need to pull RNG info
    GetRNGstate();
    rdraw = R::runif(0.0, 1.0);
    PutRNGstate();
    
    // Loop over alternatives and pick the most likely
    for (int ialt = 0; ialt < choice_set_size; ialt++){
      cumu_prob += expXB( ic*choice_set_size + ialt ) / total_prob;
      if (cumu_prob > rdraw){
        ind = ic * choice_set_size + ialt;
        rdraw = 999.0;
      }
    }
    // Set Y for this choice set to 1
    Y(ind, 0) = 1;
  }
  
  return Y;
}


// Function to evaluate marginal probability given a set of probs and observed Y
Eigen::VectorXd impsamp_calc_marginal_prob(Eigen::VectorXd Y,
                                           Eigen::MatrixXd probs,
                                           int n_choice_set){
  
  // Find the number of samples
  int n_samples_u = probs.cols();
  
  // Initialize a vector for the marginal probabilities
  Eigen::VectorXd marginal_probs = Eigen::VectorXd::Zero(n_samples_u);
  
  // see https://stackoverflow.com/questions/44268968/eigen-vector-or-matrix-componentwise-to-power
  Eigen::MatrixXd powerTerms = probs.array() // go to array domain
                                    .pow(                 // element-wise power
  Y.replicate(1, probs.cols()).array() // repeat exponents to match size of A
                                    );
  
  // Calc probs for each of the random effect samples
  for (int iu = 0; iu < n_samples_u; iu++){
    marginal_probs(iu) = powerTerms.col(iu).prod();
  }
  
  return marginal_probs;
}


 // [[Rcpp::export]]
 Eigen::MatrixXd impsamp_calc_info_given_Y(Eigen::VectorXd Y,
                                           Eigen::MatrixXd X,
                                           Eigen::VectorXd b_mean,
                                           Eigen::VectorXd var_vec,
                                           int n_choice_set,
                                           int n_samples_u,
                                           bool return_unweighted,
                                           Eigen::VectorXd & marginal_prob_store){
   
   // Track how man_samples_Y effects are random
   int n_random_effect = 0;
   for (int i = 0; i < var_vec.size(); i++){
     if (var_vec(i) > 0){
       n_random_effect++;
     }
   }
   // Leftover are fixed effects
   int n_fixed_effect = var_vec.size() - n_random_effect;
   
   // Number of beta terms
   int effect_dim = b_mean.size();
   
   // Draw random effect samples and calculate corresponding effects (beta)
   // each row is an effect and each column is a sample
   Eigen::MatrixXd u = Eigen::MatrixXd::Zero(n_random_effect, n_samples_u);
   Eigen::MatrixXd beta_i = Eigen::MatrixXd::Zero(effect_dim, n_samples_u);
   Eigen::VectorXd sd_vec = Eigen::VectorXd::Zero(n_random_effect);
   GetRNGstate();
   for (int ip = 0; ip < n_random_effect; ip++){
     sd_vec(ip) = sqrt(var_vec(ip));
     for (int iu = 0; iu < n_samples_u; iu++){
       u(ip, iu) = R::rnorm(0.0, 1.0);
     }
   }
   PutRNGstate();
   Eigen::MatrixXd sd_mat = sd_vec.asDiagonal();
   
   // Use the random draws to get the samples of the coefficients
   for (int iu = 0; iu < n_samples_u; iu++){
     beta_i.block(0,iu,  n_random_effect, 1) = 
       b_mean.segment(0, n_random_effect) + sd_mat *  u.col(iu);
     beta_i.block(n_random_effect, iu, n_fixed_effect, 1) = 
       b_mean.segment(n_random_effect, n_fixed_effect);
   }
   
   // Calculate the response probabilities
   Eigen::MatrixXd response_probabilities = Eigen::MatrixXd::Zero(Y.size(), n_samples_u);
   for (int iu = 0; iu < n_samples_u; iu++){
     response_probabilities.col(iu) = calculate_response_probabilities(X, beta_i.col(iu), n_choice_set);
   }
   
   // Calculate the marginal probabilities, each column corresponds to a u_i
   Eigen::VectorXd marginal_probs = impsamp_calc_marginal_prob(Y, response_probabilities, n_choice_set);
   
   // Average of the probabilities from above over n_samples_u is the marginal prob
   double marg_prob = marginal_probs.sum() / n_samples_u; // P(Y = y)
   
   Eigen::VectorXd f1 = (response_probabilities.matrix() * marginal_probs.matrix()).array() * (1.0 / n_samples_u) * (1.0 / marg_prob);
   
   // score function for the regression coefficients
   Eigen::VectorXd score11 = X.matrix().transpose() * (Y.array() - f1.array()).matrix();
   
   // Upper P x P block of the information matrix
   Eigen::MatrixXd I11 = (score11 * score11.transpose()).array() * marg_prob;  
   
   // U22 (for the lower 2 x 2 block)
   Eigen::MatrixXd U22 = Eigen::MatrixXd::Zero(n_random_effect, n_samples_u);
   for (int iU = 0; iU < n_samples_u; iU++){
     for (int iP = 0; iP < n_random_effect; iP++){
       U22(iP, iU) = pow(sd_vec(iP) * u(iP, iU), 2) / pow(sqrt(var_vec(iP)), 3);
     }
   }
   
   Eigen::VectorXd f2 = (U22.matrix() * marginal_probs.matrix()).array() * (1.0 / n_samples_u) * (1.0 / marg_prob);
   //Eigen::VectorXd f2 = (U22.matrix() * marginal_probs.matrix()).array() * (1.0 / n_samples_u);
   Eigen::VectorXd score2 = sd_vec.array().inverse() * (-1.0) + f2.array();
   
   Eigen::MatrixXd I12 = (score11.matrix() * score2.matrix().transpose()).array() * marg_prob;
   Eigen::MatrixXd I22 = (score2.matrix()  * score2.matrix().transpose()).array() * marg_prob;
   
   Eigen::MatrixXd information_matrix = Eigen::MatrixXd::Zero( effect_dim + n_random_effect,
                                                       effect_dim + n_random_effect );
   
   // Plug in components of information matrix calculated above
   information_matrix.block(0, 0, effect_dim, effect_dim) = I11;
   information_matrix.block(0, effect_dim, effect_dim, n_random_effect) = I12;
   information_matrix.block(effect_dim, 0, n_random_effect, effect_dim) = I12.transpose();
   information_matrix.block(effect_dim, effect_dim, n_random_effect, n_random_effect) = I22.block(0, 0, n_random_effect, n_random_effect);
   
   // Store the marginal probability
   marginal_prob_store(0) = marg_prob;
   if (return_unweighted == true){
     information_matrix = information_matrix / marg_prob;
   }
   
   // Return information matrix to main function
   return information_matrix;
 }


//' Approximate the information matrix using importance sampling based on a complete enumeration of possible responses
 //'
 //' @description Approximates the variance-covariance matrix for the model parameters 
 //' using importance sampling. This is not intended to be directly used,
 //' instead a user should use the \code{\link{varcov_approx_PML}} function with 
 //' the "importance" option.
 //'
 //' @param Y an enumeration of all possible responses. Can be obtained using
 //'  the \code{\link{gen_all_choice_seq}} function
 //' @param X The model matrix
 //' @param b_mean The vector of effect means
 //' @param var_vec The vector of effect variances
 //' @param n_choice_set The number of choice sets in the discrete choice experiment
 //' @param n_samples_u Number of samples of random effects to use
 //' @return A matrix containing the variance covariance matrix of the model parameters
 //' @keywords internal
 // [[Rcpp::export]]
 Eigen::MatrixXd importance_sample_all_Y(Eigen::MatrixXd Y,
                                        Eigen::MatrixXd X,
                                        Eigen::VectorXd b_mean,
                                        Eigen::VectorXd var_vec, 
                                        int n_choice_set, 
                                        int n_samples_u){
   
   // Track how man_samples_Y effects are random
   int n_random_effect = 0;
   for (int i = 0; i < var_vec.size(); i++){
     if (var_vec(i) > 0){
       n_random_effect++;
     }
   }
   
   // Allocate space for information matrix and marginal probs
   Eigen::MatrixXd FI = Eigen::MatrixXd::Zero(b_mean.size() + n_random_effect,
                                              b_mean.size() + n_random_effect);
   Eigen::VectorXd marginal_prob_store = Eigen::VectorXd::Zero(1);
   
   // Get the number of unique combinations of Y
   int n_samples_Y = Y.cols();
   
   // This is for debugging
   double total_probability = 0.0;
   
   // Loop over configurations of Y and eval contribution to information
   // matrix for each
   for (int iY = 0; iY < n_samples_Y; iY++){
     FI += impsamp_calc_info_given_Y(Y.col(iY), X, b_mean, var_vec, 
                                     n_choice_set, n_samples_u, false, 
                                     marginal_prob_store);
     total_probability += marginal_prob_store(0);
   }
   
   // Return variance covariance
   return FI.inverse();
   
 }




//' Approximate the information matrix using importance sampling based on random samples of possible responses
 //'
 //' @description Approximates the variance-covariance matrix for the model parameters 
 //' using importance sampling. This is not intended to be directly used,
 //' instead a user should use the \code{\link{varcov_approx_PML}} function with 
 //' the "importance" option.
 //'
 //' @param X The model matrix
 //' @param b_mean The vector of effect means
 //' @param var_vec The vector of effect variances
 //' @param n_choice_set The number of choice sets in the discrete choice experiment
 //' @param n_samples_u Number of samples of random effects to use
 //' @param n_samples_Y Number of samples of possible response vectors to use
 //' @return A matrix containing the variance covariance matrix of the model parameters
 //' @keywords internal
 // [[Rcpp::export]]
Eigen::MatrixXd importance_sample_random_Y(     Eigen::MatrixXd X,
                                      Eigen::VectorXd b_mean,
                                      Eigen::VectorXd var_vec, 
                                      int n_choice_set, 
                                      int n_samples_u, 
                                      int n_samples_Y){
  
  int n_random_effect = 0;
  for (int i = 0; i < var_vec.size(); i++){
    if (var_vec(i) > 0){
      n_random_effect++;
    }
  }
  
  Eigen::MatrixXd FI = Eigen::MatrixXd::Zero(b_mean.size() + n_random_effect, b_mean.size() + n_random_effect);
  Eigen::MatrixXd Y  = Eigen::MatrixXd::Zero(X.rows(), 1);
  Eigen::MatrixXd epsilon = Eigen::MatrixXd::Zero(X.rows(), 1);
  Eigen::VectorXd expXB  = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd beta_n = Eigen::VectorXd::Zero(b_mean.size());
  
  Eigen::VectorXd XB    = Eigen::VectorXd::Ones(X.rows());
  Eigen::VectorXd Py    = Eigen::VectorXd::Zero(X.rows());
  
  Eigen::VectorXd marginal_prob_store = Eigen::VectorXd::Zero(1);
  
  double total_probability = 0.0;
  double total_probability_check = 0.0;
  
  for (int iY = 0; iY < n_samples_Y; iY++){
    
    // Get a draw of Y to use
    for (int p = 0; p < b_mean.size(); p++){
      if (p < n_random_effect){
        beta_n(p) = b_mean(p) + sqrt(var_vec(p)) * R::rnorm(0.0, 1.0);
      } else {
        beta_n(p) = b_mean(p);
      }
    }
    XB    = X * beta_n;
    expXB = XB.array().exp();
    Y = impsamp_draw_response_vector(expXB, X, n_choice_set);
    
    // Use importance sampling at this value of Y
    FI += impsamp_calc_info_given_Y(Y.col(0), X, b_mean, var_vec, n_choice_set,
                                    n_samples_u, true, marginal_prob_store);
    
  }
  FI = FI / n_samples_Y;
  
  // Invert to get variance-covariance matrix
  FI = FI.inverse();
  
  return FI;
  
}