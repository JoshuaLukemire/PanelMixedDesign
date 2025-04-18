// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "approxutils.h"

// [[Rcpp::depends(RcppEigen)]]

// Helper function to get the matrix of difference codings
Eigen::MatrixXd get_difference_codings(Eigen::MatrixXd X,
                                     int nChoiceSet,
                                     int choice_set_size){
  
  // Storage for difference codings
  Eigen::MatrixXd XDC = Eigen::MatrixXd::Zero(nChoiceSet * (choice_set_size - 1), X.cols());
  
  // Calculate difference relative to last row
  for (int i = 0; i < nChoiceSet; i++){
    for (int j = 0; j < choice_set_size - 1; j++){
      XDC.row(i*(choice_set_size-1) + j) = X.row(i*choice_set_size + j) - X.row( (i+1)*choice_set_size - 1);
    }
  }
  
  return XDC;
}



//' Approximation to the information matrix using PQL
//'
//' @description Approximates the information matrix for the model parameters 
//' using partial quasilikelihood. This is not intended to be directly used,
//' instead a user should use the \code{\link{varcov_approx_PML}} function with 
//' the "PQL" option.
//'
//' @param n_samples Number of samples to draw
//' @param X The model matrix
//' @param b_mean The vector of effect means
//' @param var_vec The vector of effect variances
//' @param nChoiceSet The number of choice sets in the discrete choice experiment
//' @return A matrix containing the PQL approximation to the variance
//' covariance matrix of the model parameters
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd PQL_approx_vcov(int n_samples,
                          Eigen::MatrixXd X,
                          Eigen::VectorXd b_mean,
                          Eigen::VectorXd var_vec,
                          int nChoiceSet){
  
  // Determine some aspects of the design problem from the inputs
  int n_alternative_total  = X.rows(); // total number of alternatives
  int n_beta               = b_mean.size();
  int choice_set_size      = n_alternative_total / nChoiceSet;
  
  // Count how many effects are random
  int n_random_effect = 0;
  for (int i = 0; i < var_vec.size(); i++){
    if (var_vec(i) > 0){
      n_random_effect++;
    }
  }
  int n_fixed_effect = var_vec.size() - n_random_effect;
  
  // Allocate space for final precision about the model parameters
  Eigen::MatrixXd vcov_inverse = Eigen::MatrixXd::Zero(n_beta + n_random_effect,
                                                       n_beta + n_random_effect);
  
  // Inverse of prior covariance
  Eigen::MatrixXd Omega = var_vec.array().inverse().matrix().asDiagonal();
  
  // Precision matrix only for random effects 
  Eigen::MatrixXd OmegaRE = Omega.block(0, 0, n_random_effect, n_random_effect);
  Eigen::VectorXd SigmaR  = var_vec;

  // Vector to store the standard deviations for each term
  Eigen::VectorXd sd_vec = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::MatrixXd u      = Eigen::MatrixXd::Zero(n_random_effect, n_samples);
  for (int ip = 0; ip < n_random_effect; ip++){
    sd_vec(ip) = sqrt(var_vec(ip));
    for (int iu = 0; iu < n_samples; iu++){
      u(ip, iu) = R::rnorm(0.0, 1.0);
    }
  }
  Eigen::MatrixXd sd_mat = sd_vec.asDiagonal();
  
  // Get the betas for each sample
  Eigen::MatrixXd beta_start = Eigen::MatrixXd::Zero(b_mean.size(), n_samples); // each column is a sample of beta
  for (int iy = 0; iy < n_samples; iy++){
    beta_start.block(0, iy, n_random_effect, 1)            = b_mean.segment(0, n_random_effect) + sd_mat *  u.col(iy);
    beta_start.block(n_random_effect, iy, n_fixed_effect, 1) = b_mean.segment(n_random_effect, n_fixed_effect);
  }
  
  // Calculate the epsilon terms (get added to linear predictor later)
  Eigen::MatrixXd epsilon = Eigen::MatrixXd::Zero(X.rows(), n_samples);
  for (int iy = 0; iy < n_samples; iy++){
    epsilon.col(iy) = sample_epsilon(n_alternative_total);
  }
  
  //
  // Calculate the linear predictor for each Y, draw the most likely Y
  Eigen::VectorXd XB    = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd expXB = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd prob  = Eigen::VectorXd::Zero(n_alternative_total);
  Eigen::MatrixXd Yall  = Eigen::MatrixXd::Zero(X.rows(), n_samples);
  int yind              = 0; // for tracking which y is selected
  
  // Allocate space for intermediate matrix terms
  Eigen::MatrixXd I11           = Eigen::MatrixXd::Zero(n_beta, n_beta);
  Eigen::MatrixXd I11Temp       = Eigen::MatrixXd::Zero(n_beta, n_beta);
  Eigen::MatrixXd I12           = Eigen::MatrixXd::Zero(n_beta, n_random_effect);
  Eigen::MatrixXd I22           = Eigen::MatrixXd::Zero(n_beta, n_beta);
  Eigen::MatrixXd beta          = Eigen::MatrixXd::Zero(n_beta, 1);
  
  // Allocate space for terms in gradient descent
  int max_iter                 = 100;
  bool converged               = false;
  double stepsize              = 0.0;
  double tol                   = 0.000000001;
  int n_samples_converged      = n_samples; // number of samples that converge for finding tilde(u)
  double iterchange            = 0.0; // track convergence 
  Eigen::VectorXd Y            = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd u0_mle       = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd u1_mle       = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::MatrixXd delta        = Eigen::MatrixXd::Zero(n_alternative_total, n_alternative_total);
  Eigen::MatrixXd deltaDX      = Eigen::MatrixXd::Zero(nChoiceSet * (choice_set_size - 1),
                                                       nChoiceSet * (choice_set_size - 1));
  Eigen::MatrixXd deltaDXInv   = Eigen::MatrixXd::Zero(nChoiceSet * (choice_set_size - 1),
                                                       nChoiceSet * (choice_set_size - 1));
  Eigen::MatrixXd Vn           = Eigen::MatrixXd::Zero(nChoiceSet * (choice_set_size - 1),
                                                       nChoiceSet * (choice_set_size - 1));
  Eigen::MatrixXd VnInv        = Eigen::MatrixXd::Zero(nChoiceSet * (choice_set_size - 1),
                                                       nChoiceSet * (choice_set_size - 1));
  Eigen::VectorXd gradient     = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::MatrixXd H            = Eigen::MatrixXd::Zero(n_random_effect, n_random_effect);
  Eigen::MatrixXd Hinv         = Eigen::MatrixXd::Zero(n_random_effect, n_random_effect);

  // Loop over number of samples
  for (int iy = 0; iy < n_samples; iy++){
    XB    = X * beta_start.col(iy) + epsilon.col(iy);
    expXB = XB.array().exp();
    // As a starting value, pick the most likely Y for each choice set
    for (int ichoice = 0; ichoice < nChoiceSet; ichoice++){
      yind=0;
      // Find the selected value for this choice set
      for (int ialt = 1; ialt < choice_set_size; ialt++){
        if (expXB(ichoice*choice_set_size + ialt) > expXB(ichoice*choice_set_size + yind)) {
          yind = ialt;
        }
      }
      Yall(ichoice*choice_set_size + yind, iy) = 1;
    }
  }
  
  // Get the model matrix for differences
  Eigen::MatrixXd XDX =  get_difference_codings(X, nChoiceSet, choice_set_size);

  // design matrix only corresponding to the random effects.
  Eigen::MatrixXd Xr = X.block(0, 0, X.rows(), n_random_effect);
  
  // Loop over samples of Y/U
  for (int iy = 0; iy < n_samples; iy++){
    
    // Draws for Y and u
    Y         = Yall.col(iy);
    u0_mle    = u.col(iy);
    
    // Reset convergence indicator
    converged = false;
    
    // Gradient descent loop
    for(int inewton_iter = 0; inewton_iter < max_iter; inewton_iter++){
      
      // Update beta using current value for U0
      beta.block(0, 0, n_random_effect, 1)              = b_mean.segment(0, n_random_effect) + u0_mle;
      beta.block(n_random_effect, 0, n_fixed_effect, 1) = b_mean.segment(n_random_effect, n_fixed_effect);
      
      
      // Calculate response probabilities under U0_MLE
      prob = calculate_response_probabilities(X, beta, nChoiceSet);
  
      // calculate diag(p) - pp'
      delta = calculate_delta(prob, nChoiceSet);
      
      // from all rand effects version
      gradient = Xr.transpose() * (Y - prob) -  OmegaRE * u0_mle;
      
      H = -(Xr.transpose() * delta * Xr) - OmegaRE; 
      
      //stepsize = 1.0 / (1.0 + 10.0 * sqrt(gradient.array().pow(2).sum() ) );
      stepsize = 1.0;
      
      Hinv = H.inverse();
      
      u1_mle = u0_mle - stepsize*(Hinv * gradient);
      
      iterchange = 0.0;
      for (int ib = 0; ib < u1_mle.size(); ib++){
        iterchange += sqrt(pow(u1_mle(ib) - u0_mle(ib), 2));
      }
      
      if ( iterchange < tol ){
        converged = true;
        inewton_iter = max_iter + 1;
      } else {
        u0_mle = u1_mle;
      }
    } // end of newton iterations
    
    
    if (converged == true){
      
      // Grab the elements of the Delta matrix corresponding to difference matrix (XDX)
      for (int ics = 0; ics < nChoiceSet; ics++){
        deltaDX.block(ics*(choice_set_size-1), ics*(choice_set_size-1),
                      choice_set_size-1, choice_set_size-1) = delta.block(ics*choice_set_size, ics*choice_set_size,
                      choice_set_size-1, choice_set_size-1);
      }
      
      // stacked V term for precision calculations
      deltaDXInv = deltaDX.inverse();
      Vn = deltaDXInv + XDX * var_vec.matrix().asDiagonal() * XDX.transpose();
      VnInv = Vn.inverse();
      
      // Increment components of precision matrix
      I11Temp = XDX.transpose() * VnInv * XDX;
      I11 += I11Temp;
      I22 += (0.5 * SigmaR.array()).sqrt().matrix().asDiagonal() * 
        I11Temp.array().pow(2).matrix() * 
        SigmaR.array().sqrt().matrix().asDiagonal();
      
      
    } else { // case where outer newton did not converge, decrease n_samples_converged counter
      n_samples_converged -= 1;
    }
    
  } // end of loop over Y
  
  vcov_inverse.block(0, 0, n_beta, n_beta) = I11;
  vcov_inverse.block(n_beta, n_beta, n_random_effect, n_random_effect) = I22.block(0, 0, n_random_effect, n_random_effect);
  // Normalize by number of replicates of Y
  for (int i = 0; i < n_beta + n_random_effect; i++){
    for (int j = 0; j < n_beta + n_random_effect; j++){
      vcov_inverse(i,j) *= 1.0 / n_samples_converged;
    }
  }
  
  // Return variance-covariance matrix
  return(vcov_inverse.inverse());
  
}
