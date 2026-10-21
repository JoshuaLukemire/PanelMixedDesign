// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "approxutils.h"
#include <limits>

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




// Penalized negative log-likelihood for the current response and random effects.
// Also update prob so that the objective and derivatives use the same point.
double calculate_pql_objective(Eigen::VectorXd Y,
                               Eigen::MatrixXd X,
                               Eigen::VectorXd beta,
                               Eigen::VectorXd u,
                               Eigen::MatrixXd OmegaRE,
                               int nChoiceSet,
                               Eigen::VectorXd & prob){

  int choice_set_size = X.rows() / nChoiceSet;
  Eigen::VectorXd XB = X * beta;
  Eigen::VectorXd XB_cs = Eigen::VectorXd::Zero(choice_set_size);
  Eigen::VectorXd expXB_cs = Eigen::VectorXd::Zero(choice_set_size);
  double objective = 0.5 * u.dot(OmegaRE * u);
  double max_XB = 0.0;
  double sum_expXB = 0.0;

  for (int ichoice = 0; ichoice < nChoiceSet; ichoice++){
    // Center utilities within each choice set before exponentiating.
    XB_cs = XB.segment(ichoice*choice_set_size, choice_set_size);
    max_XB = XB_cs.maxCoeff();
    XB_cs = XB_cs.array() - max_XB;
    expXB_cs = XB_cs.array().exp();
    sum_expXB = expXB_cs.sum();
    prob.segment(ichoice*choice_set_size, choice_set_size) = expXB_cs / sum_expXB;
    objective += log(sum_expXB) -
      Y.segment(ichoice*choice_set_size, choice_set_size).dot(XB_cs);
  }

  return objective;
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
  
  if (n_samples <= 0 || nChoiceSet <= 0 || X.rows() == 0 ||
      X.rows() % nChoiceSet != 0 || X.rows() / nChoiceSet < 2 ||
      X.cols() != b_mean.size() || var_vec.size() != b_mean.size() ||
      !X.allFinite() || !b_mean.allFinite() || !var_vec.allFinite() ||
      (var_vec.array() < 0).any() || !(var_vec.array() > 0).any()){
    Rcpp::stop("PQL requires consistent finite inputs, positive sample/choice-set counts, and at least one random effect.");
  }

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
  if (!(var_vec.head(n_random_effect).array() > 0).all()){
    Rcpp::stop("PQL requires random-effect columns first; use varcov_approx_PML().");
  }

  
  // Allocate space for final precision about the model parameters
  Eigen::MatrixXd vcov_inverse = Eigen::MatrixXd::Zero(n_beta + n_random_effect,
                                                       n_beta + n_random_effect);
  
  // Precision matrix only for random effects
  Eigen::MatrixXd OmegaRE = var_vec.head(n_random_effect).array().inverse().matrix().asDiagonal();
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
  int max_iter                 = 200;
  bool converged               = false;
  double stepsize              = 0.0;
  double tol                   = 0.00000001; // tolerance on the gradient
  int max_backtrack            = 50;
  bool accepted               = false;
  double objective            = 0.0;
  double objective_new        = 0.0;
  double decrease             = 0.0;
  double objective_tolerance  = 0.0;
  Eigen::VectorXd beta_new    = b_mean;
  Eigen::VectorXd prob_new    = Eigen::VectorXd::Zero(n_alternative_total);
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
    // As a starting value, pick the most likely Y for each choice set
    for (int ichoice = 0; ichoice < nChoiceSet; ichoice++){
      yind=0;
      // Find the selected value for this choice set
      for (int ialt = 1; ialt < choice_set_size; ialt++){
        if (XB(ichoice*choice_set_size + ialt) > XB(ichoice*choice_set_size + yind)) {
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
    
    // Reset convergence indicator. Retry from zero if needed.
    converged = false;
    for (int iattempt = 0; iattempt < 2 && !converged; iattempt++){
      if (iattempt == 1){
        u0_mle.setZero();
      }

      for(int inewton_iter = 0; inewton_iter < max_iter; inewton_iter++){
        beta.block(0, 0, n_random_effect, 1) = b_mean.segment(0, n_random_effect) + u0_mle;
        beta.block(n_random_effect, 0, n_fixed_effect, 1) = b_mean.segment(n_random_effect, n_fixed_effect);

        objective = calculate_pql_objective(Y, X, beta, u0_mle, OmegaRE, nChoiceSet, prob);
        gradient = Xr.transpose() * (Y - prob) - OmegaRE * u0_mle;
        if (!std::isfinite(objective) || !gradient.allFinite()){
          break;
        }

        // A small damped step alone is not enough to establish convergence.
        if (gradient.cwiseAbs().maxCoeff() < tol){
          converged = true;
          break;
        }

        delta = calculate_delta(prob, nChoiceSet);
        H = -(Xr.transpose() * delta * Xr) - OmegaRE;
        Hinv = H.inverse();
        decrease = -gradient.dot(Hinv * gradient);
        if (!std::isfinite(decrease) || decrease <= 0.0){
          break;
        }

        // Try the full Newton step, then halve it until the objective decreases.
        stepsize = 1.0;
        accepted = false;
        objective_tolerance = 16.0 * std::numeric_limits<double>::epsilon() *
          (1.0 + fabs(objective));
        for (int ibacktrack = 0; ibacktrack < max_backtrack; ibacktrack++){
          u1_mle = u0_mle - stepsize * (Hinv * gradient);
          beta_new = b_mean;
          beta_new.head(n_random_effect) += u1_mle;
          objective_new = calculate_pql_objective(Y, X, beta_new, u1_mle, OmegaRE, nChoiceSet, prob_new);

          if (std::isfinite(objective_new) &&
              objective_new <= objective - 0.0001 * stepsize * decrease + objective_tolerance){
            u0_mle = u1_mle;
            accepted = true;
            break;
          }
          stepsize *= 0.5;
        }
        if (!accepted){
          break;
        }
      }

      // Check the final accepted point and refresh probabilities for information.
      beta = b_mean;
      beta.block(0, 0, n_random_effect, 1) += u0_mle;
      objective = calculate_pql_objective(Y, X, beta, u0_mle, OmegaRE, nChoiceSet, prob);
      gradient = Xr.transpose() * (Y - prob) - OmegaRE * u0_mle;
      converged = std::isfinite(objective) && gradient.allFinite() &&
        gradient.cwiseAbs().maxCoeff() < tol;
    }

    if (!converged){
      Rcpp::stop("PQL mode solver failed for draw %d of %d; covariance was not computed.",
                 iy + 1, n_samples);
    }
    delta = calculate_delta(prob, nChoiceSet);

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
      I22 += 2.0 * (SigmaR.array()).sqrt().matrix().asDiagonal() * 
        I11Temp.array().pow(2).matrix() * 
        SigmaR.array().sqrt().matrix().asDiagonal();
      
      
    }
    
  } // end of loop over Y
  
  vcov_inverse.block(0, 0, n_beta, n_beta) = I11;
  vcov_inverse.block(n_beta, n_beta, n_random_effect, n_random_effect) = I22.block(0, 0, n_random_effect, n_random_effect);
  // Normalize by number of replicates of Y
  for (int i = 0; i < n_beta + n_random_effect; i++){
    for (int j = 0; j < n_beta + n_random_effect; j++){
      vcov_inverse(i,j) *= 1.0 / n_samples;
    }
  }
  
  // Return variance-covariance matrix
  return(vcov_inverse.inverse());
  
}
