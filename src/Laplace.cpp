// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "approxutils.h"

// [[Rcpp::depends(RcppEigen)]]


//
// Newtons method to find MLE(u0)
//
Eigen::MatrixXd laplace_apply_newton_method(Eigen::MatrixXd Yall,
                                            Eigen::MatrixXd X,
                                            Eigen::MatrixXd u,
                                            Eigen::VectorXd b_mean,
                                            Eigen::MatrixXd Omega,
                                            Eigen::VectorXd Sigma_R,
                                            int n_random_effect,
                                            int n_choice_set)
{
  
  // Determine some dimensions/aspects of the design problem from the inputs
  int n_parameter         = b_mean.size();
  int n_fixed_effect      = n_parameter - n_random_effect;
  int n_alternative_total = X.rows();
  int choice_set_size     = n_alternative_total / n_choice_set;
  int n_samples           = Yall.cols();
  int n_samples_converged = Yall.cols(); // actual number of valid samples we have to work with
  
  
  // Inverse of random effects variance and corresponding subblock (columns)
  // of the model matrix (XR)
  Eigen::MatrixXd OmegaRE = Omega.block(0, 0, n_random_effect, n_random_effect);
  Eigen::MatrixXd XR = X.block(0, 0, X.rows(), n_random_effect);
  
  // Allocate space for terms in maximization
  Eigen::VectorXd u0_mle     = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd mstart     = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd u1_mle     = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd ustart_mle = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd Epy        = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd Eu2        = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd numerator  = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd numerator2 = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd prob       = Eigen::VectorXd::Zero(n_alternative_total);
  Eigen::VectorXd norm_prob  = Eigen::VectorXd::Zero(n_alternative_total);
  Eigen::MatrixXd delta      = Eigen::MatrixXd::Zero(n_alternative_total, n_alternative_total);
  Eigen::VectorXd gradient   = Eigen::VectorXd::Zero(X.cols());
  Eigen::VectorXd Ej         = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::MatrixXd H          = Eigen::MatrixXd::Zero(n_random_effect, n_random_effect);
  Eigen::MatrixXd Hinv       = Eigen::MatrixXd::Zero(n_random_effect, n_random_effect);
  Eigen::MatrixXd Hd         = Eigen::MatrixXd::Zero(n_random_effect, n_random_effect); // stores inverse H
  
  // Intermediate terms related to maximization
  Eigen::VectorXd SigmaRE          = Sigma_R.segment(0, n_random_effect);
  Eigen::VectorXd inv_Sigma_R_term = (2.0 * 3.14159265 * SigmaRE);
  inv_Sigma_R_term                 = inv_Sigma_R_term.array().rsqrt(); // reciprical square root
  Eigen::VectorXd exp_u_ratio_term = Eigen::VectorXd::Zero(n_random_effect);
  
  // Allocate space for a single sample (temporary variable)
  Eigen::VectorXd Y = Eigen::VectorXd::Zero(Yall.rows());
  
  // Rename the mean as the coefficient
  Eigen::VectorXd beta   = b_mean;
  
  // Model Matrix for a given choice set
  Eigen::MatrixXd X_cs      = Eigen::MatrixXd::Zero(choice_set_size, X.cols());
  Eigen::MatrixXd X_cs_rdm  = Eigen::MatrixXd::Zero(choice_set_size, n_random_effect);
  
  // Linear predictor for given choice set
  Eigen::VectorXd XB       = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd expXB    = Eigen::VectorXd::Zero(X.rows());
  
  // Variables related to inverse of variance-covariance
  Eigen::VectorXd score1        = Eigen::VectorXd::Zero(n_parameter);
  Eigen::VectorXd score2        = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::MatrixXd I11           = Eigen::MatrixXd::Zero(n_parameter, n_parameter);
  Eigen::MatrixXd I12           = Eigen::MatrixXd::Zero(n_parameter, n_random_effect);
  Eigen::MatrixXd I22           = Eigen::MatrixXd::Zero(n_random_effect, n_random_effect);
  
  // Scalar terms related to maximization
  double Hd_det   = 0.0;
  double denom    = 0.0;
  
  // Terms related to algorithm
  bool converged    = false;
  double tol        = 0.000000001;
  double stepsize   = 0.0;
  int max_iter      = 100;
  double addc       = 10000.0;
  double mleterm    = 0.0;
  double iterchange = 0.0;
  int i1            = 0;
  
  // stores p for a choice set, used in coverged = true loop
  Eigen::VectorXd pj = Eigen::VectorXd::Zero(choice_set_size);
  
  // Loop over all samples of Y
  for (int iy = 0; iy < n_samples; iy++){
    
    // Isolate Y, starting value for u
    Y         = Yall.col(iy);
    u0_mle    = u.col(iy);
    converged = false;
    
    // Find u1 mle via newton's method
    for(int inewton_iter = 0; inewton_iter < max_iter; inewton_iter++){
      
      beta.segment(0, n_random_effect)              = b_mean.segment(0, n_random_effect) + u0_mle;
      beta.segment(n_random_effect, n_fixed_effect) = b_mean.segment(n_random_effect, n_fixed_effect);
      
      // Calculate response probabilities under U0_MLE
      prob = calculate_response_probabilities(X, beta, n_choice_set);
      
      // Calculate diag(p) - pp' under U0_MLE
      delta = calculate_delta(prob, n_choice_set);
      
      // Corresponding gradient
      gradient = XR.transpose() * (Y - prob) -  OmegaRE * u0_mle;
      
      // Hess
      H = -(XR.transpose() * delta * XR) - OmegaRE;
      
      //stepsize = 1.0 / (1.0 + 10.0 * sqrt(gradient.array().pow(2).sum() ) );
      stepsize = 1.0;
      
      // Inverse Hess
      Hinv = H.inverse();
      
      // Increment mle
      u1_mle = u0_mle - stepsize*(Hinv * gradient);
      
      // Evaluate change at this iteration
      iterchange = 0.0;
      for (int ib = 0; ib < u1_mle.size(); ib++){
        iterchange += sqrt(pow(u1_mle(ib) - u0_mle(ib), 2));
      }
      
      // Check if converged
      if ( iterchange < tol ){
        converged = true;
        inewton_iter = max_iter + 1;
        u0_mle = u1_mle;
      } else {
        u0_mle = u1_mle;
      }
      
    } // end of newton iterations
    
    // If newton's algorithm didn't converge, move on to next sample
    if (converged == true){
      
      // Negative inv hess and corresponding determinant
      Hd     = -H.inverse();
      Hd_det = Hd.determinant();
      
      // u^2 / sigma
      exp_u_ratio_term = (
        -1.0*(u0_mle.array().pow(2)) * 
          (2.0 * SigmaRE.array()).inverse().array()
      ).exp();
      
      // Denominator
      denom = sqrt(Hd_det) *
        prob.array().pow(Y.array()).prod() * 
        (inv_Sigma_R_term.array() * exp_u_ratio_term.array()).prod();
      
      // The next nested loops cycle through each alternative one at a time
      // goal is to obtain numerator needed to evaluate E_u[p | y^(i)]
      mstart = u0_mle;
      for (int j1 = 0; j1 < X.rows(); j1++){
        
        // Start at u0_mle from first newton application
        u0_mle = mstart;
        
        // Isolate information from this choice set
        i1       = (j1 / choice_set_size); // indexes which choice set this alternative is a member of
        X_cs     = X.block(i1*choice_set_size, 0, choice_set_size, n_parameter);
        X_cs_rdm = X.block(i1*choice_set_size, 0, choice_set_size, n_random_effect);
        
        // Inner newton, calculate the numerator for ratio in laplace approximation
        // (not actually trying to get a new u here)
        for (int i_n = 0; i_n < max_iter; i_n++){
          
          beta.segment(0, n_random_effect)              = b_mean.segment(0, n_random_effect) + u0_mle;
          beta.segment(n_random_effect, n_fixed_effect) = b_mean.segment(n_random_effect, n_fixed_effect);
          
          // Response probabilities at current version of u
          prob = calculate_response_probabilities(X, beta, n_choice_set);
          
          // diag(p) - pp' at current version of u
          delta = calculate_delta(prob, n_choice_set);
          
          // Extract probabilities corresponding to current choice set
          pj = prob.block(i1*choice_set_size, 0, choice_set_size, 1).array();
          pj(j1 % choice_set_size) = -1.0 + pj(j1 % choice_set_size);
          pj *= -1.0;
          
          // Evaluate grad
          gradient = X_cs_rdm.transpose() * pj +
            XR.transpose() * (Y - prob) -
            OmegaRE * u0_mle;
          
          // Evaluate hess and its inverse
          H = -(
            X_cs_rdm.transpose() * 
              delta.block(i1*choice_set_size, i1*choice_set_size,
                          choice_set_size, choice_set_size) *
                            X_cs_rdm
          ) - (XR.transpose() * delta * XR ) - OmegaRE;
          Hinv = H.inverse();
          
          //stepsize = 1 / (1 + 10 * sqrt(gradient.array().pow(2).sum() ) );
          stepsize = 1.0;
          
          // Incr u0_mle
          u1_mle = u0_mle - stepsize*(Hinv * gradient);
          
          // Evaluate change at this iteration
          iterchange = 0.0;
          for (int ib = 0; ib < u1_mle.size(); ib++){
            iterchange += sqrt(pow(u1_mle(ib) - u0_mle(ib), 2));
          }
          
          // Check for convergence
          if ( iterchange < tol ){
            i_n = max_iter + 1;
            u0_mle = u1_mle;
          } else {
            u0_mle = u1_mle;
          }
        } // end of inner (converged) newton loop
        
        Hd     = -H.inverse();
        Hd_det = Hd.determinant();
        
        // Numerator for alternative j1
        exp_u_ratio_term = (-1.0*(u0_mle.array().pow(2)) * (2.0 * SigmaRE.array()).inverse()).exp();
        numerator(j1) = sqrt(Hd_det) * prob(j1) *
          prob.array().pow(Y.array()).prod() *
          (inv_Sigma_R_term.array() * exp_u_ratio_term.array()).prod();
        
      } // end of loop over individual alternatives (j1)
      
      // Second version, calculating numerator for random effect terms
      // Goal is to obtain numerator needed to evaluate E_u[u^3 / sigma^3 | y^(i)]
      // requires u_1j, for j indexing random effects
      for (int irandeff = 0; irandeff < n_random_effect; irandeff++){
        
        // Start at previous maximizer
        u0_mle = mstart; 
        
        // Zero out Ej and then isolate the current effect (1)
        Ej           = Eigen::VectorXd::Zero(Ej.size());
        Ej(irandeff) = 1.0;
        
        // Find u_1j, the maximizer of
        // log( (u_ij^2+csigmaj^3) / (sigmaj^3) ) + log P(Y = y^(i) | u) + log f_sigma(u)
        for (int newton_iter_2 = 0; newton_iter_2 < max_iter; newton_iter_2++){
          
          // Set beta
          beta.segment(0, n_random_effect)              = b_mean.segment(0, n_random_effect) + u0_mle;
          beta.segment(n_random_effect, n_fixed_effect) = b_mean.segment(n_random_effect, n_fixed_effect);
          
          // Calculate probabilities at beta
          prob = calculate_response_probabilities(X, beta, n_choice_set);
          
          // Calculate diag(p) - pp' at this value
          delta = calculate_delta(prob, n_choice_set);
          
          // Evaluate gradient
          gradient = (XR.transpose() * (Y - prob)) - 
            OmegaRE * u0_mle;
          
          // Construct part of term (intermediate step)
          mleterm = (
            2.0 * u0_mle(irandeff) / 
              (
                  pow(u0_mle(irandeff), 2) +
                    addc * pow(sqrt(SigmaRE(irandeff)), 3) 
              ) 
          ) ;
          
          // update gradient
          gradient(irandeff) += mleterm;//(mleterm * Ej.array()).matrix();
          
          // finalize calculation of mle term
          mleterm = 2.0 * (
            addc * pow(sqrt(Sigma_R(irandeff)), 3) - pow(u0_mle(irandeff), 2)
          ) /
            pow( pow(u0_mle(irandeff), 2) + addc * pow(sqrt(SigmaRE(irandeff)), 3), 2);
          
          // Evaluate hess (partial)
          H = Ej * Ej.transpose();
          for (int iel = 0; iel < n_random_effect; iel++){
            for (int iel2 = 0; iel2 < n_random_effect; iel2++){
              H(iel, iel2) *= mleterm;
            }
          }
          // Finalize
          H    = H - XR.transpose() * delta * XR - OmegaRE;
          
          // Invert
          Hinv = H.inverse();
          
          //stepsize = 1 / (1 + 10 * sqrt(gradient.array().pow(2).sum()) );
          stepsize = 1.0;
          
          // Increment u
          u1_mle = u0_mle - stepsize*(Hinv * gradient);
          
          // Calculate change
          iterchange = 0.0;
          for (int ib = 0; ib < u1_mle.size(); ib++){
            iterchange += sqrt(pow(u1_mle(ib) - u0_mle(ib), 2));
          }
          
          // Check for convergence
          if ( iterchange < tol ){
            newton_iter_2 = max_iter + 1;
            u0_mle = u1_mle;
          } else {
            u0_mle = u1_mle;
          }
          
        } // end of newton loop (inner, random effect)
        
        // Hd
        Hd     = -H.inverse();
        Hd_det = Hd.determinant();
        
        // Evaluate the second numerator needed for laplace approximation 
        exp_u_ratio_term = (-1.0*(u0_mle.array().pow(2)) * (2.0 * SigmaRE.array()).inverse()).exp();
        numerator2(irandeff) = sqrt(Hd_det) *
          (
              (pow(u0_mle(irandeff), 2) + addc * pow(sqrt(SigmaRE(irandeff)), 3)) /
                pow(sqrt(SigmaRE(irandeff)), 3)
          ) * prob.array().pow(Y.array()).prod() *
            (
                inv_Sigma_R_term.array() * exp_u_ratio_term.array()
            ).prod();
        
        
      } // End of random effect loop for Laplace numerator calculation
      
      // Ratio for first expectation, E_u[p | y^(i)]
      for (int i = 0; i < X.rows(); i++){
        Epy(i) = numerator(i) / denom;
      }
      
      // Now need to normalize the probabilities
      norm_prob = normalize_probabilities(Epy, n_choice_set);
      
      // Ratio for second expectation, E_u[u^3 / sigma^3 | y^(i)]
      for (int i = 0; i < n_random_effect; i++){
        Eu2(i) = numerator2(i) / denom - addc;
        if (Eu2(i) < 0.0){
          Eu2(i) *= -1.0;
        }
      }
      
      // Evaluate score terms
      score1 = X.transpose() * (Y - norm_prob);
      score2 = -(SigmaRE.array().rsqrt()) + Eu2.array(); 
      
      // Increment inverse of variance covariance matrix
      I11 += score1 * score1.transpose();
      I12 += score1 * score2.transpose();
      I22 += score2 * score2.transpose();
      
    } else { // case where outer newton did not converge, decrease n_samples_converged
      n_samples_converged -= 1;
    }
    
    
  } // end of loop over Y
  
  
  // Store the final inverse variance-covariance matrix
  Eigen::MatrixXd vcov_inverse = Eigen::MatrixXd::Zero(n_parameter + n_random_effect,
                                                       n_parameter + n_random_effect);
  
  vcov_inverse.block(0, 0, n_parameter, n_parameter)                             = I11;
  vcov_inverse.block(0, n_parameter, n_parameter, n_random_effect)               = I12;
  vcov_inverse.block(n_parameter, 0, n_random_effect, n_parameter)               = I12.transpose();
  vcov_inverse.block(n_parameter, n_parameter, n_random_effect, n_random_effect) = I22;
  
  // Normalize by number of replicates of Y
  for (int i = 0; i < n_parameter + n_random_effect; i++){
    for (int j = 0; j < n_parameter + n_random_effect; j++){
      vcov_inverse(i,j) *= 1.0 / n_samples_converged;
    }
  }
  
  // Return inverse of variance-covariance based on laplace approximation
  return(vcov_inverse);
}




//' Laplace approximation to the information matrix
 //'
 //' @description Calculates the variance-covariance matrix for the model parameters 
 //' under the Laplace Approximation. This is not intended to be directly used,
 //' instead a user should use the \code{\link{varcov_approx_PML}} function with the "Laplace" option.
 //'
 //' @param n_samples Number of response samples to draw
 //' @param X The model matrix
 //' @param b_mean The vector of effect means
 //' @param var_vec The vector of effect variances
 //' @param n_choice_set The number of choice sets in the discrete choice experiment
 //' @return A matrix containing the variance covariance matrix of the model parameters based on the Laplace approximation
 //' @keywords internal
 // [[Rcpp::export]]
 Eigen::MatrixXd Laplace_approx_vcov(int n_samples,
                                     Eigen::MatrixXd X,
                                     Eigen::VectorXd b_mean,
                                     Eigen::VectorXd var_vec,
                                     int n_choice_set){
   
   // Determine some aspects of the design problem from the inputs
   int n_alternative_total  = X.rows(); // total number of alternatives
   int n_beta               = b_mean.size();
   int choice_set_size      = n_alternative_total / n_choice_set;
   
   // Track how many effects are random
   int n_random_effect = 0;
   for (int i = 0; i < var_vec.size(); i++){
     if (var_vec(i) > 0){
       n_random_effect++;
     }
   }
   
   // How many effects are fixed
   int n_fixed_effect = var_vec.size() - n_random_effect;
   
   // Allocate space for final precision about the model parameters
   Eigen::MatrixXd vcov_inverse = Eigen::MatrixXd::Zero(n_beta + n_random_effect,
                                                        n_beta + n_random_effect);
   
   
   // Inverse of prior covariance
   Eigen::MatrixXd Omega = var_vec.array().inverse().matrix().asDiagonal();
   
   Eigen::VectorXd Sigma_R = var_vec;
   
   // Sample betas for each draw of Y
   Eigen::MatrixXd beta   = Eigen::MatrixXd::Zero(b_mean.size(), n_samples); // each column is a sample of beta
   Eigen::MatrixXd u      = Eigen::MatrixXd::Zero(n_random_effect, n_samples);
   Eigen::VectorXd sd_vec = Eigen::VectorXd::Zero(n_random_effect);
   
   // Draw samples of u
   for (int ip = 0; ip < n_random_effect; ip++){
     sd_vec(ip) = sqrt(var_vec(ip));
     for (int iu = 0; iu < n_samples; iu++){
       u(ip, iu) = R::rnorm(0.0, sd_vec(ip));
     }
   }
   Eigen::MatrixXd sd_mat = sd_vec.asDiagonal();
   
   // Get the corresponding samples of beta (b + u)
   for (int iy = 0; iy < n_samples; iy++){
     beta.block(0, iy, n_random_effect, 1)              = b_mean.segment(0, n_random_effect) + u.col(iy);
     beta.block(n_random_effect, iy, n_fixed_effect, 1) = b_mean.segment(n_random_effect, n_fixed_effect);
   }
   
   // Calculate the epsilon terms (get added to linear predictor later)
   Eigen::MatrixXd epsilon = Eigen::MatrixXd::Zero(X.rows(), n_samples);
   for (int iy = 0; iy < n_samples; iy++){
     epsilon.col(iy) = sample_epsilon(X.rows());
   }
   
   // Allocate intermediate quantities for drawing samples of the outcome vector
   Eigen::VectorXd XB    = Eigen::VectorXd::Zero(X.rows());
   Eigen::VectorXd expXB = Eigen::VectorXd::Zero(X.rows());
   Eigen::VectorXd prob  = Eigen::VectorXd::Zero(n_alternative_total);
   Eigen::MatrixXd Y     = Eigen::MatrixXd::Zero(X.rows(), n_samples);
   int yind              = 0; // for tracking which y is selected
   
   // Loop over n_samples
   for (int iy = 0; iy < n_samples; iy++){
     
     // Evaluate linear predictor
     XB    = X * beta.col(iy) + epsilon.col(iy);
     expXB = XB.array().exp();
     
     // As a starting value, pick the most likely Y for each choice set
     for (int ichoice = 0; ichoice < n_choice_set; ichoice++){
       yind=0;
       // Find the selected value for this choice set
       for (int ialt = 1; ialt < choice_set_size; ialt++){
         if (expXB(ichoice*choice_set_size + ialt) > expXB(ichoice*choice_set_size + yind)) {
           yind = ialt;
         }
       }
       Y(ichoice*choice_set_size + yind, iy) = 1;
     }
   } // /iy
   
   // Nested optimizations to approximate inverse variance covariance
   vcov_inverse = laplace_apply_newton_method(Y,
                                              X,
                                              u,
                                              b_mean,
                                              Omega,
                                              Sigma_R,
                                              n_random_effect,
                                              n_choice_set);
   
   
   // Invert to get variance-covariance matrix
   vcov_inverse = vcov_inverse.inverse();
   
   return vcov_inverse;
 }
