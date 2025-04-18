// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "approxutils.h"


// [[Rcpp::depends(RcppEigen)]]

//' Approximation to the variance-covariance matrix for the panel mixed logit model using the method of simulated moments
//'
//' @description Approximates the variance covariance matrix for the model parameters 
//' using the method of simulated moments. This is not intended to be directly used,
//' instead a user should use the \code{\link{varcov_approx_PML}} function with 
//' the "MSM" option.
//'
//' @param n_samples Number of response samples to draw
//' @param X The model matrix
//' @param b_mean The vector of effect means
//' @param var_vec The vector of effect variances
//' @param n_choice_set The number of choice sets in the discrete choice experiment
//' @return A matrix containing the MSM evaluated approximation of the variance
//' covariance matrix of the model parameters
 // [[Rcpp::export]]
Eigen::MatrixXd MSM_approx_vcov(int n_samples,
                          Eigen::MatrixXd X,
                          Eigen::VectorXd b_mean,
                          Eigen::VectorXd var_vec,
                          int n_choice_set){
  
  
  // Determine some dimensions based on input
  int n_alternative   = X.rows();
  int n_beta          = b_mean.size();
  int choice_set_size = n_alternative / n_choice_set;
  
  // Count how many effects are random
  int n_random_effect = 0;
  for (int i = 0; i < var_vec.size(); i++){
    if (var_vec(i) > 0){
      n_random_effect++;
    }
  }
  
  // Inverse of prior covariance and vector version of variance
  Eigen::MatrixXd Omega  = var_vec.array().inverse().matrix().asDiagonal();
  Eigen::VectorXd SigmaR = var_vec;
  
  // Sample betas (via u) for each draw of Y
  Eigen::MatrixXd beta   = Eigen::MatrixXd::Zero(b_mean.size(), n_samples); // each column is a sample of beta
  Eigen::VectorXd beta_i = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::MatrixXd u      = Eigen::MatrixXd::Zero(n_random_effect, n_samples);
  Eigen::VectorXd sd_vec = var_vec.array().sqrt();
  Eigen::MatrixXd sd_mat = sd_vec.asDiagonal();
  for (int ip = 0; ip < n_random_effect; ip++){
    for (int iu = 0; iu < n_samples; iu++){
      u(ip, iu) = R::rnorm(0.0, 1.0);
    }
  }
  for (int iy = 0; iy < n_samples; iy++){
    beta.col(iy) = b_mean + sd_mat *  u.col(iy);
  }
  
  // Calculate the epsilon terms (get added to linear predictor later when drawing Y)
  Eigen::MatrixXd epsilon = Eigen::MatrixXd::Zero(X.rows(), n_samples);
  for (int iy = 0; iy < n_samples; iy++){
    epsilon.col(iy) = sample_epsilon(n_alternative);
  }
  
  // Calculate the linear predictor for each Y, draw the most likely Y
  Eigen::VectorXd XB    = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd expXB = Eigen::VectorXd::Zero(X.rows());
  Eigen::MatrixXd prob  = Eigen::MatrixXd::Zero(n_alternative, n_samples);
  Eigen::MatrixXd Yall  = Eigen::MatrixXd::Zero(X.rows(), n_samples);
  int yind              = 0; // for tracking which y is selected
  
  // Allocate space for sub-blocks of inverse variance-covariance
  Eigen::MatrixXd I11           = Eigen::MatrixXd::Zero(n_beta, n_beta);
  Eigen::MatrixXd I11Temp       = Eigen::MatrixXd::Zero(n_beta, n_beta);
  Eigen::MatrixXd I12           = Eigen::MatrixXd::Zero(n_beta, n_random_effect);
  Eigen::MatrixXd I22           = Eigen::MatrixXd::Zero(n_random_effect, n_random_effect);
  
  // Allocate space for expectations in method of simulated moments
  Eigen::VectorXd ep    = Eigen::VectorXd::Zero(n_alternative);
  Eigen::MatrixXd epu   = Eigen::MatrixXd::Zero(n_alternative, b_mean.size());
  Eigen::VectorXd epp   = Eigen::VectorXd::Zero(pow(n_alternative, 2));
  Eigen::MatrixXd eppu  = Eigen::MatrixXd::Zero(pow(n_alternative, 2), b_mean.size());
  Eigen::VectorXd eppp  = Eigen::VectorXd::Zero(pow(n_alternative, 3));
  Eigen::MatrixXd epppu = Eigen::MatrixXd::Zero(pow(n_alternative, 3), b_mean.size());
  
  // Generate the "observed" outcome for each set of samples
  for (int iy = 0; iy < n_samples; iy++){
    // Linear predictor based on u
    XB    = X * beta.col(iy) + epsilon.col(iy);
    expXB = XB.array().exp();
    
    // Assign Y
    for (int ichoice = 0; ichoice < n_choice_set; ichoice++){
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
  
  // Random samples for u
  for (int ip = 0; ip < n_random_effect; ip++){
    for (int iu = 0; iu < n_samples; iu++){
      u(ip, iu) = R::rnorm(0.0, 1.0);
    }
  }
  
  // Loop over each sample
  for (int iy = 0; iy < n_samples; iy++){
    // Effects based on this draw of u
    beta_i = b_mean + sd_mat *  u.col(iy);
    
    // Evaluate response probabilities
    prob.col(iy) = calculate_response_probabilities(X, beta_i, n_choice_set);
  }
  
  // Loop over each sample
  for (int iy = 0; iy < n_samples; iy++){
    ep  += prob.col(iy);
    epu += prob.col(iy) * u.col(iy).transpose();
    // Double loop over all combinations of alternatives
    for (int ialt = 0; ialt < n_alternative; ialt++){
      for (int jalt = 0; jalt < n_alternative; jalt++){
        epp(ialt*n_alternative + jalt) += (prob(ialt, iy)*prob(jalt, iy));
        eppu.row(ialt*n_alternative + jalt) += ((prob(ialt, iy)*prob(jalt, iy)) * u.col(iy).array()).matrix();
        
        // Triple loop, cycle through alternatives again
        for (int kalt = 0; kalt < n_alternative; kalt++){
          eppp( (ialt*n_alternative + jalt)*n_alternative + kalt ) += (prob(ialt, iy)*prob(jalt, iy)*prob(kalt, iy));
          epppu.row( (ialt*n_alternative + jalt)*n_alternative + kalt ) += ((prob(ialt, iy)*prob(jalt, iy)*prob(kalt, iy)) * 
            u.col(iy).array()).matrix();
        }
      } // jalt
    } // ialt
  }
  
  // Turn the sums into averages to get expectation terms
  ep    = ep.array() / n_samples;
  epu   = epu.array() / n_samples;
  epp   = epp.array() / n_samples;
  eppu  = eppu.array() / n_samples;
  eppp  = eppp.array() / n_samples;
  epppu = epppu.array() / n_samples;
  
  // Allocate space for U terms
  Eigen::VectorXd U11temp = Eigen::VectorXd::Zero(b_mean.size());
  Eigen::MatrixXd U11     = Eigen::MatrixXd::Zero(n_alternative, b_mean.size());
  Eigen::VectorXd U12temp = Eigen::VectorXd::Zero(b_mean.size());
  Eigen::MatrixXd U12     = Eigen::MatrixXd::Zero(n_alternative, b_mean.size());
  Eigen::MatrixXd U21     = Eigen::MatrixXd::Zero( pow(n_alternative, 2), b_mean.size());
  Eigen::MatrixXd U22     = Eigen::MatrixXd::Zero( pow(n_alternative, 2), b_mean.size());
  Eigen::VectorXd U21temp = Eigen::VectorXd::Zero(b_mean.size());
  Eigen::VectorXd U22temp = Eigen::VectorXd::Zero(b_mean.size());
  
  // Tracking choice set and alternative rows
  int ics = 0;
  int ics2 = 0;
  int ind1 = 0;
  int ind2 = 0;
  
  // selected elements from ep vecs
  double eppp1 = 0.0;
  double eppp2 = 0.0;
  
  // Loop over alternatives
  for (int ialt = 0; ialt < n_alternative; ialt++){
    
    // Zero out temp variables
    U11temp = Eigen::VectorXd::Zero(b_mean.size());
    U12temp = Eigen::VectorXd::Zero(b_mean.size());
    
    // Increment choice set counter
    if (ialt == (ics+1)*choice_set_size){
      ics++;
    }
    
    // Calculate intermediate terms
    for (int jalt = 0; jalt < choice_set_size; jalt++){
      U11temp = (U11temp.transpose().array() + epp(ialt*n_alternative + ics*choice_set_size + jalt) * X.row(ics*choice_set_size + jalt).array()).transpose();
      for (int iparam = 0; iparam < b_mean.size(); iparam++){
        U12temp(iparam) = U12temp(iparam) + X(ics*choice_set_size + jalt, iparam) * eppu(ialt*n_alternative + ics*choice_set_size + jalt, iparam);
      } // iparam
    } // jalt
    
    // Evaluate U11
    U11.row(ialt) = ep(ialt) * X.row(ialt).array() - U11temp.transpose().array();
    
    // Evaluate U12
    for (int iparam = 0; iparam < b_mean.size(); iparam++){
      U12(ialt, iparam) = X(ialt, iparam) * epu(ialt, iparam) - U12temp(iparam);
    } // iparam
    
    
    ics2 = 0; // reset second choice set counter
    for (int ialt2 = 0; ialt2 < n_alternative; ialt2++){
      // Increment second choice set counter
      if (ialt2 == (ics2+1)*choice_set_size){
        ics2++;
      }
      
      // Case where comparing alternatives in the same choice set
      if (ics == ics2){
        if (ialt == ialt2){
          U21.row(ialt*n_alternative + ialt2) = U11.row(ialt);
          U22.row(ialt*n_alternative + ialt2) = U12.row(ialt);
        } else {
          U21.row(ialt*n_alternative + ialt2) = U21.row(ialt*n_alternative + ialt2).array() * 0.0;
          U22.row(ialt*n_alternative + ialt2) = U22.row(ialt*n_alternative + ialt2).array() * 0.0;
        }
        // Otherwise the choice sets don't match
      } else {
        // Zero out temporary variables
        U21temp = Eigen::VectorXd::Zero(b_mean.size());
        U22temp = Eigen::VectorXd::Zero(b_mean.size());
        // Loop over all alternatives
        for (int jalt = 0; jalt < choice_set_size; jalt++){
          ind1 = (ialt *n_alternative + ics * choice_set_size + jalt) * n_alternative + ialt2;
          ind2 = (ialt2*n_alternative + ics2* choice_set_size + jalt) * n_alternative + ialt;
          
          eppp1 = eppp( ind1 );
          eppp2 = eppp( ind2 );
          
          for (int iparam = 0; iparam < b_mean.size(); iparam++){
            U21temp(iparam) += eppp1 * X(ics*choice_set_size + jalt, iparam) + 
              eppp2 * X(ics2*choice_set_size + jalt, iparam);
          }
          
          // Loop over each parameter in the model
          for (int iparam = 0; iparam < b_mean.size(); iparam++){
            U22temp(iparam) = U22temp(iparam) + 
              X(ics*choice_set_size + jalt, iparam) * 
              epppu(ind1, iparam) +
              X(ics2*choice_set_size + jalt, iparam) * 
              epppu(ind2, iparam);
          } // iparam
          
        } // /jalt
        
        for (int iparam = 0; iparam < b_mean.size(); iparam++){
          U21(ialt*n_alternative + ialt2, iparam) = epp(ialt*n_alternative + ialt2, iparam) * 
            (X(ialt, iparam)+X(ialt2, iparam)) -
            U21temp(iparam);
          U22(ialt*n_alternative + ialt2, iparam) = (X(ialt, iparam) + X(ialt2, iparam)) * 
            eppu(ialt*n_alternative + ialt2, iparam) -
            U22temp(iparam);
        }
      } // else
    } // ialt2
  } // ialt
  
  Eigen::MatrixXd XX = Eigen::MatrixXd::Zero( pow(n_alternative, 2), b_mean.size());
  for (int iparam = 0; iparam < b_mean.size(); iparam++){
    for (int ialt = 0; ialt < n_alternative; ialt++){
      for (int jalt = 0; jalt < n_alternative; jalt++){
        XX(ialt*n_alternative + jalt, iparam) = X(ialt, iparam) * X(jalt, iparam);
      }
    }
  }
  
  Eigen::MatrixXd DU11 = X.transpose()  * U11.matrix();
  Eigen::MatrixXd DU12 = X.transpose()  * U12;
  Eigen::MatrixXd DU21 = XX.matrix().transpose() * U21.matrix();
  Eigen::MatrixXd DU22 = XX.matrix().transpose() * U22.matrix();
  
  Eigen::MatrixXd MME = Eigen::MatrixXd::Zero(2 * b_mean.size(), n_samples);
  MME.block(0, 0, b_mean.size(), n_samples) = X.transpose() * Yall;
  MME.block(b_mean.size(), 0, b_mean.size(), n_samples) = (X.transpose() * Yall).array().pow(2).matrix();
  
  Eigen::MatrixXd VM  = Eigen::MatrixXd::Zero(2 * b_mean.size(), 2 * b_mean.size());
  Eigen::MatrixXd VMi = Eigen::MatrixXd::Zero( b_mean.size() + n_random_effect, b_mean.size() + n_random_effect);
  Eigen::VectorXd EM  = Eigen::VectorXd::Zero(2 * b_mean.size());
  
  EM = MME.rowwise().sum().array() / n_samples;
  
  VM = ((MME * MME.transpose()).array() / n_samples).matrix() - EM * EM.transpose();
  VMi = VM.block(0, 0, b_mean.size() + n_random_effect, b_mean.size() + n_random_effect);
  
  Eigen::MatrixXd DU = Eigen::MatrixXd::Zero( n_beta + n_random_effect,
                                              n_beta + n_random_effect );
  
  DU.block(0, 0, n_beta, n_beta) = DU11;
  DU.block(0, n_beta, n_beta, n_random_effect) = DU12;
  DU.block(n_beta, 0, n_random_effect, n_beta) = DU21;
  DU.block(n_beta, n_beta, n_random_effect, n_random_effect) = DU22.block(0, 0, n_random_effect, n_random_effect);
  
  Eigen::MatrixXd stackedMatrix = Eigen::MatrixXd::Zero( 2*(b_mean.size() + n_random_effect), b_mean.size() + n_random_effect);
  stackedMatrix.block(0, 0, b_mean.size() + n_random_effect, b_mean.size() + n_random_effect) = DU;
  stackedMatrix.block(b_mean.size() + n_random_effect, 0, b_mean.size() + n_random_effect, b_mean.size() + n_random_effect) = VMi;
  
  return(stackedMatrix);
  
}
