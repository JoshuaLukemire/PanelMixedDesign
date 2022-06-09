// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "approxutils.h"


// [[Rcpp::depends(RcppEigen)]]

Eigen::MatrixXd laplaceSampleY(Eigen::MatrixXd expXB,
                               Eigen::MatrixXd X,
                               int nChoiceSet,
                               int ny){
  
  int choiceSetSize = X.rows() / nChoiceSet;
  Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(X.rows(), ny);
  
  
  int ind = 0;
  double maxval = 0.0;
  for (int iy = 0; iy < ny; iy++){
    
    // Determine the most likely outcome for this sampled Y
    for (int ic = 0; ic < nChoiceSet; ic++){
      maxval = 0.0;
      ind = 0;
      for (int ialt = 0; ialt < choiceSetSize; ialt++){
        // Check if this is largest
        if (expXB( ic*choiceSetSize + ialt ) > maxval){
          maxval = expXB( ic*choiceSetSize + ialt );
          ind = ic*choiceSetSize + ialt;
        }
      }
      // Set the most likely Y for this choice set to 1
      Y(ind, iy) = 1;
    }
  }
  return Y;
}

Eigen::VectorXd laplaceSampleEpsilon(int totalAlternatives){
  Eigen::VectorXd eps = Eigen::VectorXd::Zero(totalAlternatives);
  GetRNGstate();
  for (int ii = 0; ii < totalAlternatives; ii++){
    eps(ii) = R::runif(0.0, 1.0);
  }
  PutRNGstate();
  eps = (eps.array().log() * (-1.0) ).log() * (-1.0);
  return eps;
}


//
// Newtons method to find MLE(u0)
//
Eigen::MatrixXd laplaceNewtonMethod(Eigen::MatrixXd Yall,
                                    Eigen::MatrixXd X,
                                    Eigen::MatrixXd u,
                                    Eigen::VectorXd b_mean,
                                    Eigen::MatrixXd Omega,
                                    Eigen::VectorXd Sigma_R,
                                    int nRandEffect,
                                    int n_choice_set)
{
  
  
  int n_parameter        = b_mean.size();
  int nFixedEffect = n_parameter - nRandEffect;
  int n_alt_total        = X.rows();
  int choice_set_size    = n_alt_total / n_choice_set;
  int ny = Yall.cols();
  int ny_converged = Yall.cols(); // actual number of valid samples we have to work with for fisher info
  
  Eigen::MatrixXd OmegaRE = Omega.block(0, 0, nRandEffect, nRandEffect);
  Eigen::MatrixXd XR = X.block(0, 0, X.rows(), nRandEffect);
  
  //
  Eigen::VectorXd u0_mle     = Eigen::VectorXd::Zero(nRandEffect);
  Eigen::VectorXd mstart     = Eigen::VectorXd::Zero(nRandEffect);
  Eigen::VectorXd u1_mle     = Eigen::VectorXd::Zero(nRandEffect);
  Eigen::VectorXd ustart_mle = Eigen::VectorXd::Zero(nRandEffect);
  Eigen::VectorXd Epy        = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd Eu2        = Eigen::VectorXd::Zero(nRandEffect);
  
  Eigen::VectorXd SigmaRE = Sigma_R.segment(0, nRandEffect);
  
  Eigen::VectorXd inv_Sigma_R_term = (2.0 * 3.14159265 * SigmaRE);
  inv_Sigma_R_term = inv_Sigma_R_term.array().rsqrt(); // reciprical square root
  Eigen::VectorXd exp_u_ratio_term = Eigen::VectorXd::Zero(nRandEffect);
  
  Eigen::VectorXd Y = Eigen::VectorXd::Zero(Yall.rows());
  
  Eigen::VectorXd beta   = b_mean;
  
  // Model Matrix for a given choice set
  Eigen::MatrixXd X_cs      = Eigen::MatrixXd::Zero(choice_set_size, X.cols());
  Eigen::MatrixXd X_cs_rdm = Eigen::MatrixXd::Zero(choice_set_size, nRandEffect);
  
  Eigen::VectorXd XB       = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd expXB    = Eigen::VectorXd::Zero(X.rows());
  
  
  Eigen::VectorXd numerator  = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd numerator2 = Eigen::VectorXd::Zero(nRandEffect);
  
  
  Eigen::VectorXd prob      = Eigen::VectorXd::Zero(n_alt_total);
  Eigen::VectorXd norm_prob = Eigen::VectorXd::Zero(n_alt_total);
  Eigen::MatrixXd delta     = Eigen::MatrixXd::Zero(n_alt_total, n_alt_total);
  Eigen::VectorXd gradient  = Eigen::VectorXd::Zero(X.cols());
  // TODO figure out exactly what Ej is
  Eigen::VectorXd Ej       = Eigen::VectorXd::Zero(nRandEffect);
  
  Eigen::MatrixXd H        = Eigen::MatrixXd::Zero(nRandEffect, nRandEffect);
  Eigen::MatrixXd Hinv     = Eigen::MatrixXd::Zero(nRandEffect, nRandEffect);
  Eigen::MatrixXd Hd       = Eigen::MatrixXd::Zero(nRandEffect, nRandEffect); // stores inverse H
  
  // Variables related to information matrix
  Eigen::VectorXd score1        = Eigen::VectorXd::Zero(n_parameter);
  Eigen::VectorXd score2        = Eigen::VectorXd::Zero(nRandEffect);
  Eigen::MatrixXd I11           = Eigen::MatrixXd::Zero(n_parameter, n_parameter);
  Eigen::MatrixXd I12           = Eigen::MatrixXd::Zero(n_parameter, nRandEffect);
  Eigen::MatrixXd I22           = Eigen::MatrixXd::Zero(nRandEffect, nRandEffect);
  
  double Hd_det = 0.0;
  
  double denom = 0.0;
  
  bool converged = false;
  double tol = 0.000000001;
  double stepsize = 0.0;
  //int newton_iter = 0;
  int max_iter = 100;
  
  double addc = 10000.0;
  
  
  // convergence loop quantites
  int i1 = 0;
  
  // stores p for a choice set, used in coverged = true loop
  Eigen::VectorXd pj = Eigen::VectorXd::Zero(choice_set_size);
  
  double mleterm = 0.0;
  
  double iterchange = 0.0;
  
  for (int iy = 0; iy < ny; iy++){
    
    Y         = Yall.col(iy);
    u0_mle    = u.col(iy);
    converged = false;
    
    for(int inewton_iter = 0; inewton_iter < max_iter; inewton_iter++){
      
      //beta = b_mean + u0_mle;
      beta.segment(0, nRandEffect)            = b_mean.segment(0, nRandEffect) + u0_mle;
      beta.segment(nRandEffect, nFixedEffect) = b_mean.segment(nRandEffect, nFixedEffect);
      
      // Calculate response probabilities under U0_MLE
      prob = calculateResponseProbs(X,
                                    beta,
                                    n_choice_set);
      
      delta = calcDelta(prob,
                        n_choice_set);
      
      gradient = XR.transpose() * (Y - prob) -  OmegaRE * u0_mle;
      
      H = -(XR.transpose() * delta * XR) - OmegaRE; // this is Wei's version
      
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
        //u0_mle = u1_mle;
      } else {
        u0_mle = u1_mle;
      }
      
    } // end of newton iterations
    
    // If newton's algorithm convetrged, calculate information matrix 
    if (converged == true){
      
      
      Hd = -H.inverse();
      Hd_det = Hd.determinant();
      
      exp_u_ratio_term = (
        -1.0*(u0_mle.array().pow(2)) * 
          (2.0 * SigmaRE.array()).inverse().array()
      ).exp();
      
      
      denom = sqrt(Hd_det) *
        prob.array().pow(Y.array()).prod() * 
        (inv_Sigma_R_term.array() * exp_u_ratio_term.array()).prod();
      
      // The next nested loops cycle through each alternative and xxx
      mstart = u0_mle;
      
      for (int j1 = 0; j1 < X.rows(); j1++){
        
        u0_mle = mstart;
        
        i1        = (j1 / choice_set_size); // indexes which choice set this alternative is a member of
        X_cs      = X.block(i1*choice_set_size, 0, choice_set_size, n_parameter);
        X_cs_rdm = X.block(i1*choice_set_size, 0, choice_set_size, nRandEffect);
        
        
        for (int i_n = 0; i_n < max_iter; i_n++){
          
          beta.segment(0, nRandEffect)            = b_mean.segment(0, nRandEffect) + u0_mle;
          beta.segment(nRandEffect, nFixedEffect) =  b_mean.segment(nRandEffect, nFixedEffect);
          
          prob = calculateResponseProbs(X,
                                        beta,
                                        n_choice_set);
          
          delta = calcDelta(prob, n_choice_set);
          
          // Extract probabilities corresponding to current choice set
          pj = prob.block(i1*choice_set_size, 0, choice_set_size, 1).array();
          pj(j1 % choice_set_size) = -1.0 + pj(j1 % choice_set_size);
          pj *= -1.0;
          
          gradient = X_cs_rdm.transpose() * pj +
            XR.transpose() * (Y - prob) -
            OmegaRE * u0_mle;
          
          H = -(
            X_cs_rdm.transpose() * 
              delta.block(i1*choice_set_size, i1*choice_set_size,
                          choice_set_size, choice_set_size) *
                            X_cs_rdm
          ) - (XR.transpose() * delta * XR ) - OmegaRE;
          
          Hinv = H.inverse();
          
          //stepsize = 1 / (1 + 10 * sqrt(gradient.array().pow(2).sum() ) );
          stepsize = 1.0;
          
          u1_mle = u0_mle - stepsize*(Hinv * gradient);
          
          iterchange = 0.0;
          for (int ib = 0; ib < u1_mle.size(); ib++){
            iterchange += sqrt(pow(u1_mle(ib) - u0_mle(ib), 2));
          }
          
          if ( iterchange < tol ){
            //converged = true;
            i_n = max_iter + 1;
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
      
      for (int irandeff = 0; irandeff < nRandEffect; irandeff++){
        
        u0_mle = mstart; // TODO verify RHS is correct
        
        Ej = Eigen::VectorXd::Zero(Ej.size());
        Ej(irandeff) = 1.0;
        
        for (int newton_iter_2 = 0; newton_iter_2 < max_iter; newton_iter_2++){
          beta.segment(0, nRandEffect)            = b_mean.segment(0, nRandEffect) + u0_mle;
          beta.segment(nRandEffect, nFixedEffect) = b_mean.segment(nRandEffect, nFixedEffect);
          
          // Calculate probabilities
          prob = calculateResponseProbs(X,
                                        beta,
                                        n_choice_set);
          // Need to calculate delta
          delta = calcDelta(prob, n_choice_set);

          gradient = (XR.transpose() * (Y - prob)) - 
            OmegaRE * u0_mle;
          
          
          mleterm = (
            2.0 * u0_mle(irandeff) / 
              (
                  pow(u0_mle(irandeff), 2) +
                    addc * pow(sqrt(SigmaRE(irandeff)), 3) 
              ) 
          ) ;
          
          
          gradient(irandeff) += mleterm;//(mleterm * Ej.array()).matrix();
          
          
          mleterm = 2.0 * (
            addc * pow(sqrt(Sigma_R(irandeff)), 3) - pow(u0_mle(irandeff), 2)
          ) /
            pow( pow(u0_mle(irandeff), 2) + addc * pow(sqrt(SigmaRE(irandeff)), 3), 2);
          
          
          H = Ej * Ej.transpose();
          for (int iel = 0; iel < nRandEffect; iel++){
            for (int iel2 = 0; iel2 < nRandEffect; iel2++){
              H(iel, iel2) *= mleterm;
            }
          }
          
          
          H    = H - XR.transpose() * delta * XR - OmegaRE;// This is what wei had
          
          Hinv = H.inverse();
          
          //stepsize = 1 / (1 + 10 * sqrt(gradient.array().pow(2).sum()) );
          stepsize = 1.0;
          
          u1_mle = u0_mle - stepsize*(Hinv * gradient);
          
          
          iterchange = 0.0;
          for (int ib = 0; ib < u1_mle.size(); ib++){
            iterchange += sqrt(pow(u1_mle(ib) - u0_mle(ib), 2));
          }
          
          if ( iterchange < tol ){
            newton_iter_2 = max_iter + 1;
          } else {
            u0_mle = u1_mle;
          }
          
        } // end of yet another newton loop
        
        Hd     = -H.inverse();
        Hd_det = Hd.determinant();
        
        
        ////////////////
        
        
        exp_u_ratio_term = (-1.0*(u0_mle.array().pow(2)) * (2.0 * SigmaRE.array()).inverse()).exp();
        
        numerator2(irandeff) = sqrt(Hd_det) *
          (
              (pow(u0_mle(irandeff), 2) + addc * pow(sqrt(SigmaRE(irandeff)), 3)) /
                pow(sqrt(SigmaRE(irandeff)), 3)
          ) * prob.array().pow(Y.array()).prod() *
            (
                inv_Sigma_R_term.array() * exp_u_ratio_term.array()
            ).prod();
        
        
      } // This is the end of the irandeff loop
      
      for (int i = 0; i < X.rows(); i++){
        Epy(i) = numerator(i) / denom;
      }
      
      // Now need to normalize the probabilities
      norm_prob = normalizeProbabilities(Epy, n_choice_set);
      
      
      for (int i = 0; i < nRandEffect; i++){
        Eu2(i) = numerator2(i) / denom - addc;
        if (Eu2(i) < 0.0){
          Eu2(i) *= -1.0;
        }
      }
      
      
      score1 = X.transpose() * (Y - norm_prob);
      score2 = -(SigmaRE.array().rsqrt()) + Eu2.array(); 
      

      I11 += score1 * score1.transpose();
      I12 += score1 * score2.transpose();
      I22 += score2 * score2.transpose();
      
      
    } else { // case where outer newton did not converge, decrease ny_converged
      ny_converged -= 1;
    }
    
    
  } // end of loop over Y
  
  
  // Store the final information matrix
  Eigen::MatrixXd infoApprox = Eigen::MatrixXd::Zero(n_parameter + nRandEffect,
                                                     n_parameter + nRandEffect);
  
  infoApprox.block(0, 0, n_parameter, n_parameter) = I11;
  infoApprox.block(0, n_parameter, n_parameter, nRandEffect) = I12;
  infoApprox.block(n_parameter, 0, nRandEffect, n_parameter) = I12.transpose();
  infoApprox.block(n_parameter, n_parameter, nRandEffect, nRandEffect) = I22;
  // Normalize by number of replicates of Y
  for (int i = 0; i < n_parameter + nRandEffect; i++){
    for (int j = 0; j < n_parameter + nRandEffect; j++){
      infoApprox(i,j) *= 1.0 / ny_converged;
    }
  }
  
  return(infoApprox);
  
}




//' Laplace approximation to the information matrix
//'
//' @description Approximates the information matrix for the model parameters 
//' using a Laplace Approximation. This is not intended to be directly used,
//' instead a user should use the \code{\link{PMLInfoApprox}} function with the "Laplace" option.
//'
//' @param ny Number of response samples to draw
//' @param X The model matrix
//' @param b_mean The vector of effect means
//' @param var_vec The vector of effect variances
//' @param n_choice_set The number of choice sets in the discrete choice experiment
//' @return A matrix containing the Laplace approximation to the variance covariance matrix of the model parameters
// [[Rcpp::export]]
Eigen::MatrixXd laplaceApproximation(int ny,
                                     Eigen::MatrixXd X,
                                     Eigen::VectorXd b_mean,
                                     Eigen::VectorXd var_vec,
                                     int n_choice_set){
  
  
  int nRandEffect = 0;
  for (int i = 0; i < var_vec.size(); i++){
    if (var_vec(i) > 0){
      nRandEffect++;
    }
  }
  int nFixedEffect = var_vec.size() - nRandEffect;
  
  Eigen::MatrixXd infoApprox = Eigen::MatrixXd::Zero( b_mean.size() + nRandEffect,
                                                      b_mean.size() + nRandEffect );
  
  int n_alt_total = X.rows();
  int choice_set_size = n_alt_total / n_choice_set;
  
  // Inverse of prior covariance
  Eigen::MatrixXd Omega = var_vec.array().inverse().matrix().asDiagonal();
  
  Eigen::VectorXd Sigma_R = var_vec;
  
  // Sample betas for each draw of Y
  Eigen::MatrixXd beta = Eigen::MatrixXd::Zero(b_mean.size(), ny); // each column is a sample of beta
  Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nRandEffect, ny);
  Eigen::VectorXd sd_vec = Eigen::VectorXd::Zero(nRandEffect);
  
  GetRNGstate();
  for (int ip = 0; ip < nRandEffect; ip++){
    sd_vec(ip) = sqrt(var_vec(ip));
    for (int iu = 0; iu < ny; iu++){
      u(ip, iu) = norm_rand();
    }
  }
  PutRNGstate();
  
  
  Eigen::MatrixXd sd_mat = sd_vec.asDiagonal();
  
  
  
  for (int iy = 0; iy < ny; iy++){
    beta.block(0, iy, nRandEffect, 1)            = b_mean.segment(0, nRandEffect) + sd_mat *  u.col(iy);
    beta.block(nRandEffect, iy, nFixedEffect, 1) = b_mean.segment(nRandEffect, nFixedEffect);
  }
  
  // Calculate the epsilon terms (get added to linear predictor later)
  Eigen::MatrixXd epsilon = Eigen::MatrixXd::Zero(X.rows(), ny);
  for (int iy = 0; iy < ny; iy++){
    epsilon.col(iy) = laplaceSampleEpsilon(X.rows());
  }
  
  ///*
  
  //
  // Calculate the linear predictor for each Y, draw the most likely Y
  Eigen::VectorXd XB    = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd expXB = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd prob  = Eigen::VectorXd::Zero(n_alt_total);
  Eigen::MatrixXd Y     = Eigen::MatrixXd::Zero(X.rows(), ny);
  int yind              = 0; // for tracking which y is selected
  
  
  // Loop over ny
  for (int iy = 0; iy < ny; iy++){
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
  // */
  
  infoApprox = laplaceNewtonMethod(Y,
                                   X,
                                   u,
                                   b_mean,
                                   Omega,
                                   Sigma_R,
                                   nRandEffect,
                                   n_choice_set);
  
  
  
  return infoApprox;
  
}
