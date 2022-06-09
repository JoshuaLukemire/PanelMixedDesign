// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "approxutils.h"

// [[Rcpp::depends(RcppEigen)]]

// Function to get the matrix of difference codings
Eigen::MatrixXd getDifferenceCodings(Eigen::MatrixXd X,
                                     int nChoiceSet,
                                     int choiceSetSize){
  
  Eigen::MatrixXd XDC = Eigen::MatrixXd::Zero(nChoiceSet * (choiceSetSize - 1), X.cols());
  
  for (int i = 0; i < nChoiceSet; i++){
    for (int j = 0; j < choiceSetSize - 1; j++){
      XDC.row(i*(choiceSetSize-1) + j) = X.row(i*choiceSetSize + j) - X.row( (i+1)*choiceSetSize - 1);
    }
  }
  
  return XDC;
  
}



//' Approximation to the information matrix using PQL
//'
//' @description Approximates the information matrix for the model parameters 
//' using the partial quasi likelihood. This is not intended to be directly used,
//' instead a user should use the \code{\link{PMLInfoApprox}} function with 
//' the "PQL" option.
//'
//' @param nu Number of samples to draw
//' @param X The model matrix
//' @param b_mean The vector of effect means
//' @param var_vec The vector of effect variances
//' @param nChoiceSet The number of choice sets in the discrete choice experiment
//' @return A matrix containing the PQL approximation to the variance
//' covariance matrix of the model parameters
// [[Rcpp::export]]
Eigen::MatrixXd PQLApprox(int nu,
                          Eigen::MatrixXd X,
                          Eigen::VectorXd b_mean,
                          Eigen::VectorXd var_vec,
                          int nChoiceSet){
  
  // Determine some aspects of the design problem from the inputs
  int nAlternative  = X.rows();
  int nBeta         = b_mean.size();
  int choiceSetSize = nAlternative / nChoiceSet;
  
  int nRandEffect = 0;
  for (int i = 0; i < var_vec.size(); i++){
    if (var_vec(i) > 0){
      nRandEffect++;
    }
  }
  int nFixedEffect = var_vec.size() - nRandEffect;
  
  Eigen::MatrixXd fisherInfo = Eigen::MatrixXd::Zero( nBeta + nRandEffect,
                                                      nBeta + nRandEffect );
  
  // Inverse of prior covariance
  Eigen::MatrixXd Omega = var_vec.array().inverse().matrix().asDiagonal();
  
  // Precision matrix only for random effects 
  Eigen::MatrixXd OmegaRE = Omega.block(0, 0, nRandEffect, nRandEffect);
  
  Eigen::VectorXd SigmaR = var_vec;
  
  // Sample betas for each draw of Y
  std::random_device rd; 
  std::mt19937 gen(rd()); 
  std::normal_distribution<double> distribution(0.0, 1.0);
  
  Eigen::MatrixXd beta_start = Eigen::MatrixXd::Zero(b_mean.size(), nu); // each column is a sample of beta
  Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nRandEffect, nu);
  
  Eigen::VectorXd sd_vec = Eigen::VectorXd::Zero(nRandEffect);

  for (int ip = 0; ip < nRandEffect; ip++){
    sd_vec(ip) = sqrt(var_vec(ip));
    for (int iu = 0; iu < nu; iu++){
      u(ip, iu) = distribution(gen);
    }
  }
  
  Eigen::MatrixXd sd_mat = sd_vec.asDiagonal();
  
  for (int iy = 0; iy < nu; iy++){
    beta_start.block(0, iy, nRandEffect, 1) = b_mean.segment(0, nRandEffect) + sd_mat *  u.col(iy);
    beta_start.block(nRandEffect, iy, nFixedEffect, 1) = b_mean.segment(nRandEffect, nFixedEffect);
  }
  
  // Calculate the epsilon terms (get added to linear predictor later)
  Eigen::MatrixXd epsilon = Eigen::MatrixXd::Zero(X.rows(), nu);
  for (int iy = 0; iy < nu; iy++){
    epsilon.col(iy) = sampleEpsilon(nAlternative);
  }
  
  //
  // Calculate the linear predictor for each Y, draw the most likely Y
  Eigen::VectorXd XB    = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd expXB = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd prob  = Eigen::VectorXd::Zero(nAlternative);
  Eigen::MatrixXd Yall  = Eigen::MatrixXd::Zero(X.rows(), nu);
  int yind              = 0; // for tracking which y is selected
  
  Eigen::MatrixXd I11           = Eigen::MatrixXd::Zero(nBeta, nBeta);
  Eigen::MatrixXd I11Temp       = Eigen::MatrixXd::Zero(nBeta, nBeta);
  Eigen::MatrixXd I12           = Eigen::MatrixXd::Zero(nBeta, nRandEffect);
  Eigen::MatrixXd I22           = Eigen::MatrixXd::Zero(nBeta, nBeta);
  Eigen::MatrixXd beta          = Eigen::MatrixXd::Zero(nBeta, 1);
  
  // Loop over nu
  for (int iy = 0; iy < nu; iy++){
    XB    = X * beta_start.col(iy) + epsilon.col(iy);
    expXB = XB.array().exp();
    
    // As a starting value, pick the most likely Y for each choice set
    for (int ichoice = 0; ichoice < nChoiceSet; ichoice++){
      yind=0;
      // Find the selected value for this choice set
      for (int ialt = 1; ialt < choiceSetSize; ialt++){
        if (expXB(ichoice*choiceSetSize + ialt) > expXB(ichoice*choiceSetSize + yind)) {
          yind = ialt;
        }
      }
      Yall(ichoice*choiceSetSize + yind, iy) = 1;
    }
  }
  
  // Get the model matrix for differences
  Eigen::MatrixXd XDX =  getDifferenceCodings(X, nChoiceSet, choiceSetSize);

  // design matrix only corresponding to the random effects.
  Eigen::MatrixXd Xr = X.block(0, 0, X.rows(), nRandEffect);
  
  // Initialize quantities
  Eigen::VectorXd Y = Eigen::VectorXd::Zero(X.rows());
  bool converged = false;
  Eigen::VectorXd u0_mle     = Eigen::VectorXd::Zero(nRandEffect);
  Eigen::VectorXd u1_mle     = Eigen::VectorXd::Zero(nRandEffect);
  int max_iter = 100;
  Eigen::MatrixXd delta     = Eigen::MatrixXd::Zero(nAlternative, nAlternative);
  Eigen::MatrixXd deltaDX   = Eigen::MatrixXd::Zero(nChoiceSet * (choiceSetSize - 1),
                                                    nChoiceSet * (choiceSetSize - 1));
  Eigen::MatrixXd deltaDXInv   = Eigen::MatrixXd::Zero(nChoiceSet * (choiceSetSize - 1),
                                                       nChoiceSet * (choiceSetSize - 1));
  Eigen::MatrixXd Vn   = Eigen::MatrixXd::Zero(nChoiceSet * (choiceSetSize - 1),
                                               nChoiceSet * (choiceSetSize - 1));
  Eigen::MatrixXd VnInv   = Eigen::MatrixXd::Zero(nChoiceSet * (choiceSetSize - 1),
                                                  nChoiceSet * (choiceSetSize - 1));
  Eigen::VectorXd gradient  = Eigen::VectorXd::Zero(nRandEffect);
  Eigen::MatrixXd H         = Eigen::MatrixXd::Zero(nRandEffect, nRandEffect);
  Eigen::MatrixXd Hinv      = Eigen::MatrixXd::Zero(nRandEffect, nRandEffect);
  double stepsize = 0.0;
  double tol = 0.000000001;
  int nu_converged = nu;
  
  double iterchange = 0.0; // track convergence 
  
  // Loop over samples of Y/U
  for (int iy = 0; iy < nu; iy++){
    
    Y         = Yall.col(iy);
    u0_mle    = u.col(iy);
    converged = false;
    
    //std::cout << "========" << std::endl;
    
    for(int inewton_iter = 0; inewton_iter < max_iter; inewton_iter++){
      
      // Update beta using current value for U0
      beta.block(0, 0, nRandEffect, 1)            = b_mean.segment(0, nRandEffect) + u0_mle;
      beta.block(nRandEffect, 0, nFixedEffect, 1) = b_mean.segment(nRandEffect, nFixedEffect);
      
      
      // Calculate response probabilities under U0_MLE
      prob = calculateResponseProbs(X,
                                    beta,
                                    nChoiceSet);
      
      delta = calcDelta(prob, nChoiceSet);
      
      
      // from all rand effects verison
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
      
      //std::cout << iterchange << std::endl;
      
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
        deltaDX.block(ics*(choiceSetSize-1), ics*(choiceSetSize-1),
                      choiceSetSize-1, choiceSetSize-1) = delta.block(ics*choiceSetSize, ics*choiceSetSize,
                      choiceSetSize-1, choiceSetSize-1);
      }
      
      deltaDXInv = deltaDX.inverse();
      
      Vn = deltaDXInv + XDX * var_vec.matrix().asDiagonal() * XDX.transpose();
      VnInv = Vn.inverse();
      
      I11Temp = XDX.transpose() * VnInv * XDX;
      I11 += I11Temp;
      I22 += (4.0 * SigmaR.array()).sqrt().matrix().asDiagonal() * 
        I11Temp.array().pow(2).matrix() * 
        SigmaR.array().sqrt().matrix().asDiagonal();
      
      
    } else { // case where outer newton did not converge, decrease nu_converged counter
      nu_converged -= 1;
    }
    
    
  } // end of loop over Y
  
  //std::cout << nu_converged << std::endl;
  
  
  fisherInfo.block(0, 0, nBeta, nBeta) = I11;
  fisherInfo.block(0, nBeta, nBeta, nRandEffect) = I12;
  fisherInfo.block(nBeta, 0, nRandEffect, nBeta) = I12.transpose();
  fisherInfo.block(nBeta, nBeta, nRandEffect, nRandEffect) = I22.block(0, 0, nRandEffect, nRandEffect);
  // Normalize by number of replicates of Y
  for (int i = 0; i < nBeta + nRandEffect; i++){
    for (int j = 0; j < nBeta + nRandEffect; j++){
      fisherInfo(i,j) *= 1.0 / nu_converged;
    }
  }
  
  
  return(fisherInfo);
  
}
