// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "approxutils.h"


// [[Rcpp::depends(RcppEigen)]]

//' Approximation to the information matrix using the method of simulated moments
//'
//' @description Approximates the information matrix for the model parameters 
//' using the method of simulated moments. This is not intended to be directly used,
//' instead a user should use the \code{\link{PMLInfoApprox}} function with 
//' the "MSM" option.
//'
//' @param ny Number of response samples to draw
//' @param X The model matrix
//' @param b_mean The vector of effect means
//' @param var_vec The vector of effect variances
//' @param nChoiceSet The number of choice sets in the discrete choice experiment
//' @return A matrix containing the MSM evaluated of the variance
//' covariance matrix of the model parameters
// [[Rcpp::export]]
Eigen::MatrixXd MSMApprox(int ny,
                          Eigen::MatrixXd X,
                          Eigen::VectorXd b_mean,
                          Eigen::VectorXd var_vec,
                          int nChoiceSet){
  
  
  // Determine some aspects of the design problem from the inputs
  int nAlternative = X.rows();
  int nBeta        = b_mean.size();
  int choiceSetSize = nAlternative / nChoiceSet;
  
  int nRandEffect = 0;
  for (int i = 0; i < var_vec.size(); i++){
    if (var_vec(i) > 0){
      nRandEffect++;
    }
  }
  
  Eigen::MatrixXd fisherInfo = Eigen::MatrixXd::Zero( nBeta + nRandEffect,
                                                      nBeta + nRandEffect );
  
  // Inverse of prior covariance 
  Eigen::MatrixXd Omega = var_vec.array().inverse().matrix().asDiagonal();
  
  Eigen::VectorXd SigmaR = var_vec;
  
  // Sample betas for each draw of Y
  std::random_device rd; 
  std::mt19937 gen(rd()); 
  std::normal_distribution<double> distribution(0.0, 1.0);
  Eigen::MatrixXd beta = Eigen::MatrixXd::Zero(b_mean.size(), ny); // each column is a sample of beta
  Eigen::VectorXd beta_i = Eigen::VectorXd::Zero(nRandEffect);
  Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nRandEffect, ny);
  Eigen::VectorXd sd_vec = var_vec.array().sqrt();
  Eigen::MatrixXd sd_mat = sd_vec.asDiagonal();
  for (int ip = 0; ip < nRandEffect; ip++){
    for (int iu = 0; iu < ny; iu++){
      u(ip, iu) = distribution(gen);
    }
  }
  
  for (int iy = 0; iy < ny; iy++){
    beta.col(iy) = b_mean + sd_mat *  u.col(iy);
  }
  
  // Calculate the epsilon terms (get added to linear predictor later)
  Eigen::MatrixXd epsilon = Eigen::MatrixXd::Zero(X.rows(), ny);
  for (int iy = 0; iy < ny; iy++){
    epsilon.col(iy) = sampleEpsilon(nAlternative);
  }
  
  // Calculate the linear predictor for each Y, draw the most likely Y
  Eigen::VectorXd XB    = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd expXB = Eigen::VectorXd::Zero(X.rows());
  Eigen::MatrixXd prob  = Eigen::MatrixXd::Zero(nAlternative, ny);
  Eigen::MatrixXd Yall  = Eigen::MatrixXd::Zero(X.rows(), ny);
  int yind              = 0; // for tracking which y is selected
  
  Eigen::MatrixXd I11           = Eigen::MatrixXd::Zero(nBeta, nBeta);
  Eigen::MatrixXd I11Temp       = Eigen::MatrixXd::Zero(nBeta, nBeta);
  Eigen::MatrixXd I12           = Eigen::MatrixXd::Zero(nBeta, nRandEffect);
  Eigen::MatrixXd I22           = Eigen::MatrixXd::Zero(nRandEffect, nRandEffect);
  
  // MSM-specific terms
  Eigen::VectorXd ep = Eigen::VectorXd::Zero(nAlternative);
  Eigen::MatrixXd epu = Eigen::MatrixXd::Zero(nAlternative, b_mean.size());
  Eigen::VectorXd epp = Eigen::VectorXd::Zero(pow(nAlternative, 2));
  Eigen::MatrixXd eppu = Eigen::MatrixXd::Zero(pow(nAlternative, 2), b_mean.size());
  Eigen::VectorXd eppp = Eigen::VectorXd::Zero(pow(nAlternative, 3));
  Eigen::MatrixXd epppu = Eigen::MatrixXd::Zero(pow(nAlternative, 3), b_mean.size());
  
  // Loop over ny
  for (int iy = 0; iy < ny; iy++){
    XB    = X * beta.col(iy) + epsilon.col(iy);
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
  
  
  // New set of u's
  for (int ip = 0; ip < nRandEffect; ip++){
    for (int iu = 0; iu < ny; iu++){
      u(ip, iu) = distribution(gen);
    }
  }
  for (int iy = 0; iy < ny; iy++){
    

    beta_i = b_mean + sd_mat *  u.col(iy);
    
    prob.col(iy) = calculateResponseProbs(X,
             beta_i,
             nChoiceSet);
    
    
  }
  
  
  for (int iy = 0; iy < ny; iy++){
    ep  += prob.col(iy);
    epu += prob.col(iy) * u.col(iy).transpose();
    
    
    // Double loop over all combinations of alternatives
    for (int ialt = 0; ialt < nAlternative; ialt++){
      for (int jalt = 0; jalt < nAlternative; jalt++){
        epp(ialt*nAlternative + jalt) += (prob(ialt, iy)*prob(jalt, iy));
        eppu.row(ialt*nAlternative + jalt) += ((prob(ialt, iy)*prob(jalt, iy)) * u.col(iy).array()).matrix();
        
        // Triple loop, cycle through alternatives again
        for (int kalt = 0; kalt < nAlternative; kalt++){
          eppp( (ialt*nAlternative + jalt)*nAlternative + kalt ) += (prob(ialt, iy)*prob(jalt, iy)*prob(kalt, iy));
          epppu.row( (ialt*nAlternative + jalt)*nAlternative + kalt ) += ((prob(ialt, iy)*prob(jalt, iy)*prob(kalt, iy)) * 
            u.col(iy).array()).matrix();
        }
        
      } // jalt
    } // ialt
  }
  // Turn the sums into averages to get expectation terms
  ep = ep.array() / ny;
  epu = epu.array() / ny;
  epp = epp.array() / ny;
  eppu = eppu.array() / ny;
  eppp = eppp.array() / ny;
  epppu = epppu.array() / ny;
  
  Eigen::VectorXd U11temp = Eigen::VectorXd::Zero(b_mean.size());
  Eigen::MatrixXd U11     = Eigen::MatrixXd::Zero(nAlternative, b_mean.size());
  Eigen::VectorXd U12temp = Eigen::VectorXd::Zero(b_mean.size());
  Eigen::MatrixXd U12     = Eigen::MatrixXd::Zero(nAlternative, b_mean.size());
  Eigen::MatrixXd U21     = Eigen::MatrixXd::Zero( pow(nAlternative, 2), b_mean.size());
  Eigen::MatrixXd U22     = Eigen::MatrixXd::Zero( pow(nAlternative, 2), b_mean.size());
  Eigen::VectorXd U21temp = Eigen::VectorXd::Zero(b_mean.size());
  Eigen::VectorXd U22temp = Eigen::VectorXd::Zero(b_mean.size());
  
  int ics = 0;
  int ics2 = 0;
  int ind1 = 0;
  int ind2 = 0;
  
  // selected elements from ep vecs
  double eppp1 = 0.0;
  double eppp2 = 0.0;
  
  for (int ialt = 0; ialt < nAlternative; ialt++){
    
    // Zero out temp variables
    U11temp = Eigen::VectorXd::Zero(b_mean.size());
    U12temp = Eigen::VectorXd::Zero(b_mean.size());
    
    // Increment choice set counter
    if (ialt == (ics+1)*choiceSetSize){
      ics++;
    }
    
    for (int jalt = 0; jalt < choiceSetSize; jalt++){
      U11temp = (U11temp.transpose().array() + epp(ialt*nAlternative + ics*choiceSetSize + jalt) * X.row(ics*choiceSetSize + jalt).array()).transpose();
      for (int iparam = 0; iparam < b_mean.size(); iparam++){
        U12temp(iparam) = U12temp(iparam) + X(ics*choiceSetSize + jalt, iparam) * eppu(ialt*nAlternative + ics*choiceSetSize + jalt, iparam);
      } // iparam
    } // jalt
    
    
    U11.row(ialt) = ep(ialt) * X.row(ialt).array() - U11temp.transpose().array();
    
    for (int iparam = 0; iparam < b_mean.size(); iparam++){
      U12(ialt, iparam) = X(ialt, iparam) * epu(ialt, iparam) - U12temp(iparam);
    } // iparam
    
    ics2 = 0; // reset second choice set counter
    for (int ialt2 = 0; ialt2 < nAlternative; ialt2++){
      // Increment second choice set counter
      if (ialt2 == (ics2+1)*choiceSetSize){
        ics2++;
      }
      
      // Case where comparing alternatives in the same choice set
      if (ics == ics2){
        
        if (ialt == ialt2){
          U21.row(ialt*nAlternative + ialt2) = U11.row(ialt);
          U22.row(ialt*nAlternative + ialt2) = U12.row(ialt);
        } else {
          U21.row(ialt*nAlternative + ialt2) = U21.row(ialt*nAlternative + ialt2).array() * 0.0;
          U22.row(ialt*nAlternative + ialt2) = U22.row(ialt*nAlternative + ialt2).array() * 0.0;
        }
        // Otherwise the choice sets don't match
      } else {
        // Zero out temporary variables
        U21temp = Eigen::VectorXd::Zero(b_mean.size());
        U22temp = Eigen::VectorXd::Zero(b_mean.size());
        // Loop over all alternatives
        for (int jalt = 0; jalt < choiceSetSize; jalt++){
          ind1 = (ialt *nAlternative + ics * choiceSetSize + jalt) * nAlternative + ialt2;
          ind2 = (ialt2*nAlternative + ics2* choiceSetSize + jalt) * nAlternative + ialt;
          
          
          
          eppp1 = eppp( ind1 );
          eppp2 = eppp( ind2 );
          
          for (int iparam = 0; iparam < b_mean.size(); iparam++){
            U21temp(iparam) += eppp1 * X(ics*choiceSetSize + jalt, iparam) + 
              eppp2 * X(ics2*choiceSetSize + jalt, iparam);
          }
          
          
          
          // Loop over each parameter in the model
          for (int iparam = 0; iparam < b_mean.size(); iparam++){
            U22temp(iparam) = U22temp(iparam) + 
              X(ics*choiceSetSize + jalt, iparam) * 
              epppu(ind1, iparam) +
              X(ics2*choiceSetSize + jalt, iparam) * 
              epppu(ind2, iparam);
          } // iparam
          
        } // /jalt
        
        for (int iparam = 0; iparam < b_mean.size(); iparam++){
          U21(ialt*nAlternative + ialt2, iparam) = epp(ialt*nAlternative + ialt2, iparam) * 
            (X(ialt, iparam)+X(ialt2, iparam)) -
            U21temp(iparam);
          U22(ialt*nAlternative + ialt2, iparam) = (X(ialt, iparam) + X(ialt2, iparam)) * 
            eppu(ialt*nAlternative + ialt2, iparam) -
            U22temp(iparam);
        }
      } // else
    } // ialt2
  } // ialt
  
  Eigen::MatrixXd XX = Eigen::MatrixXd::Zero( pow(nAlternative, 2), b_mean.size());
  for (int iparam = 0; iparam < b_mean.size(); iparam++){
    for (int ialt = 0; ialt < nAlternative; ialt++){
      for (int jalt = 0; jalt < nAlternative; jalt++){
        XX(ialt*nAlternative + jalt, iparam) = X(ialt, iparam) * X(jalt, iparam);
      }
    }
  }
  
  Eigen::MatrixXd DU11 = X.transpose()  * U11.matrix();
  Eigen::MatrixXd DU12 = X.transpose()  * U12;
  Eigen::MatrixXd DU21 = XX.matrix().transpose() * U21.matrix();
  Eigen::MatrixXd DU22 = XX.matrix().transpose() * U22.matrix();
  
  Eigen::MatrixXd MME = Eigen::MatrixXd::Zero(2 * b_mean.size(), ny);
  MME.block(0, 0, b_mean.size(), ny) = X.transpose() * Yall;
  MME.block(b_mean.size(), 0, b_mean.size(), ny) = (X.transpose() * Yall).array().pow(2).matrix();
  
  Eigen::MatrixXd VM  = Eigen::MatrixXd::Zero(2 * b_mean.size(), 2 * b_mean.size());
  Eigen::MatrixXd VMi = Eigen::MatrixXd::Zero( b_mean.size() + nRandEffect, b_mean.size() + nRandEffect);
  Eigen::VectorXd EM  = Eigen::VectorXd::Zero(2 * b_mean.size());
  
  EM = MME.rowwise().sum().array() / ny;
  
  VM = ((MME * MME.transpose()).array() / ny).matrix() - EM * EM.transpose();
  VMi = VM.block(0, 0, b_mean.size() + nRandEffect, b_mean.size() + nRandEffect);
  
  Eigen::MatrixXd DU = Eigen::MatrixXd::Zero( nBeta + nRandEffect,
                                              nBeta + nRandEffect );
  
  DU.block(0, 0, nBeta, nBeta) = DU11;
  DU.block(0, nBeta, nBeta, nRandEffect) = DU12;
  DU.block(nBeta, 0, nRandEffect, nBeta) = DU21;
  DU.block(nBeta, nBeta, nRandEffect, nRandEffect) = DU22.block(0, 0, nRandEffect, nRandEffect);
  
  Eigen::MatrixXd stackedMatrix = Eigen::MatrixXd::Zero( 2*(b_mean.size() + nRandEffect), b_mean.size() + nRandEffect);
  stackedMatrix.block(0, 0, b_mean.size() + nRandEffect, b_mean.size() + nRandEffect) = DU;
  stackedMatrix.block(b_mean.size() + nRandEffect, 0, b_mean.size() + nRandEffect, b_mean.size() + nRandEffect) = VMi;
  return(stackedMatrix);
  
}
