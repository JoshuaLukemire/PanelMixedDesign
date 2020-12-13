// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "approxutils.h"


// [[Rcpp::depends(RcppEigen)]]

Eigen::VectorXd impsamp_calc_marginal_prob(Eigen::VectorXd Y, Eigen::MatrixXd probs,
                                           int nChoiceSet){
  // Find the number of samples
  int nU = probs.cols();
  
  // Initialize a vector for the marginal probabilities
  Eigen::VectorXd margProbs = Eigen::VectorXd::Zero(nU);
  
  // see https://stackoverflow.com/questions/44268968/eigen-vector-or-matrix-componentwise-to-power
  Eigen::MatrixXd powerTerms = probs.array() // go to array domain
                                    .pow(                 // element-wise power
  Y.replicate(1, probs.cols()).array() // repeat exponents to match size of A
                                    );
  
  // Calc probs for each of the random effect samples
  for (int iu = 0; iu < nU; iu++){
    margProbs(iu) = powerTerms.col(iu).prod();
  }
  
  return margProbs;
  
}



Eigen::MatrixXd impsamp_calc_info_given_Y(Eigen::VectorXd Y, Eigen::MatrixXd X, Eigen::VectorXd b_mean,
                                          Eigen::VectorXd var_vec, int nChoiceSet, int nU){
  
  int nRandEffect = 0;
  for (int i = 0; i < var_vec.size(); i++){
    if (var_vec(i) > 0){
      nRandEffect++;
    }
  }
  int nFixedEffect = var_vec.size() - nRandEffect;
  
  // Number of beta terms
  int effect_dim = b_mean.size();
  
  // Draw random effect samples and calculate corresponding effects (beta)
  // each row is an effect and each column is a sample
  std::random_device rd; 
  std::mt19937 gen(rd()); 
  
  //std::default_random_engine generator(); // todo random gen
  std::normal_distribution<double> distribution(0.0, 1.0);
  Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nRandEffect, nU);
  Eigen::MatrixXd beta_i = Eigen::MatrixXd::Zero(effect_dim, nU);
  Eigen::VectorXd sd_vec = Eigen::VectorXd::Zero(nRandEffect);
  for (int ip = 0; ip < nRandEffect; ip++){
    sd_vec(ip) = sqrt(var_vec(ip));
    for (int iu = 0; iu < nU; iu++){
      u(ip, iu) = distribution(gen);
    }
  }
  Eigen::MatrixXd sd_mat = sd_vec.asDiagonal();
  
  for (int iu = 0; iu < nU; iu++){
    beta_i.block(0,           iu,  nRandEffect, 1) = b_mean.segment(0,           nRandEffect) + sd_mat *  u.col(iu);
    beta_i.block(nRandEffect, iu, nFixedEffect, 1) = b_mean.segment(nRandEffect, nFixedEffect);
  }
  
  // Calculate the response probabilities
  Eigen::MatrixXd responseProbs = Eigen::MatrixXd::Zero(Y.size(), nU);
  for (int iu = 0; iu < nU; iu++){
    responseProbs.col(iu) = calculateResponseProbs(X, beta_i.col(iu), nChoiceSet);
  }

  // Calculate the marginal probabilities
  Eigen::VectorXd margProbs = impsamp_calc_marginal_prob(Y, responseProbs, nChoiceSet);
  
  // Average of the probabilities from above over nU is the marginal prob
  double marg_prob = margProbs.sum() / nU;
  
  Eigen::VectorXd f1 = (responseProbs.matrix() * margProbs.matrix()).array() * (1.0 / nU) * (1.0 / marg_prob);
  
  // score function for the regression coefficients
  Eigen::VectorXd score11 = X.matrix().transpose() * (Y.array() - f1.array()).matrix();
  
  // Upper P x P block of the information matrix
  Eigen::MatrixXd I11 = (score11 * score11.transpose()).array() * marg_prob;  
  
  // U22 (for the lower 2 x 2 block)
  Eigen::MatrixXd U22 = Eigen::MatrixXd::Zero(nRandEffect, nU);
  for (int iU = 0; iU < nU; iU++){
    for (int iP = 0; iP < nRandEffect; iP++){
      U22(iP, iU) = pow(sd_vec(iP) * u(iP, iU), 2) / pow(sqrt(var_vec(iP)), 3);
    }
  }
  
  Eigen::VectorXd f2 = (U22.matrix() * margProbs.matrix()).array() * (1.0 / nU) * (1.0 / marg_prob);
  Eigen::VectorXd score2 = sd_vec.array().inverse() * (-1.0) + f2.array();
  
  Eigen::MatrixXd I12 = (score11.matrix() * score2.matrix().transpose()).array() * marg_prob;
  Eigen::MatrixXd I22 = (score2.matrix()  * score2.matrix().transpose()).array() * marg_prob;
  
  Eigen::MatrixXd fisherInfo = Eigen::MatrixXd::Zero( effect_dim + nRandEffect,
                                                      effect_dim + nRandEffect );
  
  fisherInfo.block(0, 0, effect_dim, effect_dim) = I11;
  fisherInfo.block(0, effect_dim, effect_dim, nRandEffect) = I12;
  fisherInfo.block(effect_dim, 0, nRandEffect, effect_dim) = I12.transpose();
  fisherInfo.block(effect_dim, effect_dim, nRandEffect, nRandEffect) = I22.block(0, 0, nRandEffect, nRandEffect);
  
  return fisherInfo;
}

//' Approximate the information matrix using importance sampling
//'
//' @description Approximates the information matrix for the model parameters 
//' using importance sampling. This is not intended to be directly used,
//' instead a user should use the \code{\link{PMLInfoApprox}} function with 
//' the "importance" option.
//'
//' @param Y an enumeration of all possible responses. Can be obtained using
//'  the \code{\link{gen_all_choice_seq}} function
//' @param X The model matrix
//' @param b_mean The vector of effect means
//' @param var_vec The vector of effect variances
//' @param nChoiceSet The number of choice sets in the discrete choice experiment
//' @param nU Number of samples to use to evalutate the information matrix
//' @return A matrix containing the MSM evaluated of the variance
//' covariance matrix of the model parameters
// [[Rcpp::export]]
Eigen::MatrixXd importanceSampleFixedY(Eigen::MatrixXd Y,
                                       Eigen::MatrixXd X,
                                       Eigen::VectorXd b_mean,
                                       Eigen::VectorXd var_vec, int nChoiceSet, int nU){
  
  int nRandEffect = 0;
  for (int i = 0; i < var_vec.size(); i++){
    if (var_vec(i) > 0){
      nRandEffect++;
    }
  }
  
  Eigen::MatrixXd FI = Eigen::MatrixXd::Zero(b_mean.size() + nRandEffect, b_mean.size() + nRandEffect);
  
  int nY = Y.cols();
  
  for (int iY = 0; iY < nY; iY++){
    FI += impsamp_calc_info_given_Y(Y.col(iY), X, b_mean, var_vec, nChoiceSet, nU);
  }
  
  return FI;
  
}
