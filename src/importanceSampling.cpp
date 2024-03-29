// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "approxutils.h"


Eigen::MatrixXd impsamp_SampleY(Eigen::VectorXd expXB,
                               Eigen::VectorXd X,
                               int nChoiceSet){
  
  int choiceSetSize = X.rows() / nChoiceSet;
  Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(X.rows(), 1);
  
  int ind = 0;
  double maxval = 0.0;
  
  double total_prob = 0.0;
  double cumu_prob = 0.0;
  double rdraw = 0.0;

    // Determine the most likely outcome for this sampled Y
    for (int ic = 0; ic < nChoiceSet; ic++){
      
      total_prob = 0.0;
      for (int ialt = 0; ialt < choiceSetSize; ialt++){
          total_prob += expXB( ic*choiceSetSize + ialt );
      }

      ind = 0;
      double cumu_prob = 0.0;
      
      GetRNGstate();
      rdraw = R::runif(0.0, 1.0);
      PutRNGstate();
      
      for (int ialt = 0; ialt < choiceSetSize; ialt++){
        cumu_prob += expXB( ic*choiceSetSize + ialt ) / total_prob;
        if (cumu_prob > rdraw){
          ind = ic*choiceSetSize + ialt;
          rdraw = 999.0;
        }
      }
      // Set the most likely Y for this choice set to 1
      Y(ind, 0) = 1;
    }

  return Y;
}






Eigen::VectorXd impsamp_SampleEpsilon(int totalAlternatives){
  Eigen::VectorXd eps = Eigen::VectorXd::Zero(totalAlternatives);
  GetRNGstate();
  for (int ii = 0; ii < totalAlternatives; ii++){
    eps(ii) = R::runif(0.0, 1.0);
  }
  PutRNGstate();
  eps = (eps.array().log() * (-1.0) ).log() * (-1.0);
  return eps;
}



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


// [[Rcpp::export]]
Eigen::MatrixXd impsamp_calc_info_given_Y(Eigen::VectorXd Y,
                                          Eigen::MatrixXd X,
                                          Eigen::VectorXd b_mean,
                                          Eigen::VectorXd var_vec,
                                          int nChoiceSet,
                                          int nU,
                                          bool return_unweighted,
                                          Eigen::VectorXd & marginal_prob_store){
  
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

  // Calculate the marginal probabilities, each column corresponds to a u_i
  Eigen::VectorXd margProbs = impsamp_calc_marginal_prob(Y, responseProbs, nChoiceSet);
  
  // Average of the probabilities from above over nU is the marginal prob
  double marg_prob = margProbs.sum() / nU; // P(Y = y)
  
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
  //Eigen::VectorXd f2 = (U22.matrix() * margProbs.matrix()).array() * (1.0 / nU);
  Eigen::VectorXd score2 = sd_vec.array().inverse() * (-1.0) + f2.array();
  
  Eigen::MatrixXd I12 = (score11.matrix() * score2.matrix().transpose()).array() * marg_prob;
  Eigen::MatrixXd I22 = (score2.matrix()  * score2.matrix().transpose()).array() * marg_prob;
  
  Eigen::MatrixXd fisherInfo = Eigen::MatrixXd::Zero( effect_dim + nRandEffect,
                                                      effect_dim + nRandEffect );
  
  fisherInfo.block(0, 0, effect_dim, effect_dim) = I11;
  fisherInfo.block(0, effect_dim, effect_dim, nRandEffect) = I12;
  fisherInfo.block(effect_dim, 0, nRandEffect, effect_dim) = I12.transpose();
  fisherInfo.block(effect_dim, effect_dim, nRandEffect, nRandEffect) = I22.block(0, 0, nRandEffect, nRandEffect);
  
  marginal_prob_store(0) = marg_prob;
  if (return_unweighted == true){
    fisherInfo = fisherInfo / marg_prob;
  }
  
  return fisherInfo;
}



// [[Rcpp::export]]
Eigen::MatrixXd R_impsamp_calc_info_given_Y(Eigen::VectorXd Y,
                                          Eigen::MatrixXd X,
                                          Eigen::VectorXd b_mean,
                                          Eigen::VectorXd var_vec,
                                          int nChoiceSet,
                                          int nU,
                                          bool return_unweighted,
                                          Eigen::Map<Eigen::VectorXd> & marginal_prob_store){
  
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
  
  // Calculate the marginal probabilities, each column corresponds to a u_i
  Eigen::VectorXd margProbs = impsamp_calc_marginal_prob(Y, responseProbs, nChoiceSet);
  
  // Average of the probabilities from above over nU is the marginal prob
  double marg_prob = margProbs.sum() / nU; // P(Y = y)
  
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
  
  marginal_prob_store(0) = marg_prob;
  if (return_unweighted == true){
    fisherInfo = fisherInfo / marg_prob;
  }
  
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
  Eigen::VectorXd marginal_prob_store = Eigen::VectorXd::Zero(1);
  
  int nY = Y.cols();
  
  for (int iY = 0; iY < nY; iY++){
    FI += impsamp_calc_info_given_Y(Y.col(iY), X, b_mean, var_vec, nChoiceSet, nU, false, marginal_prob_store);
  }
  
  return FI;
  
}




// [[Rcpp::export]]
Eigen::MatrixXd importanceSample(     Eigen::MatrixXd X,
                                       Eigen::VectorXd b_mean,
                                       Eigen::VectorXd var_vec, int nChoiceSet, int nU, int nY){
  
  int nRandEffect = 0;
  for (int i = 0; i < var_vec.size(); i++){
    if (var_vec(i) > 0){
      nRandEffect++;
    }
  }
  
  Eigen::MatrixXd FI = Eigen::MatrixXd::Zero(b_mean.size() + nRandEffect, b_mean.size() + nRandEffect);
  Eigen::MatrixXd Y  = Eigen::MatrixXd::Zero(X.rows(), 1);
  Eigen::MatrixXd epsilon = Eigen::MatrixXd::Zero(X.rows(), 1);
  //Eigen::VectorXd XB    = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd expXB = Eigen::VectorXd::Zero(X.rows());
  
  Eigen::VectorXd XB    = Eigen::VectorXd::Ones(X.rows());
  Eigen::VectorXd Py    = Eigen::VectorXd::Zero(X.rows());
  
  Eigen::VectorXd marginal_prob_store = Eigen::VectorXd::Zero(1);
  
  double total_probability = 0.0;

  for (int iY = 0; iY < nY; iY++){
    epsilon.col(0) = impsamp_SampleEpsilon(X.rows());
    //XB    = X * (b_mean + epsilon.col(0));
    XB    = X * b_mean + epsilon.col(0);
    expXB = XB.array().exp();
    
    Y  = impsamp_SampleY(expXB, X, nChoiceSet);
    Py = normalizeProbabilities(expXB, nChoiceSet);
    //FI += impsamp_calc_info_given_Y(Y.col(0), X, b_mean, var_vec, nChoiceSet, nU, false, marginal_prob_store) /
    //  Py.array().pow(Y.col(0).array()).prod() / nY;
    FI += impsamp_calc_info_given_Y(Y.col(0), X, b_mean, var_vec, nChoiceSet, nU, false, marginal_prob_store);
    
    total_probability += marginal_prob_store(0);
  }
  
  FI = FI / total_probability;
  
  return FI;
  
}