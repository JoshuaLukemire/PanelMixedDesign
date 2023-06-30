// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "approxutils.h"

// each column is a response vector
void impsamp_calc_conditional_prob_matrix_inplace(Eigen::MatrixXd & condProbs,
                                           const Eigen::MatrixXd & Y,
                                           const Eigen::MatrixXd & probs,
                                           int nChoiceSet){
  // Find the number of samples
  int nU = probs.cols();
  int nY = Y.cols();
  
  // Initialize a matrix for the conditional probabilities
  //Eigen::MatrixXd condProbs = Eigen::MatrixXd::Zero(nU, nY);
  
  for (int iY = 0; iY < nY; iY++){
    
    // see https://stackoverflow.com/questions/44268968/eigen-vector-or-matrix-componentwise-to-power
    Eigen::MatrixXd powerTerms = probs.array() // go to array domain
                                      .pow(                 // element-wise power
    Y.col(iY).replicate(1, probs.cols()).array() // repeat exponents to match size of A
                                      );
    
    // Calc probs for each of the random effect samples
    for (int iu = 0; iu < nU; iu++){
      condProbs(iu, iY) = powerTerms.col(iu).prod();
    }
    
  }
  
}

// Need this version if adding Y terms "one at a time"
void impsamp_calc_conditional_prob_vector_inplace(Eigen::VectorXd & condProbs,
                                                  const Eigen::VectorXd & Y,
                                                  const Eigen::MatrixXd & probs,
                                                  int nChoiceSet){
  // Find the number of samples
  int nU = probs.cols();
  
    // see https://stackoverflow.com/questions/44268968/eigen-vector-or-matrix-componentwise-to-power
    Eigen::MatrixXd powerTerms = probs.array() // go to array domain
                                      .pow(                 // element-wise power
    Y.replicate(1, probs.cols()).array() // repeat exponents to match size of A
                                      );
    
    // Calc probs for each of the random effect samples
    for (int iu = 0; iu < nU; iu++){
      condProbs(iu) = powerTerms.col(iu).prod();
    }
  
}

// [[Rcpp::export]]
void impsamp_calc_conditional_prob_vector_inplace_R(Eigen::Map<Eigen::VectorXd> & condProbs,
                                                  const Eigen::VectorXd & Y,
                                                  const Eigen::MatrixXd & probs,
                                                  int nChoiceSet){
  // Find the number of samples
  int nU = probs.cols();
  
  // see https://stackoverflow.com/questions/44268968/eigen-vector-or-matrix-componentwise-to-power
  Eigen::MatrixXd powerTerms = probs.array() // go to array domain
                                    .pow(                 // element-wise power
  Y.replicate(1, probs.cols()).array() // repeat exponents to match size of A
                                    );
  
  // Calc probs for each of the random effect samples
  for (int iu = 0; iu < nU; iu++){
    condProbs(iu) = powerTerms.col(iu).prod();
  }
  
}



void impsamp_update_delta_inplace(Eigen::MatrixXd & delta,
                                  const Eigen::VectorXd & p,
                                  int n_choice_set){
  
  int n_alt_total = p.size();
  int choice_set_size = n_alt_total / n_choice_set;
  
  //Eigen::MatrixXd delta = Eigen::MatrixXd::Zero(n_alt_total, n_alt_total);
  
  delta *= 0.0;
  
  for (int ichoice = 0; ichoice < n_choice_set; ichoice++){
    
    // Diagonal portion of delta
    for (int ialt = 0; ialt < choice_set_size; ialt++){
      delta(ichoice*choice_set_size + ialt, ichoice*choice_set_size + ialt) = p(ichoice*choice_set_size + ialt);
    }
    
    delta.block(ichoice*choice_set_size, ichoice*choice_set_size, choice_set_size, choice_set_size) -=
      (p.segment(ichoice*choice_set_size, choice_set_size) *
      p.segment(ichoice*choice_set_size, choice_set_size).transpose());
    
  } // /loop over choice sets

}


// [[Rcpp::export]]
void impsamp_update_easy_terms_inplace(Eigen::Map<Eigen::MatrixXd> & Delta,
                                                  Eigen::Map<Eigen::MatrixXd> & PPt,
                                                  Eigen::Map<Eigen::MatrixXd> & PRatiosT,
                                                  Eigen::Map<Eigen::MatrixXd> & responseProbs,
                                                  Eigen::Map<Eigen::MatrixXd> & U22,
                                                  const Eigen::MatrixXd & X,
                                                  const Eigen::VectorXd & b_mean,
                                                  const Eigen::VectorXd & var_vec,
                                                  int nChoiceSet,
                                                  int nU){
  
  int nRandEffect = 0;
  for (int i = 0; i < var_vec.size(); i++){
    if (var_vec(i) > 0){
      nRandEffect++;
    }
  }
  int nFixedEffect = var_vec.size() - nRandEffect;
  
  // Number of beta terms
  int effect_dim = b_mean.size();
  
  // Temorary variables
  Eigen::MatrixXd delta_iu = Eigen::MatrixXd::Zero(X.rows(), X.rows());
  
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
  //Eigen::MatrixXd responseProbs = Eigen::MatrixXd::Zero(X.rows(), nU);
  for (int iu = 0; iu < nU; iu++){
    responseProbs.col(iu) = calculateResponseProbs(X, beta_i.col(iu), nChoiceSet);
    
    // Delta Term
    impsamp_update_delta_inplace(delta_iu, responseProbs.col(iu), nChoiceSet);
    Delta += delta_iu;
    
    // PPt term
    PPt += responseProbs.col(iu) * responseProbs.col(iu).transpose();
    
    // UUt term
    for (int iP = 0; iP < nRandEffect; iP++){
      U22(iP, iu) = pow(sd_vec(iP) * u(iP, iu), 2) / pow(sqrt(var_vec(iP)), 3);
    }
    PRatiosT += responseProbs.col(iu) * U22.col(iu).transpose();
    
  }
  
  Delta    *= (1.0 / nU);
  PPt      *= (1.0 / nU);
  PRatiosT *= (1.0 / nU);
  
}
