#include "approxutils.h"

//
// sampleRandomBeta
//
Eigen::MatrixXd sampleRandomBeta(int n_u,
                                 Eigen::VectorXd b_mean,
                                 Eigen::VectorXd var_vec){
  
  // Number of beta terms
  int effect_dim = b_mean.size();
  
  // Draw random effect samples and calculate corresponding effects (beta)
  // each row is an effect and each column is a sample
  GetRNGstate();
  Eigen::MatrixXd u = Eigen::MatrixXd::Zero(effect_dim, n_u);
  for (int ip = 0; ip < effect_dim; ip++){
    for (int iu = 0; iu < n_u; iu++){
      u(ip, iu) = R::rnorm(0.0, 1.0);
    }
  }
  PutRNGstate();
  
  // Generate the observed betas, note that multiplying by sqrt var_vec (0) will fix
  // things up for any fixed effects (u will be 0)
  Eigen::VectorXd sd_vec = var_vec.array().sqrt();
  Eigen::MatrixXd sd_mat = sd_vec.asDiagonal();
  Eigen::MatrixXd beta_i = sd_mat *  u;
  beta_i.colwise() += b_mean;
  
  return beta_i;
}

//
// calculateResponseProbabilitiesRaw
//
Eigen::VectorXd calculateResponseProbs(Eigen::MatrixXd X,
                                       Eigen::VectorXd beta,
                                       int n_choice){
  
  // Calculate the linear predictor
  Eigen::VectorXd XB = X * beta;
  Eigen::VectorXd expXB = XB.array().exp();
  
  // Total number of alternatives participant sees
  int n_alt_total = XB.rows();
  int choice_set_size = n_alt_total / n_choice;
  
  // Initialize probability vector
  Eigen::VectorXd prob = Eigen::VectorXd::Zero(n_alt_total);
  
  double prob_sum = 0.0;
  for (int ichoice = 0; ichoice < n_choice; ichoice++){
    
    // Calculate the total sum of choice probabilities for this set
    prob_sum = 0.0;
    for (int ialt = 0; ialt < choice_set_size; ialt++){
      prob_sum += expXB( ichoice*choice_set_size + ialt );
    }
    
    // Now normalize
    for (int ialt = 0; ialt < choice_set_size; ialt++){
      prob(ichoice*choice_set_size + ialt) = expXB( ichoice*choice_set_size + ialt ) / prob_sum;
    }
    
  }
  
  return prob;
}

Eigen::VectorXd normalizeProbabilities(Eigen::VectorXd probs,
                                       int n_choice){
  int n_alt_total = probs.rows();
  int choice_set_size = n_alt_total / n_choice;
  Eigen::VectorXd norm_probs = Eigen::VectorXd::Zero(probs.size());
  
  double prob_sum = 0.0;
  for (int ichoice = 0; ichoice < n_choice; ichoice++){
    prob_sum = 0.0;
    for (int ialt = 0; ialt < choice_set_size; ialt++){
      prob_sum += probs( ichoice*choice_set_size + ialt );
    }
    // Now normalize
    for (int ialt = 0; ialt < choice_set_size; ialt++){
      norm_probs(ichoice*choice_set_size + ialt) = probs( ichoice*choice_set_size + ialt ) / prob_sum;
    }
  }
  return(norm_probs);
}



Eigen::VectorXd sampleEpsilon(int totalAlternatives){
  Eigen::VectorXd eps = Eigen::VectorXd::Zero(totalAlternatives);
  //std::cout << "adjusting epsilon" << std::endl;
  //eps = (eps.array() + 1.0) / 2.0; 
  GetRNGstate();
  for (int ii = 0; ii < totalAlternatives; ii++){
    eps(ii) = R::runif(0.0, 1.0);
  }
  PutRNGstate();
  eps = (eps.array().log() * (-1.0) ).log() * (-1.0);
  return eps;
}


//
// Function to calculate Delta = diag(p) - pp'
//
Eigen::MatrixXd calcDelta(Eigen::VectorXd p,
                          int n_choice_set){
  
  int n_alt_total = p.size();
  int choice_set_size = n_alt_total / n_choice_set;
  
  Eigen::MatrixXd delta = Eigen::MatrixXd::Zero(n_alt_total, n_alt_total);
  
  for (int ichoice = 0; ichoice < n_choice_set; ichoice++){
    
    // Diagonal portion of delta
    for (int ialt = 0; ialt < choice_set_size; ialt++){
      delta(ichoice*choice_set_size + ialt, ichoice*choice_set_size + ialt) = p(ichoice*choice_set_size + ialt);
    }
    
    // Subtract off pp'
    //delta(Eigen::seq(ichoice*choice_set_size, ichoice*choice_set_size+choice_set_size-1),xxx
    //      Eigen::seq(ichoice*choice_set_size, ichoice*choice_set_size+choice_set_size-1)) -= xxx
    //        p(Eigen::seq(ichoice*choice_set_size, ichoice*choice_set_size+choice_set_size-1)) *xxx
    //        p(Eigen::seq(ichoice*choice_set_size, ichoice*choice_set_size+choice_set_size-1));xxx
    
    delta.block(ichoice*choice_set_size, ichoice*choice_set_size, choice_set_size, choice_set_size) -=
      (p.segment(ichoice*choice_set_size, choice_set_size) *
      p.segment(ichoice*choice_set_size, choice_set_size).transpose());
    
  } // /loop over choice sets
  
  return delta;
  
}

