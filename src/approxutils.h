#ifndef __PMLH__
#define __PMLH__

#include <RcppEigen.h>
#include <random>

Eigen::MatrixXd sampleRandomBeta(int n_u,
                                 Eigen::VectorXd b_mean,
                                 Eigen::VectorXd var_vec);

Eigen::VectorXd calculateResponseProbs(Eigen::MatrixXd X,
                                       Eigen::VectorXd beta,
                                       int n_choice);

Eigen::VectorXd normalizeProbabilities(Eigen::VectorXd probs,
                                       int n_choice);

Eigen::VectorXd sampleEpsilon(int totalAlternatives);

Eigen::MatrixXd calcDelta(Eigen::VectorXd p, int n_choice_set);


#endif // __PMLH__