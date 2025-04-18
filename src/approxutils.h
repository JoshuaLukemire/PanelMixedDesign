#ifndef __PMLH__
#define __PMLH__

#include <RcppEigen.h>
#include <random>

Eigen::MatrixXd sampleRandomBeta(int n_u,
                                 Eigen::VectorXd b_mean,
                                 Eigen::VectorXd var_vec);

Eigen::VectorXd calculate_response_probabilities(Eigen::MatrixXd X,
                                       Eigen::VectorXd beta,
                                       int n_choice);

Eigen::VectorXd normalize_probabilities(Eigen::VectorXd probs,
                                       int n_choice);

Eigen::VectorXd sample_epsilon(int n_alternative_total);

Eigen::MatrixXd calculate_delta(Eigen::VectorXd p, int n_choice_set);


#endif // __PMLH__