// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "approxutils.h"
#include <limits>

// [[Rcpp::depends(RcppEigen)]]


// Evaluate the negative log integrand for the requested Laplace integral.
// integral_type: 0 = denominator, 1 = probability numerator, 2 = shifted square.
double calculate_laplace_objective(Eigen::VectorXd Y,
                                   Eigen::MatrixXd X,
                                   Eigen::VectorXd beta,
                                   Eigen::VectorXd u,
                                   Eigen::MatrixXd OmegaRE,
                                   Eigen::VectorXd SigmaRE,
                                   int n_choice_set,
                                   double addc,
                                   int integral_type,
                                   int integral_index,
                                   Eigen::VectorXd & prob){
  int choice_set_size = X.rows() / n_choice_set;
  Eigen::VectorXd XB = X * beta;
  Eigen::VectorXd XB_cs = Eigen::VectorXd::Zero(choice_set_size);
  Eigen::VectorXd log_prob = Eigen::VectorXd::Zero(X.rows());
  double objective = 0.5 * u.dot(OmegaRE * u);
  double max_XB = 0.0;
  double log_sum_expXB = 0.0;

  for (int ichoice = 0; ichoice < n_choice_set; ichoice++){
    XB_cs = XB.segment(ichoice*choice_set_size, choice_set_size);
    max_XB = XB_cs.maxCoeff();
    XB_cs = XB_cs.array() - max_XB;
    log_sum_expXB = log(XB_cs.array().exp().sum());
    log_prob.segment(ichoice*choice_set_size, choice_set_size) =
      XB_cs.array() - log_sum_expXB;
  }
  prob = log_prob.array().exp();
  objective -= Y.dot(log_prob);

  if (integral_type == 1){
    objective -= log_prob(integral_index);
  }
  if (integral_type == 2){
    objective -= log(pow(u(integral_index), 2) /
      pow(sqrt(SigmaRE(integral_index)), 3) + addc);
  }
  return objective;
}

// All integral ratios share the same constants, so only this log determinant
// and the negative log integrand are needed. Require a local maximum.
double calculate_laplace_log_determinant(Eigen::MatrixXd H,
                                         int iy,
                                         int integral_type,
                                         int integral_index){
  if (!H.allFinite()){
    Rcpp::stop("Laplace Hessian is not finite for draw %d, integral type %d, index %d.",
               iy + 1, integral_type, integral_index + 1);
  }
  Eigen::LLT<Eigen::MatrixXd> H_factor(-H);
  if (H_factor.info() != Eigen::Success){
    Rcpp::stop("Laplace mode has invalid curvature for draw %d, integral type %d, index %d.",
               iy + 1, integral_type, integral_index + 1);
  }
  Eigen::MatrixXd chol_H = H_factor.matrixU();
  return 2.0 * chol_H.diagonal().array().log().sum();
}

// Newton mode searches and Laplace information calculation.
Eigen::MatrixXd laplace_apply_newton_method(Eigen::MatrixXd Yall,
                                            Eigen::MatrixXd X,
                                            Eigen::MatrixXd u,
                                            Eigen::VectorXd b_mean,
                                            Eigen::MatrixXd Omega,
                                            Eigen::VectorXd Sigma_R,
                                            int n_random_effect,
                                            int n_choice_set)
{
  
  // Determine some dimensions/aspects of the design problem from the inputs
  int n_parameter         = b_mean.size();
  int n_fixed_effect      = n_parameter - n_random_effect;
  int n_alternative_total = X.rows();
  int choice_set_size     = n_alternative_total / n_choice_set;
  int n_samples           = Yall.cols();
  
  
  // Inverse of random effects variance and corresponding subblock (columns)
  // of the model matrix (XR)
  Eigen::MatrixXd OmegaRE = Omega.block(0, 0, n_random_effect, n_random_effect);
  Eigen::MatrixXd XR = X.block(0, 0, X.rows(), n_random_effect);
  
  // Allocate space for terms in maximization
  Eigen::VectorXd u0_mle     = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd mstart     = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd u1_mle     = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd Epy        = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd Eu2        = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd log_numerator  = Eigen::VectorXd::Zero(X.rows());
  Eigen::VectorXd log_numerator2 = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd prob       = Eigen::VectorXd::Zero(n_alternative_total);
  Eigen::VectorXd norm_prob  = Eigen::VectorXd::Zero(n_alternative_total);
  Eigen::MatrixXd delta      = Eigen::MatrixXd::Zero(n_alternative_total, n_alternative_total);
  Eigen::VectorXd gradient   = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::MatrixXd H          = Eigen::MatrixXd::Zero(n_random_effect, n_random_effect);
  Eigen::MatrixXd Hinv       = Eigen::MatrixXd::Zero(n_random_effect, n_random_effect);
  
  // Intermediate terms related to maximization
  Eigen::VectorXd SigmaRE          = Sigma_R.segment(0, n_random_effect);
  
  // Allocate space for a single sample (temporary variable)
  Eigen::VectorXd Y = Eigen::VectorXd::Zero(Yall.rows());
  
  // Rename the mean as the coefficient
  Eigen::VectorXd beta   = b_mean;
  
  // Model Matrix for a given choice set
  Eigen::MatrixXd X_cs_rdm  = Eigen::MatrixXd::Zero(choice_set_size, n_random_effect);
  
  // Variables related to inverse of variance-covariance
  Eigen::VectorXd score1        = Eigen::VectorXd::Zero(n_parameter);
  Eigen::VectorXd score2        = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::MatrixXd I11           = Eigen::MatrixXd::Zero(n_parameter, n_parameter);
  Eigen::MatrixXd I12           = Eigen::MatrixXd::Zero(n_parameter, n_random_effect);
  Eigen::MatrixXd I22           = Eigen::MatrixXd::Zero(n_random_effect, n_random_effect);
  
  // Terms related to algorithm
  bool converged    = false;
  double tol        = 0.000000000001; // shifted-square ratios need accurate modes
  double stepsize   = 0.0;
  int max_iter      = 200;
  double addc       = 10000.0;
  double mleterm    = 0.0;
  int i1            = 0;
  
  // stores p for a choice set, used in coverged = true loop
  Eigen::VectorXd pj = Eigen::VectorXd::Zero(choice_set_size);
  
  int max_backtrack = 50;
  bool accepted = false;
  double objective = 0.0;
  double objective_new = 0.0;
  double objective_tolerance = 0.0;
  double decrease = 0.0;
  double log_denom = 0.0;
  double log_H_det = 0.0;
  Eigen::VectorXd direction = Eigen::VectorXd::Zero(n_random_effect);
  Eigen::VectorXd beta_new = b_mean;
  Eigen::VectorXd prob_new = Eigen::VectorXd::Zero(n_alternative_total);

  // Integral types: 0 = denominator, 1 = probability, 2 = shifted square.
  for (int iy = 0; iy < n_samples; iy++){
    Y = Yall.col(iy);
    u0_mle = u.col(iy);

    // Denominator mode.
    converged = false;
    for (int iattempt = 0; iattempt < 2 && !converged; iattempt++){
      if (iattempt == 1){
        u0_mle.setZero();
      }
      for (int inewton_iter = 0; inewton_iter < max_iter; inewton_iter++){
        beta = b_mean;
        beta.head(n_random_effect) += u0_mle;
        objective = calculate_laplace_objective(Y, X, beta, u0_mle, OmegaRE,
          SigmaRE, n_choice_set, addc, 0, -1, prob);
        delta = calculate_delta(prob, n_choice_set);
        gradient = XR.transpose() * (Y - prob) - OmegaRE * u0_mle;
        H = -(XR.transpose() * delta * XR) - OmegaRE;
        if (!std::isfinite(objective) || !gradient.allFinite()){
          break;
        }
        if (gradient.cwiseAbs().maxCoeff() < tol){
          converged = true;
          break;
        }

        Hinv = H.inverse();
        direction = -(Hinv * gradient);
        decrease = gradient.dot(direction);
        if (!direction.allFinite() || !std::isfinite(decrease) || decrease <= 0.0){
          // A numerator Hessian need not be negative definite away from its mode.
          // Use a bounded gradient direction if Newton is not a descent direction.
          direction = gradient / (1.0 + gradient.norm());
          decrease = gradient.dot(direction);
        }
        if (!std::isfinite(decrease) || decrease <= 0.0){
          break;
        }

        stepsize = 1.0;
        accepted = false;
        objective_tolerance = 16.0 * std::numeric_limits<double>::epsilon() *
          (1.0 + fabs(objective));
        for (int ibacktrack = 0; ibacktrack < max_backtrack; ibacktrack++){
          u1_mle = u0_mle + stepsize * direction;
          beta_new = b_mean;
          beta_new.head(n_random_effect) += u1_mle;
          objective_new = calculate_laplace_objective(Y, X, beta_new, u1_mle,
            OmegaRE, SigmaRE, n_choice_set, addc, 0, -1, prob_new);
          if (std::isfinite(objective_new) &&
              objective_new <= objective - 0.0001 * stepsize * decrease + objective_tolerance){
            u0_mle = u1_mle;
            accepted = true;
            break;
          }
          stepsize *= 0.5;
        }
        if (!accepted){
          break;
        }
      }

      // Re-evaluate the objective, gradient, and Hessian at the accepted point.
      beta = b_mean;
      beta.head(n_random_effect) += u0_mle;
      objective = calculate_laplace_objective(Y, X, beta, u0_mle, OmegaRE,
        SigmaRE, n_choice_set, addc, 0, -1, prob);
      delta = calculate_delta(prob, n_choice_set);
      gradient = XR.transpose() * (Y - prob) - OmegaRE * u0_mle;
      H = -(XR.transpose() * delta * XR) - OmegaRE;
      converged = std::isfinite(objective) && gradient.allFinite() &&
        gradient.cwiseAbs().maxCoeff() < tol;
    }
    if (!converged){
      Rcpp::stop("Laplace mode solver failed for draw %d, integral type %d, index %d; covariance was not computed.",
                 iy + 1, 0, -1 + 1);
    }

    log_H_det = calculate_laplace_log_determinant(H, iy, 0, -1);
    log_denom = -objective - 0.5 * log_H_det;
    mstart = u0_mle;

    // Probability numerator modes.
    for (int j1 = 0; j1 < X.rows(); j1++){
      u0_mle = mstart;
      i1 = j1 / choice_set_size;
      X_cs_rdm = X.block(i1*choice_set_size, 0, choice_set_size, n_random_effect);
      converged = false;
      for (int iattempt = 0; iattempt < 2 && !converged; iattempt++){
        if (iattempt == 1){
          u0_mle.setZero();
        }
        for (int inewton_iter = 0; inewton_iter < max_iter; inewton_iter++){
          beta = b_mean;
          beta.head(n_random_effect) += u0_mle;
          objective = calculate_laplace_objective(Y, X, beta, u0_mle, OmegaRE,
            SigmaRE, n_choice_set, addc, 1, j1, prob);
          delta = calculate_delta(prob, n_choice_set);
          pj = prob.segment(i1*choice_set_size, choice_set_size);
          pj(j1 % choice_set_size) -= 1.0;
          pj *= -1.0;
          gradient = X_cs_rdm.transpose() * pj +
            XR.transpose() * (Y - prob) - OmegaRE * u0_mle;
          H = -(X_cs_rdm.transpose() *
            delta.block(i1*choice_set_size, i1*choice_set_size,
                        choice_set_size, choice_set_size) * X_cs_rdm) -
            XR.transpose() * delta * XR - OmegaRE;
          if (!std::isfinite(objective) || !gradient.allFinite()){
            break;
          }
          if (gradient.cwiseAbs().maxCoeff() < tol){
            converged = true;
            break;
          }

          Hinv = H.inverse();
          direction = -(Hinv * gradient);
          decrease = gradient.dot(direction);
          if (!direction.allFinite() || !std::isfinite(decrease) || decrease <= 0.0){
            // A numerator Hessian need not be negative definite away from its mode.
            // Use a bounded gradient direction if Newton is not a descent direction.
            direction = gradient / (1.0 + gradient.norm());
            decrease = gradient.dot(direction);
          }
          if (!std::isfinite(decrease) || decrease <= 0.0){
            break;
          }

          stepsize = 1.0;
          accepted = false;
          objective_tolerance = 16.0 * std::numeric_limits<double>::epsilon() *
            (1.0 + fabs(objective));
          for (int ibacktrack = 0; ibacktrack < max_backtrack; ibacktrack++){
            u1_mle = u0_mle + stepsize * direction;
            beta_new = b_mean;
            beta_new.head(n_random_effect) += u1_mle;
            objective_new = calculate_laplace_objective(Y, X, beta_new, u1_mle,
              OmegaRE, SigmaRE, n_choice_set, addc, 1, j1, prob_new);
            if (std::isfinite(objective_new) &&
                objective_new <= objective - 0.0001 * stepsize * decrease + objective_tolerance){
              u0_mle = u1_mle;
              accepted = true;
              break;
            }
            stepsize *= 0.5;
          }
          if (!accepted){
            break;
          }
        }

        // Re-evaluate the objective, gradient, and Hessian at the accepted point.
        beta = b_mean;
        beta.head(n_random_effect) += u0_mle;
        objective = calculate_laplace_objective(Y, X, beta, u0_mle, OmegaRE,
          SigmaRE, n_choice_set, addc, 1, j1, prob);
        delta = calculate_delta(prob, n_choice_set);
        pj = prob.segment(i1*choice_set_size, choice_set_size);
        pj(j1 % choice_set_size) -= 1.0;
        pj *= -1.0;
        gradient = X_cs_rdm.transpose() * pj +
          XR.transpose() * (Y - prob) - OmegaRE * u0_mle;
        H = -(X_cs_rdm.transpose() *
          delta.block(i1*choice_set_size, i1*choice_set_size,
                      choice_set_size, choice_set_size) * X_cs_rdm) -
          XR.transpose() * delta * XR - OmegaRE;
        converged = std::isfinite(objective) && gradient.allFinite() &&
          gradient.cwiseAbs().maxCoeff() < tol;
      }
      if (!converged){
        Rcpp::stop("Laplace mode solver failed for draw %d, integral type %d, index %d; covariance was not computed.",
                   iy + 1, 1, j1 + 1);
      }

      log_H_det = calculate_laplace_log_determinant(H, iy, 1, j1);
      log_numerator(j1) = -objective - 0.5 * log_H_det;
    }

    // Shifted-square numerator modes. Keep the original addc approximation.
    for (int irandeff = 0; irandeff < n_random_effect; irandeff++){
      u0_mle = mstart;
      converged = false;
      for (int iattempt = 0; iattempt < 2 && !converged; iattempt++){
        if (iattempt == 1){
          u0_mle.setZero();
        }
        for (int inewton_iter = 0; inewton_iter < max_iter; inewton_iter++){
          beta = b_mean;
          beta.head(n_random_effect) += u0_mle;
          objective = calculate_laplace_objective(Y, X, beta, u0_mle, OmegaRE,
            SigmaRE, n_choice_set, addc, 2, irandeff, prob);
          delta = calculate_delta(prob, n_choice_set);
          gradient = XR.transpose() * (Y - prob) - OmegaRE * u0_mle;
          mleterm = 2.0 * u0_mle(irandeff) /
            (pow(u0_mle(irandeff), 2) + addc * pow(sqrt(SigmaRE(irandeff)), 3));
          gradient(irandeff) += mleterm;
          mleterm = 2.0 * (addc * pow(sqrt(SigmaRE(irandeff)), 3) -
            pow(u0_mle(irandeff), 2)) /
            pow(pow(u0_mle(irandeff), 2) + addc * pow(sqrt(SigmaRE(irandeff)), 3), 2);
          H = -(XR.transpose() * delta * XR) - OmegaRE;
          H(irandeff, irandeff) += mleterm;
          if (!std::isfinite(objective) || !gradient.allFinite()){
            break;
          }
          if (gradient.cwiseAbs().maxCoeff() < tol){
            converged = true;
            break;
          }

          Hinv = H.inverse();
          direction = -(Hinv * gradient);
          decrease = gradient.dot(direction);
          if (!direction.allFinite() || !std::isfinite(decrease) || decrease <= 0.0){
            // A numerator Hessian need not be negative definite away from its mode.
            // Use a bounded gradient direction if Newton is not a descent direction.
            direction = gradient / (1.0 + gradient.norm());
            decrease = gradient.dot(direction);
          }
          if (!std::isfinite(decrease) || decrease <= 0.0){
            break;
          }

          stepsize = 1.0;
          accepted = false;
          objective_tolerance = 16.0 * std::numeric_limits<double>::epsilon() *
            (1.0 + fabs(objective));
          for (int ibacktrack = 0; ibacktrack < max_backtrack; ibacktrack++){
            u1_mle = u0_mle + stepsize * direction;
            beta_new = b_mean;
            beta_new.head(n_random_effect) += u1_mle;
            objective_new = calculate_laplace_objective(Y, X, beta_new, u1_mle,
              OmegaRE, SigmaRE, n_choice_set, addc, 2, irandeff, prob_new);
            if (std::isfinite(objective_new) &&
                objective_new <= objective - 0.0001 * stepsize * decrease + objective_tolerance){
              u0_mle = u1_mle;
              accepted = true;
              break;
            }
            stepsize *= 0.5;
          }
          if (!accepted){
            break;
          }
        }

        // Re-evaluate the objective, gradient, and Hessian at the accepted point.
        beta = b_mean;
        beta.head(n_random_effect) += u0_mle;
        objective = calculate_laplace_objective(Y, X, beta, u0_mle, OmegaRE,
          SigmaRE, n_choice_set, addc, 2, irandeff, prob);
        delta = calculate_delta(prob, n_choice_set);
        gradient = XR.transpose() * (Y - prob) - OmegaRE * u0_mle;
        mleterm = 2.0 * u0_mle(irandeff) /
          (pow(u0_mle(irandeff), 2) + addc * pow(sqrt(SigmaRE(irandeff)), 3));
        gradient(irandeff) += mleterm;
        mleterm = 2.0 * (addc * pow(sqrt(SigmaRE(irandeff)), 3) -
          pow(u0_mle(irandeff), 2)) /
          pow(pow(u0_mle(irandeff), 2) + addc * pow(sqrt(SigmaRE(irandeff)), 3), 2);
        H = -(XR.transpose() * delta * XR) - OmegaRE;
        H(irandeff, irandeff) += mleterm;
        converged = std::isfinite(objective) && gradient.allFinite() &&
          gradient.cwiseAbs().maxCoeff() < tol;
      }
      if (!converged){
        Rcpp::stop("Laplace mode solver failed for draw %d, integral type %d, index %d; covariance was not computed.",
                   iy + 1, 2, irandeff + 1);
      }

      log_H_det = calculate_laplace_log_determinant(H, iy, 2, irandeff);
      log_numerator2(irandeff) = -objective - 0.5 * log_H_det;
    }

    // Common Gaussian normalization constants cancel in these integral ratios.
    Epy = (log_numerator.array() - log_denom).exp();
    if (!Epy.allFinite()){
      Rcpp::stop("Laplace probability integral ratio is not finite for draw %d.", iy + 1);
    }
    norm_prob = normalize_probabilities(Epy, n_choice_set);
    for (int i = 0; i < n_random_effect; i++){
      Eu2(i) = exp(log_numerator2(i) - log_denom) - addc;
      if (Eu2(i) < 0.0){
        Eu2(i) *= -1.0;
      }
    }
    if (!norm_prob.allFinite() || !Eu2.allFinite()){
      Rcpp::stop("Laplace integral ratio is not finite for draw %d.", iy + 1);
    }

    score1 = X.transpose() * (Y - norm_prob);
    score2 = -(SigmaRE.array().rsqrt()) + Eu2.array();
    I11 += score1 * score1.transpose();
    I12 += score1 * score2.transpose();
    I22 += score2 * score2.transpose();
  }

  // Store the final inverse variance-covariance matrix
  Eigen::MatrixXd vcov_inverse = Eigen::MatrixXd::Zero(n_parameter + n_random_effect,
                                                       n_parameter + n_random_effect);
  
  vcov_inverse.block(0, 0, n_parameter, n_parameter)                             = I11;
  vcov_inverse.block(0, n_parameter, n_parameter, n_random_effect)               = I12;
  vcov_inverse.block(n_parameter, 0, n_random_effect, n_parameter)               = I12.transpose();
  vcov_inverse.block(n_parameter, n_parameter, n_random_effect, n_random_effect) = I22;
  
  // Normalize by number of replicates of Y
  for (int i = 0; i < n_parameter + n_random_effect; i++){
    for (int j = 0; j < n_parameter + n_random_effect; j++){
      vcov_inverse(i,j) *= 1.0 / n_samples;
    }
  }
  
  // Return inverse of variance-covariance based on laplace approximation
  return(vcov_inverse);
}




//' Laplace approximation to the information matrix
 //'
 //' @description Calculates the variance-covariance matrix for the model parameters 
 //' under the Laplace Approximation. This is not intended to be directly used,
 //' instead a user should use the \code{\link{varcov_approx_PML}} function with the "Laplace" option.
 //'
 //' @param n_samples Number of response samples to draw
 //' @param X The model matrix
 //' @param b_mean The vector of effect means
 //' @param var_vec The vector of effect variances
 //' @param n_choice_set The number of choice sets in the discrete choice experiment
 //' @return A matrix containing the variance covariance matrix of the model parameters based on the Laplace approximation
 //' @keywords internal
 // [[Rcpp::export]]
 Eigen::MatrixXd Laplace_approx_vcov(int n_samples,
                                     Eigen::MatrixXd X,
                                     Eigen::VectorXd b_mean,
                                     Eigen::VectorXd var_vec,
                                     int n_choice_set){
   
   if (n_samples <= 0 || n_choice_set <= 0 || X.rows() == 0 ||
       X.rows() % n_choice_set != 0 || X.rows() / n_choice_set < 2 ||
       X.cols() != b_mean.size() || var_vec.size() != b_mean.size() ||
       !X.allFinite() || !b_mean.allFinite() || !var_vec.allFinite() ||
       (var_vec.array() < 0).any() || !(var_vec.array() > 0).any()){
     Rcpp::stop("Laplace requires consistent finite inputs, positive sample/choice-set counts, and at least one random effect.");
   }

   // Determine some aspects of the design problem from the inputs
   int n_alternative_total  = X.rows(); // total number of alternatives
   int n_beta               = b_mean.size();
   int choice_set_size      = n_alternative_total / n_choice_set;
   
   // Track how many effects are random
   int n_random_effect = 0;
   for (int i = 0; i < var_vec.size(); i++){
     if (var_vec(i) > 0){
       n_random_effect++;
     }
   }
   
   // How many effects are fixed
   int n_fixed_effect = var_vec.size() - n_random_effect;
   if (!(var_vec.head(n_random_effect).array() > 0).all()){
     Rcpp::stop("Laplace requires random-effect columns first; use varcov_approx_PML().");
   }

   
   // Allocate space for final precision about the model parameters
   Eigen::MatrixXd vcov_inverse = Eigen::MatrixXd::Zero(n_beta + n_random_effect,
                                                        n_beta + n_random_effect);
   
   
   // Inverse of prior covariance
   Eigen::MatrixXd Omega = Eigen::MatrixXd::Zero(n_beta, n_beta);
   Omega.block(0, 0, n_random_effect, n_random_effect) =
     var_vec.head(n_random_effect).array().inverse().matrix().asDiagonal();
   
   Eigen::VectorXd Sigma_R = var_vec;
   
   // Sample betas for each draw of Y
   Eigen::MatrixXd beta   = Eigen::MatrixXd::Zero(b_mean.size(), n_samples); // each column is a sample of beta
   Eigen::MatrixXd u      = Eigen::MatrixXd::Zero(n_random_effect, n_samples);
   Eigen::VectorXd sd_vec = Eigen::VectorXd::Zero(n_random_effect);
   
   // Draw samples of u
   for (int ip = 0; ip < n_random_effect; ip++){
     sd_vec(ip) = sqrt(var_vec(ip));
     for (int iu = 0; iu < n_samples; iu++){
       u(ip, iu) = R::rnorm(0.0, sd_vec(ip));
     }
   }
   Eigen::MatrixXd sd_mat = sd_vec.asDiagonal();
   
   // Get the corresponding samples of beta (b + u)
   for (int iy = 0; iy < n_samples; iy++){
     beta.block(0, iy, n_random_effect, 1)              = b_mean.segment(0, n_random_effect) + u.col(iy);
     beta.block(n_random_effect, iy, n_fixed_effect, 1) = b_mean.segment(n_random_effect, n_fixed_effect);
   }
   
   // Calculate the epsilon terms (get added to linear predictor later)
   Eigen::MatrixXd epsilon = Eigen::MatrixXd::Zero(X.rows(), n_samples);
   for (int iy = 0; iy < n_samples; iy++){
     epsilon.col(iy) = sample_epsilon(X.rows());
   }
   
   // Allocate intermediate quantities for drawing samples of the outcome vector
   Eigen::VectorXd XB    = Eigen::VectorXd::Zero(X.rows());
   Eigen::VectorXd prob  = Eigen::VectorXd::Zero(n_alternative_total);
   Eigen::MatrixXd Y     = Eigen::MatrixXd::Zero(X.rows(), n_samples);
   int yind              = 0; // for tracking which y is selected
   
   // Loop over n_samples
   for (int iy = 0; iy < n_samples; iy++){
     
     // Evaluate linear predictor
     XB    = X * beta.col(iy) + epsilon.col(iy);
     
     // As a starting value, pick the most likely Y for each choice set
     for (int ichoice = 0; ichoice < n_choice_set; ichoice++){
       yind=0;
       // Find the selected value for this choice set
       for (int ialt = 1; ialt < choice_set_size; ialt++){
         if (XB(ichoice*choice_set_size + ialt) > XB(ichoice*choice_set_size + yind)) {
           yind = ialt;
         }
       }
       Y(ichoice*choice_set_size + yind, iy) = 1;
     }
   } // /iy
   
   // Nested optimizations to approximate inverse variance covariance
   vcov_inverse = laplace_apply_newton_method(Y,
                                              X,
                                              u,
                                              b_mean,
                                              Omega,
                                              Sigma_R,
                                              n_random_effect,
                                              n_choice_set);
   
   
   // Invert to get variance-covariance matrix
   vcov_inverse = vcov_inverse.inverse();
   
   return vcov_inverse;
 }
