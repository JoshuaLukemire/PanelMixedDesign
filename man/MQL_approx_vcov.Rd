% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MQL.R
\name{MQL_approx_vcov}
\alias{MQL_approx_vcov}
\title{Marginal Quasi-Likelihood approach to approximate the variance covariance matrix for the
model parameters}
\usage{
MQL_approx_vcov(model_matrix, effect_mean, effect_var, n_choice_set)
}
\arguments{
\item{model_matrix}{The model matrix (effects coded), should be (n_alternative * n_choice_set) x n_beta}

\item{effect_mean}{Vector of means for the effects coded attribute effects}

\item{effect_var}{Vector of variances for the effects coded attribute effects}

\item{n_choice_set}{Number of choice sets}
}
\value{
The approximation to the variance covariance matrix for the model parameters
}
\description{
Marginal Quasi-Likelihood approach to approximate the variance covariance matrix for the
model parameters
}
\keyword{internal}
