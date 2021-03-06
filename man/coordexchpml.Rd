% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordexchpml.R
\name{coordexchpml}
\alias{coordexchpml}
\title{Search for an optimal design under the panel mixed logit model}
\usage{
coordexchpml(
  nChoiceSet,
  nAlternative,
  nLevelAttribute,
  effectMean,
  effectVar,
  approx = "pql",
  nStartDesign = 10,
  optcrit = "D"
)
}
\arguments{
\item{nChoiceSet}{Number of choice sets}

\item{nAlternative}{The number of alternatives in each choice set}

\item{nLevelAttribute}{Vector where each element is the number of levels for a factor}

\item{effectMean}{Vector of means for the effects coded attribute effects}

\item{effectVar}{Vector of variances for the effects coded attribute effects}

\item{approx}{Type of approximation to use (MQL, PQL, Importance, Laplace, or MSM)}

\item{nStartDesign}{Number of times to run the coordinate exchange algorithm}

\item{optcrit}{Optimality criterion. Either "A" or "D"}
}
\value{
Object containing the best design found by the coordinate exchange algorithm
}
\description{
Search for an optimal design under the panel mixed logit model
}
