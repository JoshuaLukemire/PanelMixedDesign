
# PanelMixedDesign <img src="https://img.shields.io/badge/status-dev-yellow" alt="status badge" align="right"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/JoshuaLukemire/PanelMixedDesign/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JoshuaLukemire/PanelMixedDesign/actions)
<!-- badges: end -->

The `PanelMixedDesign` R package provides fast approximations to the
variance-covariance matrix of model parameters for discrete choice
experiments modeled under the Panel Mixed Logit (PML) framework. These
approximations avoid computationally expensive numerical integration and
support searches for efficient experimental designs via a coordinate
exchange algorithm.

------------------------------------------------------------------------

## Installation

You can install the development version of `PanelMixedDesign` from
GitHub with:

\`\`\`r \# install.packages(“devtools”)
devtools::install_github(“JoshuaLukemire/PanelMixedDesign”)
