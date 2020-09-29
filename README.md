
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nebula

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/changgee/Nebula?branch=master&svg=true)](https://ci.appveyor.com/project/changgee/Nebula)
[![Travis build
status](https://travis-ci.org/margarethannum/Nebula.svg?branch=master)](https://travis-ci.org/margarethannum/Nebula)
<!-- badges: end -->

The Nebula package implements Network-based latent-Dirichlet subtype
analysis (Nebula) algorithm.

Practically, this can be used to incorporate biological
networks/pathways to inform clustering solutions.

Flexible sparsity parameters for multiple input data types allows for
control over which data types need sparse vs rich feature selection.

## Installation

You can install the development version from
[GitHub](https://github.com/nebula-group/nebula) with:

``` r
# install.packages("remotes")

remotes::install_github("nebula-group/nebula", build_vignettes = TRUE)

library(nebula)
```

Learn more in `vignette("nebula_tutorial", package = "nebula")` or
`?nebula`.
