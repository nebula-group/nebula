
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Nebula

<!-- badges: start -->

<!-- badges: end -->

The Nebula package implements Network-based latent-Dirichlet subtype
analysis algorithm.

Practically, this can be used to incorporate biological
networks/pathways to inform clustering solutions.

Flexible sparsity parameters for multiple input data types allows for
control over which data types need sparse vs rich feature selection.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("changgee/Nebula", ref = "development")
```

See documentation for description of parameters:

``` r
?Nebula::DMMVI
```

Arguments description:

X: n by p matrix where n is the sample size and p is the feature size

type: p-dimensional vector of feature types. Currently supports
continuous(=0) and binary(=1)

A: the adjacency matrix of the graph representing the graphical
structure of X

H: the number of clusters to be fit

eta: sparsity parameter for gamma’s

nu: smoothness parameter for gamma’s

alpha: concentration parameter for dirichlet process

mu0: mean of the non-selected continuous features

sig0: variance of the non-selected continuous features

pr0: “on” probability of the non-selected binary features

alpha\_sigma: shape parameter of the prior of residual variance(sigma^2)

beta\_sigma: rate parameter of the prior of residual variance(sigma^2)

alpha\_p: first shape parameter of the prior of the “on”
probabilities(p\_hj) of binary features

beta\_p: second shape parameter of the prior of the “on”
probabilities(p\_hj) of binary features

binit: n by H initial matrix of B, exp(B\_ih) is proportional to
Pr(z\_i=h). If NULL (default), random numbers are filled in.
adding a line
