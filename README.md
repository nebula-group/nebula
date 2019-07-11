
<!-- badges: start -->
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/margarethannum/PanCancer?branch=master&svg=true)](https://ci.appveyor.com/project/margarethannum/PanCancer)
[![Travis build status](https://travis-ci.org/margarethannum/PanCancer.svg?branch=master)](https://travis-ci.org/margarethannum/PanCancer) 
<!-- badges: end -->


# Variational Inference for Dirichlet Mixture Model with Incorporation of Graphical Knowledge

Practically, this can be used to incorporate biological networks/pathways to inform clustering solutions. 

Flexible sparsity parameters for multiple input data types allows for control over which data types need sparse vs rich feature selection. 

## Installation

You can install the production version of the package with: 

``` r
# install.packages("remotes")
remotes::install_github("changgee/PanCancer")
```

See documentation for description of parameters:

``` r
?PanCancer::DMMVI
```

Arguments description:

X: n by p matrix where n is the sample size and p is the feature size

type: p-dimensional vector of feature types. Currently supports continuous(=0) and binary(=1)

A: the adjacency matrix of the graph representing the graphical structure of X

H: the number of clusters to be fit

eta: sparsity parameter for gamma's

nu: smoothness parameter for gamma's

alpha: concentration parameter for dirichlet process

mu0: mean of the non-selected continuous features

sig0: variance of the non-selected continuous features

pr0: "on" probability of the non-selected binary features

alpha_sigma: shape parameter of the prior of residual variance(sigma^2)

beta_sigma: rate parameter of the prior of residual variance(sigma^2)

alpha_p: first shape parameter of the prior of the "on" probabilities(p_hj) of binary features

beta_p: second shape parameter of the prior of the "on" probabilities(p_hj) of binary features

binit: n by H initial matrix of B, exp(B_ih) is proportional to Pr(z_i=h). If NULL (default), random numbers are filled in.


