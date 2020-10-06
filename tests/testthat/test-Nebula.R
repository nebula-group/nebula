context("data input")

#TODO add test for correct modtype input (0 or 1)
#TODO add test for missing modtype, E, H, modeta, nu, alpha, lam
#TODO missing data test
#TODO test for binary data being of 2 types
#TODO test for continuous input numeric

set.seed(123)


test_that("M modalities have same number of samples", {
  #make vector of number of rows in each matrix M of input data
  colon_bad_modality <- colon$modal
  colon_bad_modality[[1]] <- colon_bad_modality[[1]][1:330,]
  #should error if there are not
  expect_error(nebula(
    data = colon_bad_modality,
    modtype = c(0, 1), #
    E = colon$network,
    H = 3,
    modeta = c(1, 0.2),
    nu = 1,
    alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
  ))
})

test_that("modtype incorrect", {
  #make vector of number of rows in each matrix M of input data
  #should error if there are not
  expect_error(nebula(
    data = colon$modal,
    modtype = c(0), #
    E = colon$network,
    H = 3,
    modeta = c(1, 0.2),
    nu = 1,
    alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
  ))
})

test_that("modeta incorrect", {
  #make vector of number of rows in each matrix M of input data
  #should error if there are not
  expect_error(nebula(
    data = colon$modal,
    modtype = c(0, 1), #
    E = colon$network,
    H = 3,
    modeta = c(1),
    nu = 1,
    alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
  ))
})


