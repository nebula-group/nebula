context("data and parameter inputs")

#TODO test for binary data being of 2 types
#TODO test for continuous input numeric
#TODO test for E input - four columns and appropriate entries

set.seed(123)

#make smaller data for testing purposes
data = lapply(colon$modal, function(x) x[1:100,])

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

test_that("Check for missing data", {
  #make vector of number of rows in each matrix M of input data
  colon_missing <- colon$modal
  colon_missing[[2]][1,1] <- NA
  #should error if there are not
  expect_error(nebula(
    data = colon_missing,
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
  #should error if there are not the right number
  expect_error(nebula(
    data = data,
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
  #or if modtype is missing
  expect_error(nebula(
    data = data,
   # modtype = c(0), #
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
  expect_error(nebula(
    data = data,
    # modtype = c(0), #
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
  #or incorrect entry
  expect_error(nebula(
    data = data,
    modtype = c("continuous", "binary"), #
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
    data = data,
    modtype = c(0, 1), #
    E = colon$network,
    H = 3,
    modeta = c(1), # wrong number of parameters
    nu = 1,
    alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
  ))
  # or missing
  expect_error(nebula(
    data = data,
    modtype = c(0, 1), #
    E = colon$network,
    H = 3,
   # modeta = c(1),# missing modet
    nu = 1,
    alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
  ))
})

test_that("H input", {
  # #should error if there are not correct input
  # expect_error(nebula(
  #   data = data,
  #   modtype = c(0, 1), #
  #   E = colon$network,
  #   H = 3,
  #   modeta = c(0,1),
  #   nu = 1,
  #   alpha = 1,
  #   lam = 1,
  #   alpha_sigma = 10,
  #   beta_sigma = 10,
  #   alpha_p = 1,
  #   beta_p = 1,
  # ))
  # or missing
  expect_error(nebula(
    data = data,
    modtype = c(0, 1), #
    E = colon$network,
    #H = 3, #missing
    modeta = c(0, 1),
    nu = 1,
    alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
  ))
})

test_that("nu input", {
  # #should error if there are not correct input
  # expect_error(nebula(
  #   data = data,
  #   modtype = c(0, 1), #
  #   E = colon$network,
  #   H = 3,
  #   modeta = c(0,1),
  #   nu = 1,
  #   alpha = 1,
  #   lam = 1,
  #   alpha_sigma = 10,
  #   beta_sigma = 10,
  #   alpha_p = 1,
  #   beta_p = 1,
  # ))
  # or missing
  expect_error(nebula(
    data = data,
    modtype = c(0, 1), #
    E = colon$network,
    H = 3,
    modeta = c(0, 1),
    #nu = 1, #missing
    alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
  ))
})

test_that("alpha input", {
  # #should error if there are not correct input
  # expect_error(nebula(
  #   data = data,
  #   modtype = c(0, 1), #
  #   E = colon$network,
  #   H = 3,
  #   modeta = c(0,1),
  #   nu = 1,
  #   alpha = 1,
  #   lam = 1,
  #   alpha_sigma = 10,
  #   beta_sigma = 10,
  #   alpha_p = 1,
  #   beta_p = 1,
  # ))
  # or missing
  expect_error(nebula(
    data = data,
    modtype = c(0, 1), #
    E = colon$network,
    H = 3,
    modeta = c(0, 1),
    nu = 1, #missing
    #alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
  ))
})

test_that("lam input", {
  # #should error if there are not correct input
  # expect_error(nebula(
  #   data = data,
  #   modtype = c(0, 1), #
  #   E = colon$network,
  #   H = 3,
  #   modeta = c(0,1),
  #   nu = 1,
  #   alpha = 1,
  #   lam = 1,
  #   alpha_sigma = 10,
  #   beta_sigma = 10,
  #   alpha_p = 1,
  #   beta_p = 1,
  # ))
  # or missing
  expect_error(nebula(
    data = data,
    modtype = c(0, 1), #
    E = colon$network,
    H = 3,
    modeta = c(0, 1),
    nu = 1, #missing
    alpha = 1,
    #lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
  ))
})

test_that("pr0 input", {
  # #should error if there are not correct input
  expect_error(nebula(
    data = data,
    modtype = c(0, 1), #
    E = colon$network,
    H = 3,
    modeta = c(0,1),
    nu = 1,
    alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
    pr0 = -1
  ))
  # or missing
  expect_error(nebula(
    data = data,
    modtype = c(0, 1), #
    E = colon$network,
    H = 3,
    modeta = c(0, 1),
    nu = 1, #missing
    alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
    pr0 = 2
  ))
})

test_that("sig0 input", {
  # #should error if there are not correct input
  expect_error(nebula(
    data = data,
    modtype = c(0, 1), #
    E = colon$network,
    H = 3,
    modeta = c(0,1),
    nu = 1,
    alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
    sig0 = -1
  ))
})

test_that("binit input", {
  n_bad <- 99
  H <- 3
  binit_bad <- matrix(stats::rnorm(n_bad*H,0,0.01),n_bad,H)
  # #should error if there are not correct input
  expect_error(nebula(
    data = data,
    modtype = c(0, 1), #
    E = colon$network,
    H = 3,
    modeta = c(0,1),
    nu = 1,
    alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1,
    binit = binit_bad
  ))
})

test_that("E input", {
  e_bad <- colon$network[,1:3]
  # #should error if there are not correct input
  expect_error(nebula(
    data = data,
    modtype = c(0, 1), #
    E = e_bad,
    H = 3,
    modeta = c(0,1),
    nu = 1,
    alpha = 1,
    lam = 1,
    alpha_sigma = 10,
    beta_sigma = 10,
    alpha_p = 1,
    beta_p = 1
  ))
})

