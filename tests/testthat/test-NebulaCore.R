context("inner workings of core formula")

library(nebula)

set.seed(123)

#set all data as is done in shell nebula function
#make smaller data for testing purposes
data = lapply(colon$modal, function(x) x[1:100,])


modtype = c(0, 1)
E = colon$network
H = 3
modeta = c(1, 0.2)
nu = 1
alpha = 1
lam = 1
alpha_sigma = 10
beta_sigma = 10
alpha_p = 1
beta_p = 1

# now do some prep to the data as in nebula() function, get ready
# to pass to NebulaCore()
M <- length(data)

p <- rep(0, M)

for (m in 1:M)
{
  if (m == 1) {
    n <- nrow(data[[1]])
  } else if (nrow(data[[m]]) != n) {
    stop("data matrices need to have same number of samples (rows)")
  }
  p[m] <- ncol(data[[m]])
}

P <- sum(p)
cump <- stats::diffinv(p)

n <- nrow(data[[1]])

X <- matrix(0, n, P)
mod <- matrix(FALSE, P, M) # modality
type <- rep(0, P) # variable type

E2 <- E[, c(2, 4)] + cump[E[, c(1, 3)]]
eta <- rep(0, P) # sparsity parameter for each variable
for (m in 1:M)
{
  X[, (cump[m] + 1):cump[m + 1]] <- data[[m]]
  mod[(cump[m] + 1):cump[m + 1], m] <- TRUE
  type[(cump[m] + 1):cump[m + 1]] <- modtype[m]
  eta[(cump[m] + 1):cump[m + 1]] <- modeta[m]
}

# set other parameters
alpha_sigma = 1
beta_sigma = 1
alpha_p = 1
beta_p = 1
mu0 = 0
sig0 = 20
pr0 = 0.5
binit = NULL

out <- nebula:::NebulaCore(X, type, E2, H, eta, nu, alpha, lam, alpha_sigma, beta_sigma, alpha_p, beta_p, mu0, sig0, pr0, binit = NULL)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
