#' Nebula
#'
#' Network based latent dirichlet subtype analysis
#' (ver. 20190920)
#'
#' @param data list of M data matrices, where each matrix is n samples by p_m features for modality M
#' @param modtype M-length vector of feature types for M modalities. Currently supports continuous(=0) and binary(=1)
#' @param E e by 4 matrix, each tuple (row) of which represents an edge; (m1,j1,m2,j2) variable j1 of modality m1 is connected to variable j2 of modality m2
#' @param H the number of clusters to be fit
#' @param modeta length M vector of sparsity parameters for M modalities
#' @param nu smoothness parameter for gamma's
#' @param alpha concentration parameter for dirichlet process
#' @param lam shrinkage parameter for means of selected continuous features
#' @param alpha_sigma shape parameter of the prior of residual variance(sigma^2)
#' @param beta_sigma rate parameter of the prior of residual variance(sigma^2)
#' @param alpha_p first shape parameter of the prior of the 'active' probabilities(p_hj) of binary features
#' @param beta_p second shape parameter of the prior of the 'active' probabilities(p_hj) of binary features
#' @param mu0 mean of the non-selected continuous features
#' @param sig0 variance of the non-selected continuous features
#' @param pr0 'active' probability of the non-selected binary features
#' @param binit n by H initial matrix of B, exp(B_ih) is proportional to Pr(z_i=h). If NULL (default), random numbers are filled in.
#'
#' @return A list containing clustering assignments, variable selection, and
#'   posterior probabilities
#' \itemize{
#' \item clustering cluster assignment
#' \item defvar list of M matrices; each matrix is p_m by H indicating the variable j in modality m is a defining variable for the cluster h.
#' \item clus_pr n by H matrix containing the probability that the subject i belongs to the cluster h.
#' \item defvar_pr list of M matrices; each matrix is p_m by H containing the probability that the variable j in modality m is a defining variable for the cluster h.
#' \item def_m list of M matrices; each matrix is p_m by H containing the mean of the variable j as a defining variable for the cluster h. (continuous variable only)
#' \item def_lpr list of M matrices: each matrix is p_m by H containing the log probabilities of the variable j being 'active' ad a defining variable for the cluster h. (binary variable only)
#' \item iter the number of iterations until the algorithm converges.
#' }
#' @export
#' @author Changgee Chang
#' @examples
#' # ADD EXAMPLE
#' # source("NebulaCore.R") No need to source if function is also in our package
Nebula <- function(data, modtype, E, H, modeta, nu, alpha, lam, alpha_sigma, beta_sigma, alpha_p, beta_p, mu0 = 0, sig0 = 20, pr0 = 0.5, binit = NULL) {
  M <- length(data)

  if (length(modtype) != M) {
    stop("invalid length of modtype")
  }
  if (length(modeta) != M) {
    stop("invalid length of modeta")
  }

  p <- rep(0, M)
  for (m in 1:M)
  {
    if (m == 1) {
      n <- nrow(data[[1]])
    } else if (nrow(data[[m]]) != n) {
      stop("different size of data matrices")
    }
    p[m] <- ncol(data[[m]])
  }
  P <- sum(p)
  cump <- stats::diffinv(p)

  X <- matrix(0, n, P)
  mod <- matrix(FALSE, P, M) # modality
  type <- rep(0, P) # variable type
  if (ncol(E) != 4) {
    stop("Matrix E must have 4 columns")
  }
  E <- E[, c(2, 4)] + cump[E[, c(1, 3)]]
  eta <- rep(0, P) # sparsity parameter for each variable
  for (m in 1:M)
  {
    X[, (cump[m] + 1):cump[m + 1]] <- data[[m]]
    mod[(cump[m] + 1):cump[m + 1], m] <- TRUE
    if (modtype[m] != 0 & modtype[m] != 1) {
      stop("invalid value for modality data type")
    }
    type[(cump[m] + 1):cump[m + 1]] <- modtype[m]
    eta[(cump[m] + 1):cump[m + 1]] <- modeta[m]
  }

  out <- NebulaCore(X, type, E[, c(2, 4)], H, eta, nu, alpha, lam, alpha_sigma, beta_sigma, alpha_p, beta_p, mu0, sig0, pr0, binit)

  defvar <- list()
  defvar_pr <- list()
  def_m <- list()
  def_lpr <- list()
  for (m in 1:M)
  {
    defvar[[m]] <- matrix(out$Egam[mod[, m], ] > 0.5, ncol = H)
    defvar_pr[[m]] <- matrix(out$Egam[mod[, m], ], ncol = H)
    if (modtype[m] != 0) {
      out$m[mod[, m], ] <- NA
    }
    def_m[[m]] <- matrix(out$m[mod[, m], ], ncol = H)
    if (modtype[m] != 1) {
      out$lpr[mod[, m], ] <- NA
    }
    def_lpr[[m]] <- matrix(out$lpr[mod[, m], ], ncol = H)
  }

  #
  # Outputs
  #
  # clustering: cluster assignment
  # defvar: list of M matrices; each matrix is p_m by H indicating the variable j in modality m is a defining variable for the cluster h.
  # clus_pr: n by H matrix containing the probability that the subject i belongs to the cluster h.
  # defvar_pr: list of M matrices; each matrix is p_m by H containing the probability that the variable j in modality m is a defining variable for the cluster h.
  # def_m: list of M matrices; each matrix is p_m by H containing the mean of the variable j as a defining variable for the cluster h. (continuous variable only)
  # def_lpr: list of M matrices: each matrix is p_m by H containing the log probabilities of the variable j being 'active' ad a defining variable for the cluster h. (binary variable only)
  # iter: the number of iterations until the algorithm converges.
  list(clustering = apply(out$EI, 1, which.max), defvar = defvar, clus_pr = out$EI, defvar_pr = defvar_pr, def_m = def_m, def_lpr = def_lpr, iter = out$siter)
}
