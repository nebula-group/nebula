#' NebulaCore
#'
#' Network based latent dirichlet subtype analysis
#' (ver. 20190920)
#'
#' @param X n by p matrix where n is the sample size and p is the feature size
#' @param type p-dimensional vector of feature types. Currently supports continuous(=0) and binary(=1)
#' @param E e by 2 matrix of the graph edges representing the graphical structure of X
#' @param H the number of clusters to be fit
#' @param eta length p vector of sparsity parameters
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
#' @export
#' @author Changgee Chang
#' @examples
#' # ADD EXAMPLE

NebulaCore <- function(X,type,E,H,eta,nu,alpha,lam,alpha_sigma,beta_sigma,alpha_p,beta_p,mu0,sig0,pr0,binit=NULL)
{
  n = nrow(X)
  p = ncol(X)
  idx0 = which(type==0)
  idx1 = which(type==1)
  p0 = length(idx0)
  p1 = length(idx1)

  E = E[order(E[,1],E[,2]),]
  e = nrow(E)
  nadj = rep(0,p)
  tmp = c(diff(E[,1])!=0,TRUE)
  nadj[E[tmp,1]] = diff(c(0,which(tmp)))
  Eidx = stats::diffinv(nadj)

  # Initialize
  if ( is.null(binit) )
    b = matrix(stats::rnorm(n*H,0,0.01),n,H)
  else
    b = binit
  b = b - apply(b,1,max)
  eb = exp(b)
  EI = eb/apply(eb,1,sum)

  c = matrix(0,p,H)
  cbase = matrix(0,p,H)
  Egam = 1/(1+exp(-c))

  Elgau = array(0,c(n,p,H))
  Elbin = array(0,c(n,p,H))
  lgau0 = (X[,idx0] - mu0)^2/sig0 + log(sig0)
  lbin0 = X[,idx1]*log(pr0) + (1-X[,idx1])*log(1-pr0)

  iter = 0
  while (TRUE)
  {
    iter = iter + 1

    # Step 1
    sEI = apply(EI,2,sum)
    f = 1 + sEI
    tmp = stats::diffinv(-sEI)
    g = alpha + tmp[-1] - tmp[H+1]

    tmp = digamma(f+g)
    Elw = digamma(f) - tmp
    El1w = digamma(g) - tmp
    Elw[H] = 0
    El1w[H] = -Inf
    cumEl1w = stats::diffinv(El1w)


    # Step 2
    EgamsEI = Egam*rep(sEI,each=p)
    XEI = t(X)%*%EI
    EgamXEI = Egam*XEI
    d = alpha_sigma + EgamsEI
    v = lam + EgamsEI
    m = EgamXEI/v
    r = beta_sigma - v*m^2 + Egam*(t(X^2)%*%EI)

    Elsig = log(r/2) - digamma(d/2)
    for ( h in 1:H )
      Elgau[,,h] = t((t(X) - m[,h])^2*d[,h]/r[,h] + 1/v[,h] + Elsig[,h])


    # Step 3
    s = alpha_p + EgamXEI
    t = beta_p + EgamsEI - EgamXEI
    s[idx0,] = 1
    t[idx0,] = 1

    tmp = digamma(s+t)
    Elp = digamma(s) - tmp
    El1p = digamma(t) - tmp
    for ( h in 1:H )
      Elbin[,,h] = t(t(X)*Elp[,h]+(1-t(X))*El1p[,h])


    # Step 4
    for ( h in 1:H )
    {
      cbase[idx0,h] = - t(Elgau[,idx0,h]) %*% EI[,h]/2
      cbase[idx0,h] = cbase[idx0,h] + t(lgau0) %*% EI[,h]/2 - eta[idx0]
      cbase[idx1,h] = t(Elbin[,idx1,h]) %*% EI[,h]
      cbase[idx1,h] = cbase[idx1,h] - t(lbin0) %*% EI[,h] - eta[idx1]
    }

    Egam[nadj==0,] = 1/(1+exp(-cbase[nadj==0,]))
    while (TRUE)
    {
      pc = c
      for ( j in which(nadj!=0) )
      {
        if ( nadj[j] == 1 )
          c[j,] = cbase[j,] + nu*(2*Egam[E[Eidx[j+1],2],]-1)
        else
          c[j,] = cbase[j,] + nu*apply(2*Egam[E[(Eidx[j]+1):Eidx[j+1],2],]-1,2,sum)
        Egam[j,] = 1/(1+exp(-c[j,]))
      }
      if ( max(abs(pc-c)) < 1e-4 )
        break
    }


    # Step 5
    pb = b
    for ( h in 1:H )
    {
      b[,h] = - Elgau[,idx0,h] %*% Egam[idx0,h]/2
      b[,h] = b[,h] - lgau0 %*% (1-Egam[idx0,h])/2
      b[,h] = b[,h] + Elbin[,idx1,h] %*% Egam[idx1,h]
      b[,h] = b[,h] + lbin0 %*% (1-Egam[idx1,h])
      b[,h] = b[,h] + Elw[h] + cumEl1w[h]
    }
    b = b - apply(b,1,max)
    eb = exp(b)
    EI = eb/apply(eb,1,sum)

    if ( max(abs(pb-b)) < 1e-4 )
      break
  }

  list(EI=EI,Egam=Egam,m=m,lpr=Elp,iter=iter)
}

