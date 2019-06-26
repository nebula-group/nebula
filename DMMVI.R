# Variational Inference for Dirichlet Mixture Model
# written by Changgee Chang
# ver. 20191126

DMMVI <- function(X,type,A,H,eta,nu,alpha,lam,mu0,sig0,pr0,alpha_sigma,beta_sigma,alpha_p,beta_p,binit=NULL)
{
  n = nrow(X)
  p = ncol(X)
  idx0 = which(type==0)
  idx1 = which(type==1)
  p0 = length(idx0)
  p1 = length(idx1)

  # Initialize
  if ( is.null(binit) )
    b = matrix(rnorm(n*H,0,0.01),n,H)
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
    tmp = diffinv(-sEI)
    g = alpha + tmp[-1] - tmp[H+1]
    
    tmp = digamma(f+g)
    Elw = digamma(f) - tmp
    El1w = digamma(g) - tmp
    Elw[H] = 0
    El1w[H] = -Inf
    cumEl1w = diffinv(El1w)

    
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
      cbase[idx0,h] = cbase[idx0,h] + t(lgau0) %*% EI[,h]/2
      cbase[idx1,h] = t(Elbin[,idx1,h]) %*% EI[,h]
      cbase[idx1,h] = cbase[idx1,h] - t(lbin0) %*% EI[,h]
    }
    cbase = cbase - eta
    
    while (TRUE)
    {
      pc = c
      for ( j in 1:p )
      {
        c[j,] = cbase[j,] + nu*A[j,]%*%(2*Egam-1)
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
  
  list(EI=EI,Egam=Egam,m=m,iter=iter)
}

