
## Asymptotic covariance matrices --  undirected graphs

  B.decomp <- eigen(pref.mtx)
  B.boolean <- abs(B.decomp$values) > 1e-10
  Ipq <- diag(sign(B.decomp$values[B.boolean]))
  
  B.evals <- B.decomp$values[B.boolean]
  B.evecs <- as.matrix(B.decomp$vectors[,B.boolean])
      nonzero.row.index.B <- min(which(sapply(1:K, function(xxx) min(abs(B.evecs[xxx,])) > 1e-10) == TRUE))
  B.evecs <- B.evecs %*% diag(sign(B.evecs[nonzero.row.index.B,]), nrow=d)
  # NB: set choice of sign.
  
  T.mtx <- diag(as.vector(abs(B.evals)^(1/2)), nrow=d) %*% t(B.evecs)
  # NB: columns of T.mtx are GRDPG latent positions for full rank SBM setting.
  
  delta.mtx <- T.mtx %*% diag(pi.vec) %*% t(T.mtx)
  
  adj.mtx <- diag(1/sqrt(pi.vec)) %*% solve(T.mtx)
  # NB: linear transformation (adjustment matrix) from GRDPGs.
  # NB: requires full rank transform.
  
  # FYI GRDPG model with latent positions T.mtx %*% Z_i where Z_i is a standard basis vector.
  
  # GRDPG ASE asymptotic covariance for SBMs
  cov_mtx_GRDPG_ASE <- function(vec){
    Ipq %*% solve(delta.mtx) %*%
      Reduce(`+`, lapply(1:K, function(xxx)
        pi.vec[xxx] *
          (vec %*% Ipq %*% (T.mtx %*% diag(K)[xxx,]))[1] *
          (1 - (vec %*% Ipq %*% (T.mtx %*% diag(K)[xxx,])))[1] *
          ((T.mtx %*% diag(K)[xxx,]) %*% t(T.mtx %*% diag(K)[xxx,])))) %*% 
      solve(delta.mtx) %*% Ipq
  }
  
  
  # Varimax asymptotic covariance for SBMs.
  cov_mtx_vmx <- function(vec){
    adj.mtx %*% cov_mtx_GRDPG_ASE(vec) %*% t(adj.mtx)
  }
  
  
  # GRDPG evec asymptotic covariance for SBMs.
  cov_mtx_GRDPG_evec <- function(vec){
    solve(delta.mtx)^(1/2) %*% cov_mtx_GRDPG_ASE(vec) %*% t(solve(delta.mtx)^(1/2))
  }
  
  #lapply(1:K, function(xxx) cov_mtx_GRDPG_ASE(as.vector(T.mtx %*% diag(K)[xxx,])))
  
  #lapply(1:K, function(xxx) cov_mtx_vmx(as.vector(T.mtx %*% diag(K)[xxx,])))
  
  # sparse versions
  cov_mtx_GRDPG_ASE_sparse <- function(vec){
    Ipq %*% solve(delta.mtx) %*%
      Reduce(`+`, lapply(1:K, function(xxx)
        pi.vec[xxx] *
          (vec %*% Ipq %*% (T.mtx %*% diag(K)[xxx,]))[1] *
          ((T.mtx %*% diag(K)[xxx,]) %*% t(T.mtx %*% diag(K)[xxx,])))) %*% 
      solve(delta.mtx) %*% Ipq
  }
  #
    cov_mtx_vmx_sparse <- function(vec){
    adj.mtx %*% cov_mtx_GRDPG_ASE_sparse(vec) %*% t(adj.mtx)
  }
  #
    cov_mtx_GRDPG_evec_sparse <- function(vec){
    solve(delta.mtx)^(1/2) %*% cov_mtx_GRDPG_ASE_sparse(vec) %*% t(solve(delta.mtx)^(1/2))
    }
