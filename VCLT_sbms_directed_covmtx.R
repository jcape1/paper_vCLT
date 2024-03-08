##########################################################

## VCLT -- directed SBMs -- asymptotic covariance matrices

##########################################################

  B.decomp <- svd(pref.mtx)
  B.boolean <- B.decomp$d > 1e-10
  #Ipq <- diag(sign(B.decomp$values[B.boolean]))
  
  B.svals <- B.decomp$d[B.boolean]
  B.svecs.left <- as.matrix(B.decomp$u[,B.boolean])
  B.svecs.right <- as.matrix(B.decomp$v[,B.boolean])
  
  T.mtx.left <- diag(as.vector(B.svals^(1/2)), nrow=d) %*% t(B.svecs.left)
  T.mtx.right <- diag(as.vector(B.svals^(1/2)), nrow=d) %*% t(B.svecs.right)
    
  delta.mtx.left <- T.mtx.left %*% diag(pi.vec.left) %*% t(T.mtx.left)
  delta.mtx.right <- T.mtx.right %*% diag(pi.vec.right) %*% t(T.mtx.right)
  
  adj.mtx.left <- solve( sqrtm(solve(delta.mtx.left)) %*% T.mtx.left %*% diag(sqrt(pi.vec.left)) )
  adj.mtx.right <- solve( sqrtm(solve(delta.mtx.right)) %*% T.mtx.right %*% diag(sqrt(pi.vec.right)) )
  
################

  # Left latent variable -- unscaled ASE CLT for directed SBMs
  cov_mtx_dSBM_evec_left <- function(vec){
    sqrtm(solve(delta.mtx.left)) %*% solve(delta.mtx.right) %*%
      Reduce(`+`, lapply(1:K, function(xxx)
        pi.vec.right[xxx] *
          (vec %*% (T.mtx.right %*% diag(K)[xxx,]))[1] *
          (1 - (vec %*% (T.mtx.right %*% diag(K)[xxx,])))[1] *
          ((T.mtx.right %*% diag(K)[xxx,]) %*% t(T.mtx.right %*% diag(K)[xxx,])))) %*% 
      solve(delta.mtx.right) %*% sqrtm(solve(delta.mtx.left))
  }
  
  # Left latent variable -- asymptotic covariance for varimax embedding
  cov_mtx_vmx_dSBM_left <- function(vec){
    adj.mtx.left %*% cov_mtx_dSBM_evec_left(vec) %*% t(adj.mtx.left)
  }

################  
  
  # Right latent variable -- unscaled ASE CLT for directed SBMs
  cov_mtx_dSBM_evec_right <- function(vec){
    sqrtm(solve(delta.mtx.right)) %*% solve(delta.mtx.left) %*%
      Reduce(`+`, lapply(1:K, function(xxx)
        pi.vec.left[xxx] *
          (vec %*% (T.mtx.left %*% diag(K)[xxx,]))[1] *
          (1 - (vec %*% (T.mtx.left %*% diag(K)[xxx,])))[1] *
          ((T.mtx.left %*% diag(K)[xxx,]) %*% t(T.mtx.left %*% diag(K)[xxx,])))) %*% 
      solve(delta.mtx.left) %*% sqrtm(solve(delta.mtx.right))
  }
  
  # Right latent variable -- asymptotic covariance for varimax embedding
  cov_mtx_vmx_dSBM_right <- function(vec){
    adj.mtx.right %*% cov_mtx_dSBM_evec_right(vec) %*% t(adj.mtx.right)
  }
  
  