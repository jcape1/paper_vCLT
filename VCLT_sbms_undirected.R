
## SBM generative procedure - part II

  ## (Automated) -- block-structured for use with function sample_sbm
  
  block.sizes <- table(sample(1:K, n, replace=TRUE, prob=pi.vec))
  
  Z.mtx <- t(do.call(cbind, lapply(1:K, function(x) replicate(block.sizes[x], diag(K)[x,]))))
  # NB: need standardized variance for identifiability; acts as DCSBM scaling.
  # NB: order blocks for simplicity; need to change for directed SBM code.
  
  memb.index <- unlist(lapply(1:K, function(x) replicate(block.sizes[x], x)))
  
  my.graph <- sample_sbm(n,
                         rho * pref.mtx,
                         block.sizes,
                         directed = !isSymmetric(pref.mtx),
                         loops = TRUE)
