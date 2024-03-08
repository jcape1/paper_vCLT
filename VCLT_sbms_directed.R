#######################################################################

## VCLT -- generating directed SBMs with independent left/right factors

#######################################################################

# NB: memberships are generated in an independent fashion,
# i.e., Z_{i} independent of Y_{i} for each index value i

## (Automated)
  
  memb.index.left <- sample(1:K, n, replace=TRUE, prob=pi.vec.left)
  memb.index.right <- sample(1:K, n, replace=TRUE, prob=pi.vec.right)
  
  block.sizes.left <- table(memb.index.left)
  block.sizes.right <- table(memb.index.right)
  
  Z.mtx.left <- diag(K)[memb.index.left,]
  Z.mtx.right <- diag(K)[memb.index.right,]
  
  
  # NB: Avoid saving edge probability matrix to reduce memory space
  # my.graph <- runif(n * n) <= Z.mtx.left %*% (rho*pref.mtx) %*% t(Z.mtx.right)
  # NB: the above graph is actually a logical adjacency matrix
  
  
  
