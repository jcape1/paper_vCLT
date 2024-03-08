
## Simulation model -- three block DCSBM

    n <- 2000
    rho <- 1
    
    k.pre.def <- 3
    
    pi.vec.left <- c(1/k.pre.def, 1/k.pre.def, 1/k.pre.def) # i.e., vsp latent factor matrix "Z" probabilities
    pi.vec.right <- pi.vec.left
    
    aa <- 0.2; bb <- 0.1;
    
    pref.mtx <- bb*(rep(1,k.pre.def) %o% rep(1,k.pre.def)) + (aa-bb)*diag(k.pre.def)
    
    d <- sum(svd(pref.mtx)$d > 1e-10)
    
    K <- dim(pref.mtx)[1]
    
    memb.index.left <- sample(1:K, n, replace=TRUE, prob=pi.vec.left)
    memb.index.right <- memb.index.left
    
    block.sizes.left <- table(memb.index.left)
    block.sizes.right <- block.sizes.left
    
    theta.mean <- 1/2
    mmt.const <- sqrt((1/3)*(13/48))
    Z.mtx.left <- diag(runif(n, 1/4, 3/4)) %*% diag(K)[memb.index.left,] / mmt.const
    Z.mtx.right <- Z.mtx.left