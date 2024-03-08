
## SBM generative procedure - example based on Section 5.2 in Rohe SPCA arXiv paper v1
  
  n <- 5000
  rho <- 1
  
  pi.vec <- rep(1/4,4)
  
  aa <- 0.60; bb <- 0.20; cc <- 0.10; dd <- 0.10
  ee <- 0.70; ff <- 0.05; gg <- 0.05;
  hh <- 0.60; ii <- 0.25;
  jj <- 0.60;
  
  pref.mtx <- rbind(c(aa, bb, cc, dd),
                    c(bb, ee, ff, gg),
                    c(cc, ff, hh, ii),
                    c(dd, gg, ii, jj))
  
  d <- sum(abs(eigen(pref.mtx)$values) > 1e-10)
  
  K <- dim(pref.mtx)[1]