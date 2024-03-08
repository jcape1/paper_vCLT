
## VCLT -- Preamble

  library(igraph) # package for network analysis & random graphs
  library(irlba) # package for numerical linear algebra routines
  library(mclust) # package for Gaussian mixture model clustering
  
  #library(plot3D) # option for 3D plotting
  #library(rgl); library(car) # option for 3d plotting
  
  library(plotrix)
  library(RSpectra)
  library(MASS)
  
  library(expm) # to compute matrix square roots
  
  sym <- function(mtx){
    mtx[lower.tri(mtx, diag=FALSE)] =
      t(mtx)[lower.tri(mtx, diag=FALSE)];
    mtx # << symmetrized
  }