#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat fcb3(arma::mat basisvar, arma::mat basislag, arma::vec lag) {
  
  int nrow = basisvar.n_rows;
  int ncol = basisvar.n_cols * basislag.n_cols;
  
  int max_lag = max(lag);
  
  arma::mat crossbasis(ncol, nrow);
  crossbasis = crossbasis.fill(NA_REAL); 
  
  for (int v = (max_lag); v < (basisvar.n_rows); v++) {
    crossbasis.col(v) =   arma::vectorise(basislag.t() * arma::flipud(basisvar.rows(v-max_lag, v)), 0);
  }
  
  return(crossbasis.t());
  
}






/*** R
################################################################################
# TEST ON ALTERNATIVE CROSS-BASIS COMPUTATION
################################################################################

# LOAD THE PACKAGE
library(dlnm)
  library(splines)
  library(tsModel)
  
# DEFINE DIMENSIONS
  n <- 20000
l <- 300
vardf <- 5
lagdf <- 5

# CREATE OBJECTS
x <- rnorm(n)
  lag <- c(0,l)
  basisvar <- ns(x, df=vardf)
  basislag <- ns(lag[1]:lag[2],df=lagdf,int=T)
  
################################################################################
# CREATE FUNCTIONS
  
# FUNCTION REPLICATING dlnm VERSION 2.3.4
  fcb1 <- function(basisvar, basilag, lag) {
    crossbasis <- matrix(NA_real_, nrow=nrow(basisvar), ncol=ncol(basisvar)*ncol(basislag))
    for(v in seq(length=ncol(basisvar))) {
      mat <- as.matrix(Lag(basisvar[, v], lag[1]:lag[2]))
      for(l in seq(length=ncol(basislag))) {
        crossbasis[,ncol(basislag)*(v-1)+l] <- mat %*% (basislag[,l])
      }
    }
    crossbasis
  }
  
# ALTERNATIVE FUNCTION (MODIFIED FROM SUGGESTION BY JONATHAN KENNEL)
  fcb2 <- function(basisvar,basilag,lag) {
    crossbasis <- matrix(NA_real_, nrow=nrow(basisvar), ncol=ncol(basisvar)*ncol(basislag))
    for(v in (max(lag)+1):nrow(basisvar)) {
      mat <- basisvar[v:(v-(max(lag))),]
      crossbasis[v,] <- t(basislag) %*% mat
    }
    crossbasis
  }
  
  fftfcb <- function(basisvar, baislag) {
    
    convolve(basisvar[,1], rev(basislag[,1]), type = 'o')
    
  }
################################################################################
# TEST RESULTS AND TIMING
  
# COMPUTE CROSS-BASIS AND TEST TIME
  system.time(cb <- crossbasis(x,lag,list(df=5),list(df=5)))
    system.time(cb1 <- fcb1(basisvar,basislag,lag))
    system.time(cb2 <- fcb2(basisvar,basislag,lag))
    system.time(cb3 <- fcb3(basisvar,basislag,lag))
    
    system.time(tmp <- fftfcb(basisvar = basisvar, basislag))
    
    plot(tmp[-(length(tmp):(length(tmp)-l+1))], type='l')
    points(cb3[,1], type='l', col='red')
    tail(tmp[-(length(tmp):(length(tmp)-l+1))])
    tail(cb3[,1])
    
    
    all.equal(cb1, cb2)
    all.equal(cb1, cb3)
    all.equal(cb2, cb3)
    
    
    */
