# required packages
library(matrixStats)
source('./scripts/helper.R')

# function to get the kernel matrix
kern <- function(x1, x2, nloc, distmat = NULL) {
  #' @param x1 first data-point of the covariate
  #' @param x2 second data-point of the covariate. The kernel will be calculated based on these two points.
  #' @param nloc the number of spatial points
  #' @param distmat distance matrix between the spatial points (if available)
  
  if (is.null(distmat)) {
    return(exp(-0.1 * sum((x1 - x2)^2)) * exp(-0.1 * outer(1:nloc, 1:nloc, FUN = "-")^2))
  } else {
    return(exp(-0.1 * sum((x1 - x2)^2)) * exp(-0.1 * distmat))
  }
} 


# function to estimate quantile from a given data
estimate.quantile <- function(Y, Xt, Z = NULL, quant = 0.5, x = colMedians(Xt), maxiter = 1000, tol = 1e-5) {
  
  #' @param Y data-matrix of the response variable
  #' @param Xt data-matrix of the covariate process 
  #' @param Z latitude-longitude information (if available)
  #' @param quant which quantile to work on, default is median (0.5)
  #' @param x value at which the estimate needs to be obtained, default is sample median
  #' @param maxiter maximum number of iterations to run the algorithm, default 1000
  #' @param tol tolerance for convergence
  
  # get the dimensions of the data
  ntime <- nrow(Y)
  nloc <- ncol(Y)
  
  # check if the data fits the required setting for theoretical gurantee
  alpha <- log(nloc) / log(ntime)
  if (alpha > 1/2) {
    warning("The number of locations are higher than expected, theoretical gurantees may not be accurate.")
  }
  
  # bandwidth specification
  beta <- 1.5 * alpha - 1
  nbn <- ntime^(-beta)
  
  eps <- 1e-10  # tolerance level for matrix inversion and conditioning
  
  # get distance matrix based on available info about locations
  if (is.null(Z)) {
    distmat <- NULL
  } else {
    distmat <- as.matrix(dist(Z))
  }

  # calculate the kernal values for all Xt
  kern.array <- lapply(c(1:ntime), FUN = function(j) kern(Xt[j,], x, nloc, distmat))
  kern.norm <- unlist(lapply(kern.array, FUN = function(x1) norm(x1)))
  
  # find the closest x-index for every observation
  best.index <- which.max(kern.norm)
  
  # initialization of the algorithm
  init.u <- ntime + runif(nloc)
  init.u <- init.u/inner.norm(init.u) * (2 * quant - 1)
  est.mean <- numeric(nloc)
  est.var <- matrix(0, nrow = nloc, ncol = nloc)
  
  # preliminary settings for the iterative algorithm
  niter <- 0
  converged <- FALSE
  err <- Inf
  q <- runif(nloc)
  
  # calculations that are needed several times
  idx = which((abs(1:ntime - best.index) < nbn) & (kern.norm > eps))  # apply the kernel to the support specified by the bandwidth
  kmats = kern.array[idx]
  Ymat = Y[idx,]
  mean.calc = lapply(1:length(idx),FUN = function(x1) (kmats[[x1]] %*% Ymat[x1,]))
  var.calc = lapply(1:length(idx),FUN = function(x1) tcrossprod(mean.calc[[x1]], Ymat[x1,]))
  
  wmat = Reduce('+', kmats)
  add.mean = Reduce('+', mean.calc)
  add.var = Reduce('+', var.calc)
  wmat.inv <- solve(wmat + eps*diag(nloc))
  
  # run the algorithm until convergence
  while (!converged) {
    
    est.mean = est.mean + add.mean
    est.var = est.var + add.var
    rhs.term = est.mean + est.var %*% init.u / 2
    qnew <- wmat.inv %*% rhs.term
    
    # error calculation
    err <- sqrt(sum((q - qnew)^2)) / sqrt(sum(q^2))
    q <- qnew
    niter <- niter + 1
    
    # convergence check
    if ((niter > maxiter) | (err < tol)) {
      converged <- TRUE
    }  
  }
  
  # return the estimated quantile
  return(list(
    "q" = qnew
  ))
}

estimate.sigma <- function(Y, Xt, Z = NULL, x = colMedians(Xt)) {
  
  #' @param Y data-matrix of the response variable
  #' @param Xt data-matrix of the covariate process 
  #' @param Z latitude-longitude information (if available)
  #' @param x value at which the estimate needs to be obtained, default is sample median
  
  # get the dimensions of the data and some initial parameters
  ntime <- nrow(Y)
  nloc <- ncol(Y)
  alpha <- log(nloc) / log(ntime)  
  beta <- 1.5 * alpha - 1
  nbn <- ntime^(-beta)
  
  if (is.null(Z)) {
    distmat <- NULL
  } else {
    distmat <- as.matrix(dist(Z))
  }
  kern.array = lapply(c(1:ntime), FUN = function(j) kern(Xt[j,], x, nloc, distmat))
  kern.norm <- unlist(lapply(kern.array, FUN = function(x1) norm(x1)))
  best.index <- which.max(kern.norm)
  idx = which((abs(1:ntime - best.index) < nbn) & (kern.norm > eps))  # apply the kernel to the support specified by the bandwidth
  
  # estimating sigma using nonparametric method using kernel smoothing
  var.calc <- lapply(1:length(idx),FUN = function(i) (kern.array[[idx[i]]] %*% tcrossprod(Y[idx[i],])) )
  var.est <- solve(Reduce('+', kern.array[idx]) + 1e-10*diag(nloc), Reduce('+', var.calc))  
  return(var.est)
}



# estimate the quantity eta_t (eq. (15)) given a dataset
eta_t <- function(q0, Y, Xt, Z = NULL, x = colMedians(Xt), normalize = TRUE) {
  #' @param q0 the quantile vector at which eta_t(q_0) is to be calculated 
  #' @param Y data-matrix of the response variable
  #' @param Xt data-matrix of the covariate process 
  #' @param Z latitude-longitude information (if available)
  #' @param x value at which the estimate needs to be obtained, default is sample median
  #' @param normalize should we divide by npb_n?
  
  # get the dimensions of the data and some initial parameters
  ntime <- nrow(Y)
  nloc <- ncol(Y)
  alpha <- log(nloc) / log(ntime)  
  beta <- 1.5 * alpha - 1
  nbn <- ntime^(-beta)
  eps <- 1e-10
  
  if (is.null(Z)) {
    distmat <- NULL
  } else {
    distmat <- as.matrix(dist(Z))
  }
  kern.array <- lapply(c(1:ntime), FUN = function(j) kern(Xt[j,], x, nloc, distmat))
  kern.norm <- unlist(lapply(kern.array, FUN = function(x1) norm(x1)))
  best.index <- which.max(kern.norm)
  idx = which((abs(1:ntime - best.index) < nbn) & (kern.norm > eps))  # apply the kernel to the support specified by the bandwidth
  
  mu0hat <- estimate.quantile(Y, Xt, Z, 0.5, x)[["q"]]
  muhat <- t(sapply(1:length(idx), FUN = function(i) {
    estimate.quantile(Y, Xt, Z, 0.5, Xt[idx[i], ])[["q"]]
  }))  # transposing ensures each row corresponds to a timepoint
  
  etavec <- t(sapply(1:length(idx), FUN = function(i) {
    mu_minus_q0 <- (mu0hat - q0)  # Simply (mu_0 - q_0)
    eta <- kern.array[[idx[i]]] %*% (muhat[i, ] - q0)  # eta(t_i) as in Eq. (12)
    return(normalize.vec(mu_minus_q0) - normalize.vec(eta))
  }))  # each row corresponds to a timepoint
  
  term.calc = lapply(1:length(idx),FUN = function(i) (kern.array[[idx[i]]] %*% etavec[i,] ))
  if (normalize) {
    return(Reduce('+', term.calc) / (nloc * nbn))
  } else {
    return(Reduce('+', term.calc))
  }
}

# estimate the quantity Psi_t (eq. (17)) given a dataset
Psi_t <- function(q0, Y, Xt, Z = NULL, x = colMedians(Xt), normalize = TRUE) {
  #' @param q0 the quantile vector at which Psi_t(q_0) is to be calculated 
  #' @param Y data-matrix of the response variable
  #' @param Xt data-matrix of the covariate process 
  #' @param Z latitude-longitude information (if available)
  #' @param x value at which the estimate needs to be obtained, default is sample median
  #' @param normalize should we divide by n?
  
  # get the dimensions of the data and some initial parameters
  ntime <- nrow(Y)
  nloc <- ncol(Y)
  alpha <- log(nloc) / log(ntime)  
  beta <- 1.5 * alpha - 1
  nbn <- ntime^(-beta)
  eps <- 1e-10
  
  if (is.null(Z)) {
    distmat <- NULL
  } else {
    distmat <- as.matrix(dist(Z))
  }
  kern.array = lapply(c(1:ntime), FUN = function(j) kern(Xt[j,], x, nloc, distmat))
  kern.norm <- unlist(lapply(kern.array, FUN = function(x1) norm(x1)))
  best.index <- which.max(kern.norm)
  idx = which((abs(1:ntime - best.index) < nbn) & (kern.norm > eps))  # apply the kernel to the support specified by the bandwidth
  
  # obtaining the eta(t_i) values as in Eq. (12)
  muhat <- t(sapply(1:length(idx), FUN = function(i) {
    estimate.quantile(Y, Xt, Z, 0.5, Xt[idx[i], ])[["q"]]
  }))  # transposing ensures each row corresponds to a timepoint
  eta <- t(sapply(1:length(idx), function(i) (kern.array[[idx[i]]] %*% (muhat[i, ] - q0))))
  
  # creating V_2(t_i) values as in Eq. (14)
  v2 <- lapply(1:length(idx), FUN = function(i) {
    etanorm <- inner.norm(eta[i, ])  # norm value || eta(t_i) ||
    if (etanorm > eps) {
      return((diag(nloc) - outer(eta[i, ], eta[i, ])/ etanorm^2)/etanorm)
    } else {
      return((diag(nloc) -  matrix(1/nloc, nrow = nloc, ncol = nloc)))  # handle very small eta values
    }
  })
  term.calc <- lapply(1:length(idx),FUN = function(i) (kern.array[[idx[i]]] %*% v2[[i]] %*% kern.array[[idx[i]]] ))
  if (normalize) {
    return(Reduce('+', term.calc) / (ntime))
  } else {
    return(Reduce('+', term.calc))
  }
}


Omega_t <- function(q0, Y, Xt, Z = NULL, x = colMedians(Xt), normalize = TRUE) {
  #' @param q0 the quantile vector at which Psi_t(q_0) is to be calculated 
  #' @param Y data-matrix of the response variable
  #' @param Xt data-matrix of the covariate process 
  #' @param Z latitude-longitude information (if available)
  #' @param x value at which the estimate needs to be obtained, default is sample median
  #' @param normalize should we divide by n^2p^2b_n^2?
  
  # get the dimensions of the data and some initial parameters
  ntime <- nrow(Y)
  nloc <- ncol(Y)
  alpha <- log(nloc) / log(ntime)  
  beta <- 1.5 * alpha - 1
  nbn <- ntime^(-beta)
  eps <- 1e-10
  
  if (is.null(Z)) {
    distmat <- NULL
  } else {
    distmat <- as.matrix(dist(Z))
  }
  kern.array = lapply(c(1:ntime), FUN = function(j) kern(Xt[j,], x, nloc, distmat))
  kern.norm <- unlist(lapply(kern.array, FUN = function(x1) norm(x1)))
  best.index <- which.max(kern.norm)
  idx = which((abs(1:ntime - best.index) < nbn) & (kern.norm > eps))  # apply the kernel to the support specified by the bandwidth
  
  # obtaining the eta(t_i) values as in Eq. (12)
  muhat <- t(sapply(1:length(idx), FUN = function(i) {
    estimate.quantile(Y, Xt, Z, 0.5, Xt[idx[i], ])[["q"]]
  }))  # transposing ensures each row corresponds to a timepoint
  eta <- t(sapply(1:length(idx), function(i) (kern.array[[idx[i]]] %*% (muhat[i, ] - q0))))
  
  # creating V_2(t_i) values as in Eq. (14)
  v2 <- lapply(1:length(idx), FUN = function(i) {
    etanorm <- inner.norm(eta[i, ])  # norm value || eta(t_i) ||
    if (etanorm > eps) {
      return((diag(nloc) - outer(eta[i, ], eta[i, ])/ etanorm^2)/etanorm)
    } else {
      return((diag(nloc) -  matrix(1/nloc, nrow = nloc, ncol = nloc)))  # handle very small eta values
    }
  })
  KV2K <- lapply(1:length(idx),FUN = function(i) (kern.array[[idx[i]]] %*% v2[[i]] %*% kern.array[[idx[i]]] ))
  
  # estimate the sigma
  sigmas <- lapply(1:length(idx), FUN = function(i) (estimate.sigma(Y, Xt, Z, Xt[idx[i],])) )
  
  # calculate the term K V2 V1 V2 K
  term.calc <- lapply(1:length(idx), FUN = function(i) (KV2K[[i]] %*% sigmas[[i]] %*% KV2K[[i]]) )
  if (normalize) {
    return(Reduce('+', term.calc) / (nloc^2 * nbn^2))
  } else {
    return(Reduce('+', term.calc))
  }
}

