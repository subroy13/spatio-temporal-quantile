# load other functions
source('scripts/helper.R')
source('scripts/SD_est_functions.R')


# H0: Some restrictions on Q

# function to test any general hypothesis
test.quantile <- function(Y, estimate.H0, Xt = NULL, Z = NULL, quant = 0.5, est.quant = NULL,sig.level = 0.05, verbose = TRUE) {
  #' @param Y data-matrix of the response variable
  #' @param estimate.H0 a function that takes input (Y, Xt, current_quantile, Z, quant) => quantile_under_H0 by applying the restrictions
  #' @param Xt data-matrix of the covariate variable (if not available, takes the timepoints)
  #' @param Z latitude-longitude information (if available)
  #' @param quant which quantile to work on, default is median (0.5)
  #' @param est.quant estimated quantile matrix if it has been already computed
  #' @param sig.level the significance level at which the hypothesis is to be rejected
  #' @param verbose whether to print progress bar or not
  
  # initialize some parameters
  ntime <- nrow(Y)
  nloc <- ncol(Y)
  alpha <- log(nloc) / log(ntime)  
  beta <- 1.5 * alpha - 1
  bn <- ntime^(-1-beta)
  eps <- 1e-4
  
  if (is.null(Xt)) {
    Xt <- matrix(c(1:ntime)/ntime, ncol = 1)  # only one covariate is time
  }
  
  # find the estimator at each time point
  if (is.null(est.quant)){
    quantiles <- matrix(NA, nrow = ntime, ncol = nloc)  # to hold the quantiles
    for (i in 1:ntime) {
      quantiles[i, ] <- estimate.quantile(Y, Xt, Z, quant, x = Xt[i,])[["q"]]
    }
  } else {
    quantiles <- est.quant
  }
  
  # estimate of quantiles under H0
  q0s <- estimate.H0(Y, Xt, quantiles, Z, quant)  # returns (ntime, nloc) matrix
  
  # calculate the test statistic by computing required eta, psi and omega
  testmat.sig = numeric(ntime)
  
  # # running parallel computation for saving time
  # library(doParallel)
  # library(foreach)
  # cluster <- makeCluster((detectCores() - 2)) # leaving two cores for other work in the machine 
  # registerDoParallel(cluster)
  
  # computation for the test statistic
  testmat.sig = numeric(ntime) 
  for(i in 1:ntime){
    
    # ignore normalization factors for now, leads to better numerical stability
    # compute the eta_t, psi_t, omega_t as in Eq. (15), (16) and (17)
    out = calc.for.tests(q0 = quantiles[i, ], Y, Xt, Z, x = Xt[i,], normalize = FALSE)
    eta_ts <- out$eta
    psi_ts <- out$Psi
    omega_ts <- out$Omega
    
    # the mu_q vector, mean of gaussian process
    mu_q <- solve(psi_ts + eps*diag(nloc), eta_ts)
    
    # sigma_q^{-1/2}, where sigma_q is the covariance of the gaussian process
    sigma_q.invsq <- tcrossprod(
      matsq(psi_ts) %*% matsq(omega_ts, inverse = TRUE), 
      matsq(psi_ts)
    )
    
    # create the test matrix for i-th timepoint
    testmat <- tcrossprod(sigma_q.invsq %*% (quantiles[i,] - q0s[i,] - mu_q / (nloc * bn)) )  
    
    # largest singular value of the test matrix
    testmat.sig[i] = norm(testmat, type = "2") 
  }

  # final adjustment to test statistic for normalization for the ignored normalization
  teststat <- max(testmat.sig) * (nloc * bn)  
  
  # calculate the p-value
  m_n <- exp((nloc * bn)^(-2))
  gumbel_crit <- inv.folded.crit(m_n, teststat + sqrt(2 * log(m_n)))
  pval <- 1 - exp(-exp(-gumbel_crit))
  
  return(list(
    "statistic" = teststat,
    "pvalue" = pval,
    "reject" = (pval < sig.level)
  ))
}



###########################
# H01: Q(tau, x) = beta(x) (1, 1, ..., 1)  (covariate effect is same across all spatial points)

# function to test the above hypothesis
test.type1 <- function(Y, Xt, Z = NULL, quant = 0.5, est.quant = NULL, sig.level = 0.05, verbose = TRUE) {
  #' @param Y data-matrix of the response variable
  #' @param Xt data-matrix of the covariate variable
  #' @param Z latitude-longitude information (if available)
  #' @param quant which quantile to work on, default is median (0.5)
  #' @param est.quant estimated quantile matrix if it has been already computed
  #' @param sig.level the significance level at which the hypothesis is to be rejected
  #' @param verbose whether to print progress bar or not
  
  estimate.H0 <- function(Y, Xt, current_quantiles, Z = NULL, quant = 0.5) {
    # under H0, each row of current_quantiles should belong to span of (1, 1, ... 1)
    ntime <- nrow(current_quantiles)
    nloc <- ncol(current_quantiles)
    row_mean <- rowMeans(current_quantiles)
    return(matrix(row_mean, nrow = ntime, ncol = nloc, byrow = FALSE))
  }
  
  return(test.quantile(Y, estimate.H0, Xt, Z, quant, est.quant, sig.level, verbose))
}


###########################
# H02: Q(tau, t) = tQ(tau, 1), for a single value of tau, but for all t

# function to test the above hypothesis
test.type2 <- function(Y, Z = NULL, quant = 0.5, est.quant = NULL, sig.level = 0.05, verbose = TRUE) {
  #' @param Y data-matrix of the response variable
  #' @param Z latitude-longitude information (if available)
  #' @param quant which quantile to work on, default is median (0.5)
  #' @param est.quant estimated quantile matrix if it has been already computed
  #' @param sig.level the significance level at which the hypothesis is to be rejected
  #' @param verbose whether to print progress bar or not
  
  estimate.H0 <- function(Y, Xt, current_quantiles, Z = NULL, quant = 0.5) {
    # under H0, an estimate of Q(tau, 1) = mean_t (Q(tau, t) / t)
    ntime <- nrow(current_quantiles)
    nloc <- ncol(current_quantiles)
    q1 <- colMeans(quantiles / Xt[,1])
    return(tcrossprod(Xt[,1], q1))
  }
  
  return(test.quantile(Y, estimate.H0, NULL, Z, quant, est.quant, sig.level, verbose))
}