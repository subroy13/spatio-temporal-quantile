# required packages
library(matrixStats)

# Load the estimator function
source('./scripts/estimator.R')

###########################
# H01: Q(tau, x) = beta(x) Qval(tau)




###########################
# H02: Q(tau, t) = tQ(tau, 1), for a single value of tau, but for all t

# function to test the above hypothesis
test.type2 <- function(Y, Z = NULL, quant = 0.5, sig.level = 0.05, verbose = TRUE) {
  #' @param Y data-matrix of the response variable
  #' @param Z latitude-longitude information (if available)
  #' @param quant which quantile to work on, default is median (0.5)
  #' @param sig.level the significance level at which the hypothesis is to be rejected
  #' @param verbose whether to print progress bar or not
  
  # initialize some parameters
  ntime <- nrow(Y)
  nloc <- ncol(Y)
  alpha <- log(nloc) / log(ntime)  
  beta <- 1.5 * alpha - 1
  bn <- ntime^(-1-beta)
  eps <- 1e-10

  Xt <- matrix(c(1:ntime)/ntime, ncol = 1)  # only one covariate is time

  # find the estimator at each time point
  quantiles <- matrix(NA, nrow = ntime, ncol = nloc)  # to hold the quantiles
  eta_ts <- as.list(1:ntime)
  psi_ts <- as.list(1:ntime)
  omega_ts <- as.list(1:ntime)

  pb <- txtProgressBar(min=0, max = 2 * ntime, style = 3)
  for (i in 1:ntime) {
    quantiles[i, ] <- estimate.quantile(Y, Xt, Z, quant, x = Xt[i,])[["q"]]

    # ignore normalization factors for now, leads to better numerical stability
    eta_ts[i, ] <- eta_t(quantiles[i, ], Y, Xt, Z, x = Xt[i,], normalize = FALSE)
    psi_ts[i, ] <- Psi_t(quantiles[i, ], Y, Xt, Z, x = Xt[i,], normalize = FALSE)
    omega_ts[i, ] <- Omega_t(quantiles[i, ], Y, Xt, Z, x = Xt[i, ], normalize = FALSE)
    setTxtProgressBar(pb, value = i)
  }

  # estimate of quantile under H0 will be mean(Q(tau, t)/t)
  q1 <- colMeans(quantiles / Xt[,1])

  teststat <- 0  # initial value of test statistic
  for (i in 1:ntime) {
    mu_q <- solve(psi_t[i,] + eps*diag(nloc), eta_ts[i, ])  # the mu_q vector, mean of gaussian process
    sigma_q.invsq <- tcrossprod(matsq(psi_ts[i,]) %*% matsq(omega_ts[i,], inverse = TRUE), matsq(psi_ts[i, ]))
    testmat <- tcrossprod(sigma_q.invsq %*% (quantiles[i, ] - q1 - mu_q / (nloc * bn)) )  # create the test matrix for t-th timepoint
    testmat.sig <- norm(testmat, type = "2")  # largest singular value of the test matrix
    teststat <- max(testmat.sig, teststat)   # update the test statistic by taking supremum 
    setTxtProgressBar(pb, value = i + ntime)
  }
  close(pb)

  teststat <- teststat * (nloc * bn)  # final adjustment to test statistic for normalization
  m_n <- exp((nloc * bn)^(-2))

  # calculate the p-value
  gumbel_crit <- inv.folded.crit(m_n, teststat + sqrt(2 * log(m_n)))
  pval <- 1 - exp(-exp(-gumbel_crit))

  return(list(
    "statistic" = teststat,
    "pvalue" = pval,
    "reject" = (pval < sig.level)
  ))
}
