# required packages
library(matrixStats)

# Load the estimator function
source('./scripts/estimator.R')

##########################
# H0: Some restrictions on Q

# function to test any general hypothesis
test.quantile <- function(Y, estimate.H0, Xt = NULL, Z = NULL, quant = 0.5, sig.level = 0.05, verbose = TRUE) {
  #' @param Y data-matrix of the response variable
  #' @param estimate.H0 a function that takes input (Y, Xt, current_quantile, Z, quant) => quantile_under_H0 by applying the restrictions
  #' @param Xt data-matrix of the covariate variable (if not available, takes the timepoints)
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
  
  if (is.null(Xt)) {
    Xt <- matrix(c(1:ntime)/ntime, ncol = 1)  # only one covariate is time
  }
  
  # find the estimator at each time point
  quantiles <- matrix(NA, nrow = ntime, ncol = nloc)  # to hold the quantiles
  eta_ts <- as.list(1:ntime)
  psi_ts <- as.list(1:ntime)
  omega_ts <- as.list(1:ntime)
  
  pb <- txtProgressBar(min=0, max = 2 * ntime, style = 3)
  for (i in 1:ntime) {
    quantiles[i, ] <- estimate.quantile(Y, Xt, Z, quant, x = Xt[i,])[["q"]]
    
    # ignore normalization factors for now, leads to better numerical stability
    # compute the eta_t, psi_t, omega_t as in Eq. (15), (16) and (17)
    eta_ts[[i]] <- eta_t(quantiles[i, ], Y, Xt, Z, x = Xt[i,], normalize = FALSE)
    psi_ts[[i]] <- Psi_t(quantiles[i, ], Y, Xt, Z, x = Xt[i,], normalize = FALSE)
    omega_ts[[i]] <- Omega_t(quantiles[i, ], Y, Xt, Z, x = Xt[i, ], normalize = FALSE)
    setTxtProgressBar(pb, value = i)
  }
  
  # estimate of quantiles under H0
  q0s <- estimate.H0(Y, Xt, quantiles, Z, quant)  # returns (ntime, nloc) matrix
  
  teststat <- 0  # initial value of test statistic
  for (i in 1:ntime) {
    # the mu_q vector, mean of gaussian process
    mu_q <- solve(psi_t[[i]] + eps*diag(nloc), eta_ts[[i]])
    
    # sigma_q^{-1/2}, where sigma_q is the covariance of the gaussian process
    sigma_q.invsq <- tcrossprod(
      matsq(psi_ts[[i]]) %*% matsq(omega_ts[[i]], inverse = TRUE), 
      matsq(psi_ts[[i]])
    )
    
    # create the test matrix for i-th timepoint
    testmat <- tcrossprod(sigma_q.invsq %*% (quantiles[i,] - q0s[i,] - mu_q / (nloc * bn)) )  
    
    # largest singular value of the test matrix
    testmat.sig <- norm(testmat, type = "2")  
    
    # update the test statistic by taking supremum 
    teststat <- max(testmat.sig, teststat)
    setTxtProgressBar(pb, value = i + ntime)
  }
  close(pb)
  
  # final adjustment to test statistic for normalization for the ignored normalization
  teststat <- teststat * (nloc * bn)  
  
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
test.type1 <- function(Y, Xt, Z = NULL, quant = 0.5, sig.level = 0.05, verbose = TRUE) {
  #' @param Y data-matrix of the response variable
  #' @param Xt data-matrix of the covariate variable
  #' @param Z latitude-longitude information (if available)
  #' @param quant which quantile to work on, default is median (0.5)
  #' @param sig.level the significance level at which the hypothesis is to be rejected
  #' @param verbose whether to print progress bar or not
  
  estimate.H0 <- function(Y, Xt, current_quantiles, Z = NULL, quant = 0.5) {
    # under H0, each row of current_quantiles should belong to span of (1, 1, ... 1)
    ntime <- nrow(current_quantiles)
    nloc <- ncol(current_quantiles)
    row_mean <- rowMeans(current_quantiles)
    return(matrix(row_mean, nrow = ntime, ncol = nloc, byrow = FALSE))
  }
  
  return(test.quantile(Y, estimate.H0, Xt, Z, quant, sig.level, verbose))
}


###########################
# H02: Q(tau, t) = tQ(tau, 1), for a single value of tau, but for all t

# function to test the above hypothesis
test.type2 <- function(Y, Z = NULL, quant = 0.5, sig.level = 0.05, verbose = TRUE) {
  #' @param Y data-matrix of the response variable
  #' @param Z latitude-longitude information (if available)
  #' @param quant which quantile to work on, default is median (0.5)
  #' @param sig.level the significance level at which the hypothesis is to be rejected
  #' @param verbose whether to print progress bar or not
  
  estimate.H0 <- function(Y, Xt, current_quantiles, Z = NULL, quant = 0.5) {
    # under H0, an estimate of Q(tau, 1) = mean_t (Q(tau, t) / t)
    ntime <- nrow(current_quantiles)
    nloc <- ncol(current_quantiles)
    q1 <- colMeans(quantiles / Xt[,1])
    return(tcrossprod(Xt[,1], q1))
  }
  
  return(test.quantile(Y, estimate.H0, NULL, Z, quant, sig.level, verbose))
}


############
# Now apply for the dataset

library(readr)
library(dplyr)
library(tidyr)

df <- read_csv('./data/EPM Data/my_epm1_processed.csv')
colSums(table(df$Household, df$Cluster) > 1)

clustNum <- 12   # change the cluster number here
quant_level <- 0.75  # change the quantile level here

Xfull <- df %>% 
  filter(Cluster == clustNum) %>%
  select(WeekT, DayN, HH, DayW, hh, hw, wh, ww) %>% distinct() %>%
  as.matrix()
dim(Xfull)

Yfull <- df %>% 
  filter(Cluster == clustNum) %>%
  pivot_wider(names_from = "Household", values_from = "Demand") %>%
  select(-c(DateTime, WeekT, DayN, HH, DayW, Season, DayC, hh, hw, wh, ww, Cluster, Train)) %>%
  as.matrix()
dim(Yfull)

# do the hypothesis testing at specified quantile level
test.type1(Yfull, Xfull, quant = quant_level, verbose = TRUE)
test.type2(Yfull, quant = quant_level, verbose = TRUE)
