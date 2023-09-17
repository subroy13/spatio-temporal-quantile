############################################
# Our Model Simulation
# 
# One Idea: X ~ F, then F(X) ~ U, 
# For quantile, we know Q(F(X)) = X, so, Q(U) = X ~ F, where U ~ Unif(0, 1)
############################################

# Finds the quantile for the normal 
library(mvtnorm)
library(mlVAR)

true.quant.norm <- function(p, tau = 0.5) {
    B <- 10000
    qq <- seq(-3, 3, 0.01)
    qn <- length(qq)
    
    normexp <- numeric(qn)
    for (i in 1:qn) {
        u <- (2 * qq[i] - 1) * rep(1, p)
        X <- mvtnorm::rmvnorm(B, mean = -u, sigma = diag(p))
        normexp[i] <- mean(sqrt(rowSums(X^2)) + rowSums(X) * (2 * tau - 1) / sqrt(p) )
    }
    plot(qq, normexp, type = "l")
    return(qq[which.min(normexp)])
}

#true.quant.norm(10, tau = 0.95)

# CREATE A MAPPING
normquant <- c(
    "0.5" = 0.48,
    "0.75" = 0.79,
    "0.8" = 0.88,
    "0.85" = 0.97,
    "0.9" = 1.12,
    "0.95" = 1.55,
    "0.99" = 2.85
)


sim.model1 <- function(n, p, seed = 1234,
                      quants = c(0.5, 0.75, 0.8, 0.85, 0.9, 0.95)) {
    set.seed(seed)
    # n = # of timepoints
    # p = # of locations
    
    params.VAR <- matrix(c(0.2, 0.1, -0.3, 0.4), nrow = 2, ncol = 2)
    params.slope <- matrix(runif( 3*p, min = -1, max = 1), nrow = 3, ncol = p)
    params.slope[-1, ] <- params.slope[-1, ] * 20
    params.intercept <- rnorm(p)
    Xt <- as.matrix(cbind(1:n, simulateVAR(params.VAR, Nt = n )))
    
    true_quants <- list()
    for (i in 1:length(quants)) {
        true_quants[[i]] <- matrix(0, n, p)
    }
    
    for (t in 1:n) {
        mu <- params.intercept + Xt[t, 1] * params.slope[1, ] + sin(2 * pi * Xt[t, 2]) * params.slope[2] + cos(2 * pi * Xt[t, 3]) * params.slope[3]
        # reduce error variance so that the signal is prominent
        Sigma <- exp(-0.1 * outer(1:nloc, 1:nloc, FUN = "-")^2) * mean(abs(mu)) * 0.1
        Y[t, ] <- MASS::mvrnorm(mu = mu, Sigma = Sigma)
        for (i in 1:length(quants)) {
            qq <- quants[[i]]
            true_quants[[i]][t, ] <- mu + chol(Sigma) %*% normquant[[as.character(quants)]]
        }
    }
    names(true_quants) <- as.character(quants)
    
    return(list(
        "Y" = Y,
        "X" = Xt,
        "Z" = NULL,
        "true_quants" = true_quants
    ))
}

#############################################
# Perform the simulation
source('./scripts/estimator.R')
library(pbapply)

quants <- c(0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99)
ts <- c(0, 0.25, 0.5, 0.75, 1)

calc_error <- function(true_quants, result) {
    errs <- matrix(NA, nrow = length(quants), ncol = length(ts))
    rmses <- matrix(NA, nrow = length(quants), ncol = length(ts))
    maes <- matrix(NA, nrow = length(quants), ncol = length(ts))
    
    for (i in 1:length(result)) {
        tvals <- true_quants[[1]][1 + ts*100, ]
        esvals <- result[[i]]
        for (j in 1:nrow(esvals)) {
            rmses[i, j] <- sqrt(sum((tvals[j, ] - esvals[j, ])^2))
            maes[i, j] <- sum(abs(tvals[j, ] - esvals[j, ])^2)
            errs[i, j] <- maes[i, j] / sum(abs(tvals[j, ]))
        }
    }
    return(list(
        "RMSE" = rmses,
        "MAE" = maes,
        "PError" = errs
    ))
}


simulation.inference <- function(B, seed = 1234, rh = 0.2, bh = 0.1) {
    set.seed(seed)
    seeds <- sample(1e8, size = B, replace = FALSE)   # seeds for each MC sample
    errs <- matrix(0, nrow = length(quants), ncol = length(ts))
    rmses <- matrix(0, nrow = length(quants), ncol = length(ts))
    maes <- matrix(0, nrow = length(quants), ncol = length(ts))
    
    pb <- txtProgressBar(max = B, style = 3)
    datas <- list()
    for (b in 1:B) {
        datas[[b]] <- sim.model1(200, 10, Z, seed = seeds[b])
    }
    
    for (b in 1:B) {
        # generate the data
        dat <- datas[[b]]
        
        # do the estimate for each x and tau
        results <- list()
        for (i in 1:length(quants)) {
            qq <- quants[[i]]
            results[[i]] <- t(sapply(ts, FUN = function(x) {
                return(est_nr(dat, x, tau = qq, rho = rh, bandwidth = bh)$q)
            }))
        }
        
        # calculate errors
        tmp <- calc_error(dat$true_quants, results)
        rmses <- rmses + tmp$RMSE^2
        maes <- maes + tmp$MAE
        errs <- errs + tmp$PError
        
        setTxtProgressBar(pb, value = b)
    }
    close(pb)
    
    rmses <- sqrt(rmses / B)
    errs <- errs / B
    maes <- maes / B
    
    return(list(
        "RMSE" = rmses,
        "MAE" = maes,
        "PError" = errs
    ))
}


output <- simulation.inference(B = 10, seed = 2022, rh = 0.2, bh = 1e-1)


xtable::xtable(output$MSE)


###############################
# Prediction Performance Checking
ts <- seq(0.91, 1, by = 0.01)
simulation.prediction <- function(B, seed = 1234, rh = 0.2, bh = 0.1) {
    set.seed(seed)
    seeds <- sample(1e8, size = B, replace = FALSE)   # seeds for each MC sample
    errs <- matrix(0, nrow = length(quants), ncol = length(ts))
    mses <- matrix(0, nrow = length(quants), ncol = length(ts))
    
    pb <- txtProgressBar(max = B, style = 3)
    Z <- matrix(runif(10 * 2, min = 0, max = 1), nrow = 10, ncol = 2)  # locations
    
    datas <- list()
    for (b in 1:B) {
        datas[[b]] <- sim.bayes(101, 10, Z, seed = seeds[b])
    }
    
    for (b in 1:B) {
        # generate the data
        dat <- datas[[b]]
        dat <- dat$Y[1:90,]
        
        # do the estimate for each x and tau
        results <- list()
        for (i in 1:length(quants)) {
            qq <- quants[[i]]
            results[[i]] <- t(sapply(ts, FUN = function(x) {
                return(est_nr(dat, x, tau = qq, rho = rh, bandwidth = bh)$q)
            }))
        }
        
        # calculate errors
        tmp <- calc_error(dat$true_quants, results)
        mses <- mses + tmp$MSE
        errs <- errs + tmp$PError
        
        setTxtProgressBar(pb, value = b)
    }
    close(pb)
    
    mses <- sqrt(mses / B)
    errs <- errs / B
    
    return(list(
        "MSE" = mses,
        "PError" = errs
    ))
}






