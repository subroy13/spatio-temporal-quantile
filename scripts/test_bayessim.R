############################################
# Bayesian Spatio-temporal Model Simulation
# Reference:  P. Das, S. Ghoshal paper
# 
# One Idea: X ~ F, then F(X) ~ U, 
# For quantile, we know Q(F(X)) = X, so, Q(U) = X ~ F, where U ~ Unif(0, 1)
############################################

sim.bayes <- function(n, p, Z, seed = 1234,
                      quants = c(0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99)) {
    set.seed(seed)
    # n = # of timepoints
    # p = # of locations
    
    X <- seq(0, 1, length.out = n)
    
    U <- matrix(runif(n * p), nrow = n, ncol = p)
    Y <- matrix(NA, nrow = n, ncol = p)
    
    for (t in 1:n) {
        xi1 <- (1 - (Z[, 1] + Z[, 2])/2 ) * U[t, ]^2 + Z[, 1] * log(1+U[t, ])/(2*log(2)) + Z[, 2]/2 * U[t, ]^3
        xi2 <- (1 - Z[, 2]^2) * sin(pi * U[t, ]/2) + Z[, 2]^2 * (exp(U[t, ]) - 1)/(exp(1) - 1)
        Y[t, ] <- X[t] * xi1 + (1 - X[t]) * xi2
    }
    
    true_quants <- list()
    for (i in 1:length(quants)) {
        qq <- quants[[i]]
        xi1 <- (1 - (Z[, 1] + Z[, 2])/2 ) * qq^2 + Z[, 1] * log(1+qq)/(2*log(2)) + Z[, 2]/2 * qq^3
        xi2 <- (1 - Z[, 2]^2) * sin(pi * qq/2) + Z[, 2]^2 * (exp(qq) - 1)/(exp(1) - 1)
        true_quants[[i]] <- outer(X, xi1) + outer(1-X, xi2)
    }
    names(true_quants) <- as.character(quants)
    
    return(list(
        "Y" = Y,
        "X" = X,
        "Z" = Z,
        "true_quants" = true_quants
    ))
}

#Z <- matrix(runif(10 * 2, min = 0, max = 1), nrow = 10, ncol = 2)  # locations
#dat <- sim.bayes(200, 10, Z, seed = 1)

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
        tvals <- true_quants[[as.character(quants[[i]])]][1 + ts*100, ]
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
    Z <- matrix(runif(15 * 2, min = 0, max = 1), nrow = 15, ncol = 2)  # locations
    
    datas <- list()
    for (b in 1:B) {
        datas[[b]] <- sim.bayes(101, 15, Z, seed = seeds[b])
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


output <- simulation.inference(B = 100, seed = 2022, rh = 0.2, bh = 1e-1)

format.output <- function(output) {
    mat <- matrix("", nrow = nrow(output$MAE), ncol = ncol(output$MAE))
    for (i in 1:nrow(output$MAE)) {
        for (j in 1:ncol(output$PError)) {
            mat[i,j] <- paste0(
                round(output$MAE[i,j],3), 
                " (", 
                round(output$PError[i,j]*100,2), 
                ")" 
            )
        }
    }
    return(xtable::xtable(mat))
}

format.output(output)
#xtable::xtable(output$MAE, digits = 3)
#xtable::xtable(output$PError, digits = 3)


###############################
# Prediction Performance Checking
ts <- seq(0.91, 1, by = 0.01)
simulation.prediction <- function(B, seed = 1234, rh = 0.2, bh = 0.1) {
    set.seed(seed)
    seeds <- sample(1e8, size = B, replace = FALSE)   # seeds for each MC sample
    errs <- matrix(0, nrow = length(quants), ncol = length(ts))
    rmses <- matrix(0, nrow = length(quants), ncol = length(ts))
    maes <- matrix(0, nrow = length(quants), ncol = length(ts))
    
    
    pb <- txtProgressBar(max = B, style = 3)
    Z <- matrix(runif(15 * 2, min = 0, max = 1), nrow = 15, ncol = 2)  # locations
    
    datas <- list()
    for (b in 1:B) {
        datas[[b]] <- sim.bayes(101, 15, Z, seed = seeds[b])
    }
    
    for (b in 1:B) {
        # generate the data
        dat <- datas[[b]]
        dat$Y <- dat$Y[1:90,1:10]
        dat$Z <- dat$Z[1:10,]
        dat$X <- dat$X[1:90]
        
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


output <- simulation.prediction(B = 5, seed = 2022, rh = 0.2, bh = 1e-1)
format.output(output)




























