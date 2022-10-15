########################
# Inputs:
#   - The Y matrix (n timepoints x p locations)
#   - The X matrix (n timepoints x d dimension)
#   - The quantile to obtain the estimate for (tau)
#   - The kernel function to use
########################
library(matrixStats)

fit.quantile <- function(Y, X, xnew = colMedians(X), tau = 0.5, locations = NULL,
                         maxiter = 100, eps = 1e-5, delta = 1e-8, verbose = FALSE, bn = 0.5) {
    n <- nrow(Y)
    p <- ncol(Y)
    d <- ncol(X)
    if (is.null(locations)) {
        dismat <- abs(outer(1:p, 1:p, FUN = "-"))
    } else {
        dismat <- as.matrix(dist(locations))
    }
    u <- (2 * tau - 1) * rep(1, p)/sqrt(p)
    
    err <- (X - matrix(xnew, nrow = n, ncol = d, byrow = T))
    normdiffs <- sqrt(rowSums(err^2))
    
    init.index <- which.min(normdiffs)
    init.q <- Y[init.index, ] + chol(var(Y)) %*% u
    q <- init.q
    
    for (iter in 1:maxiter) {
        weights <- numeric(n)
        num <- numeric(p)
        denom <- matrix(0, nrow = p, ncol = p)
        
        for (i in 1:n) {
            err1 <- Y[i, ] - q
            kk <- exp(-bn * dismat * normdiffs[i])
            weights[i] <- crossprod(kk %*% err1)
            weights[i] <- 1 / sqrt(weights[i] + delta)
            num <- num + (2 * weights[i] * crossprod(kk) %*% Y[i, ]) + (kk %*% u)
            denom <- denom + (2 * weights[i] * crossprod(kk))
        }
        
        qnew <- solve(denom + delta*diag(p), num)
        relerr <- sqrt(sum((qnew - q)^2) / sum(q^2))
        if (relerr < eps) {
            break
        } else {
            q <- qnew
            if (verbose) {
                message(paste0("Iteration ", iter, "\tCurrent error: ", round(relerr, 5), "\tCurrent Value: ", paste(round(qnew,3), collapse = ",") ))
            }
        }
    }

    return(list(
        "q" = q
    ))    
}

############################################
# Bayesian Spatio-temporal Model Simulation
# Reference:  P. Das, S. Ghoshal paper
# 
# One Idea: X ~ F, then F(X) ~ U, 
# For quantile, we know Q(F(X)) = X, so, Q(U) = X ~ F, where U ~ Unif(0, 1)
############################################

sim.bayes <- function(n, p, Z, seed = 1234,
                      quants = c(0.5, 0.75, 0.8, 0.85, 0.9, 0.95)) {
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
    dim(X) <- c(n, 1)
    
    return(list(
        "Y" = Y,
        "X" = X,
        "Z" = Z,
        "true_quants" = true_quants
    ))
}


# Run simulation for this setup (Inference)
n <- 200
p <- 15
set.seed(2022)
Z <- matrix(runif(p * 2, min = 0, max = 1), nrow = p, ncol = 2)  # locations
dat <- sim.bayes(n, p, Z, seed = 2022)


fit.quantile(dat$Y, dat$X, xnew = dat$X[1, ], tau = 0.75, locations = dat$Z, 
             verbose = F, bn = 5)$q












