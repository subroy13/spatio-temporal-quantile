# Estimating The Parameter
est_nr <- function(dat, x, tau = 0.75, maxiter = 1000, tol = 1e-5, rho = 0.5, bandwidth = 1e-2) {
    n <- nrow(dat$Y)
    p <- ncol(dat$Y)
    if (is.null(dim(dat$X))) {
        dim(dat$X) <- c(n, 1)    # convert to column vector
    }
    
    # Creating kernel
    if (is.null(dat$Z)) {
        d <- abs(outer(1:p, 1:p, FUN = "-"))
    } else {
        d <- as.matrix(dist(dat$Z))^2
    }
    kernmat <- function(t1, t2, bb = bandwidth) {
        tdiff <- sqrt(sum((t1 - t2)^2) )
        return(exp(-bb * (tdiff + d) ))
    } 
    
    # Demean Y
    # Y.mean <- colMeans(dat$Y)
    Y.mean <- rep(0, p)
    Y.demean <- dat$Y - matrix(Y.mean, nrow = n, ncol = p, byrow = TRUE)
    
    niter <- 0
    converged <- FALSE
    err <- Inf
    q <- runif(p)
    u <- (2 * tau - 1) * rep(1, p) / sqrt(p)
    
    kernmatList <- lapply(1:n, FUN = function(i) {
        return(kernmat(dat$X[i, ], x))
    })
    
    kernmatList.sq <- lapply(1:n, FUN = function(i) {
        return(kernmatList[[i]] %*% kernmatList[[i]])
    })
    
    kmatsum <- Reduce('+', kernmatList) %*% u
    
    while (!converged) {
        # Find next best q
        
        norms <- numeric(n)
        for (i in 1:n) {
            e <- (dat$Y[i, ] - q)
            norms[i] <- sqrt(t(e) %*% kernmatList.sq[[i]] %*% e)
        }
        denom <- Reduce('+', x = lapply(1:n, FUN = function(i) {
            return(2 * (1/norms[i]) * kernmatList.sq[[i]])
        }))
        denom <- denom + rho * diag(p)
        
        num <- Reduce('+', x = lapply(1:n, FUN = function(i) {
            2 * (1/norms[i]) * kernmatList.sq[[i]] %*% Y.demean[i, ]
        }))
        num <- num + kmatsum
        
        qnew <- MASS::ginv(denom) %*% num
        
        # Error calculation
        err <- sqrt(sum((q - qnew)^2)) / sqrt(sum(q^2))
        q <- qnew
        niter <- niter + 1
        
        # Convergence check
        if ((niter > maxiter) | (err < tol)) {
            converged <- TRUE
        }  
    }
    return (
        list(
            "q" = q + Y.mean, "niter" = niter, "err" = err
        )
    )
}

