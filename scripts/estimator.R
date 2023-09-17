####################
# Code for the estimator

library(matrixStats)

inner.norm <- function(x) { sqrt(sum(x^2)) }

inner.prod <- function(x, y) { sum(x*y) }

# True quantile
true.quantile <- function(Ydata, Xt, alpha = 0.5, x = colMedians(Xt), method = "BFGS", ...) {
    nloc <- ncol(Ydata[[1]])
    ntime <- nrow(Ydata[[1]])
    index <- which.min(rowSums((Xt - matrix(x, nrow = ntime, ncol = ncol(Xt), byrow = T) )^2))
    
    loss <- function(q, u) {
        mean(sapply(1:length(Ydata), FUN = function(i) {
            inner.norm(Ydata[[i]][index, ] - q) + inner.prod(u, Ydata[[i]][index, ] - q)            
        }))
    }
    
    init.q <- rowMedians(sapply(1:length(Ydata), FUN = function(i){ Ydata[[i]][index, ] }))
    init.u <- rep(1,nloc)/inner.norm(rep(1, nloc)) * (2 * alpha - 1)
    return(optim(init.q, fn = loss, u = init.u, method = method, ...)$par)
}


# Estimator of the quantile
estimate.quantile <- function(Y, Xt, Z = NULL, quant = 0.5, x = colMedians(Xt), maxiter = 1, tol = 1e-5) {
    
    ntime <- nrow(Y)
    nloc <- ncol(Y)
    
    alpha <- log(nloc) / log(ntime)
    if (alpha > 1/2) {
        warning("The number of locations are higher than expected, theoretical gurantees may not be accurate.")
    }
    beta <- 1.5 * alpha - 1
    nbn <- ntime^(-beta)
    
    # choose kernel with appropriate bandwidth based on p
    if (is.null(Z)) {
        kern <- function(x1, x2) {
            exp(-0.1 * sum((x1 - x2)^2)) * exp(-0.1 * outer(1:nloc, 1:nloc, FUN = "-")^2)
        }
    } else {
        d <- as.matrix(dist(Z))^2
        kern <- function(x1, x2) {
            exp(-0.1 * sum((x1 - x2)^2)) * exp(-0.1 * d)
        }
    }
    
    best.index <- which.max(apply(Xt, 1, FUN = function(x1) { norm(kern(x1, x)) } ))
    init.u <- ntime + runif(nloc)
    init.u <- init.u/inner.norm(init.u) * (2 * quant - 1)
    
    est.mean <- numeric(nloc)
    est.var <- matrix(0, nrow = nloc, ncol = nloc)

    niter <- 0
    converged <- FALSE
    err <- Inf
    q <- runif(nloc)
    while (!converged) {
        wmat <- matrix(0, nrow = nloc, ncol = nloc)
        for (t in 1:ntime) {
            if (abs(t - best.index) < nbn) {
                kmat <- kern(Xt[t, ], x)
                est.mean <- est.mean + kmat %*% Y[t, ]
                est.var <- est.var + kmat %*% tcrossprod(Y[t, ]) 
                wmat <- wmat + kmat
            }
        }
        eps <- 1e-10
        qnew <- solve(wmat + 1e-10*diag(nloc), est.mean + est.var %*% init.u / 2)

        # Error calculation
        err <- sqrt(sum((q - qnew)^2)) / sqrt(sum(q^2))
        q <- qnew
        niter <- niter + 1
        
        # Convergence check
        if ((niter > maxiter) | (err < tol)) {
            converged <- TRUE
        }  
    }

    return(list(
        "q" = qnew
    ))
}

# Same estimator, with Newton Raphson iteration
# faster for low value of p
estimate2.quantile <- function(dat, x, tau = 0.75, maxiter = 1000, tol = 1e-5, rho = 0.5) {
    n <- nrow(dat$Y)
    p <- ncol(dat$Y)
    if (is.null(dim(dat$X))) {
        dim(dat$X) <- c(n, 1)    # convert to column vector
    }
    
    # Creating kernel
    d <- as.matrix(dist(dat$Z))
    kernmat <- function(t1, t2, bandwidth = 1e-2) {
        tdiff <- sqrt(sum((t1 - t2)^2) )
        return(exp(-bandwidth * tdiff * d))
    } 
    
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
            2 * (1/norms[i]) * kernmatList.sq[[i]] %*% dat$Y[i, ]
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
            "q" = q, "niter" = niter, "err" = err
        )
    )
}