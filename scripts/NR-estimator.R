# Dummy Data Generation Function
# X_t = t
# Y_t = (1 + 2*sin(Xt), Xt + cos(Xt), Xt^2 - 4*Xt) + epsilon
library(matrixStats)


gen.data <- function(n = 100, seed = NULL) {
    # generate for 100 time points
    if (!is.null(seed)) {
        set.seed(seed)
    }
    Xt <- 1:n
    Yt <- cbind(
        Xt + 2*sin(Xt),
        Xt + cos(Xt),
        2 * Xt - 1
    )
    Yt <- Yt + 0.1 * rnorm(length(Xt) * 3)
    return(list("X" = matrix(Xt, nrow = n, ncol = 1), "Y" = Yt))
}

inner.prod <- function(x, y) { sum(x*y) }
inner.norm <- function(x) {sqrt(sum(x^2))}

dat <- gen.data(n = 50)


# Perform the NR Algorithm to obtain the estimate
NR.estimate <- function(dat, x, u = rep(0, ncol(dat$Y)), 
                        kernmat_func = NULL, tol = 1e-5, maxiter = 1e3) {
    ntime <- nrow(dat$Y)
    nloc <- ncol(dat$Y)
    if (is.null(kernmat_func)) {
        kernmat_func <- function(x1, x2) {
            d <- inner.norm(x1 - x2)^2
            k <- outer(1:nloc, 1:nloc)
            k <- exp(-0.1 * d ) * exp(-0.1 * k^2)
            return(k)
        }
    }
    if (sum(u^2) > 1) {
        stop("u must be within unit circle")
    }
    q <- colQuantiles(dat$Y, probs = 0.5, na.rm = T)
    err <- Inf
    niter <- 1
    while ((err > tol) & (niter < maxiter)) {
        #print(niter)
        scorevec <- numeric(nloc)
        Jmat <- matrix(0, nrow = nloc, ncol = nloc)
        Z <- matrix(NA, nrow = ntime, ncol = nloc)
        for (t in 1:ntime) {
            kmat <- kernmat_func(x, dat$X[t,]) / sqrt(ntime)
            Z[t, ] <- kmat %*% (dat$Y[t, ] - q)
            scorevec <- scorevec + kmat %*% (u + (Z[t, ])/ inner.norm(Z[t, ]) )
            Jmat <- Jmat + kmat %*% (
                tcrossprod(Z[t, ]) / inner.norm(Z[t, ])^3 - diag(nloc) / inner.norm(Z[t, ])
            ) %*% kmat
        }
        newq <- q - solve(Jmat + diag(nloc) * 0 , scorevec)   # add penalization
        err <- inner.norm(q-newq)
        niter <- niter + 1
        q <- newq
    }
    cat(paste("Iteration: ", niter, "\tError: ", err, "\n"))
    return(q)
}


x <- as.vector(5)
NR.estimate(dat, x)

####################################
# Using Optim

obj <- function(q, dat, x, u) {
    ntime <- nrow(dat$Y)
    nloc <- ncol(dat$Y)
    kernmat_func <- function(x1, x2) {
            d <- inner.norm(x1 - x2)
            k <- outer(1:nloc, 1:nloc)
            k <- exp(-0.1 * d ) * exp(-0.01 * k^2)
            return(k)
    }
    
    L <- 0
    for (t in 1:ntime) {
        kmat <- kernmat_func(x, dat$X[t,]) / sqrt(ntime)
        Z <- kmat %*% (dat$Y[t, ] - q)
        L <- L + inner.norm(Z) + inner.prod(Z, u)
    }
    return(log(L))
}

init.q <- colQuantiles(dat$Y, probs = 0.5, na.rm = T)
optim(init.q, fn = obj, dat = dat, x = 5, u = -rep(1,3)/3)$par





















