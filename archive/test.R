library(mvtnorm)

B <- 1000
p <- 25
qq <- seq(-2, 2, 0.01)
qn <- length(qq)

normexp <- numeric(qn)
for (i in 1:qn) {
    u <- (2 * qq[i] - 1) * rep(1, p)
    X <- mvtnorm::rmvnorm(B, mean = -u, sigma = diag(p))
    normexp[i] <- mean(sqrt(rowSums(X^2)) + rowSums(X) * (2 * 0.75 - 1) / sqrt(p) )
}

plot(qq, normexp, type = "l")
qq[which.min(normexp)]





