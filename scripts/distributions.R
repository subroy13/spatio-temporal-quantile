# Verify the distributional quantities

B <- 10000
d <- 5
X <- MASS::mvrnorm(B, mu = rep(0,d), Sigma = diag(d))

y1 <- t(apply(X, MARGIN = 1, FUN = function(x) { x / sqrt(sum(x^2)) } ))
colMeans(y1)  # E(X/||X||) = 0, symmetric


ds <- 2:25
res <- c()
for (d in ds) {
    X <- MASS::mvrnorm(B, mu = rep(0,d), Sigma = diag(d))
    y2 <- 1/sqrt(rowSums(X^2))
    res <- c(res, mean(y2))    # E(1/||X||) 
}

plot(ds, res, type = "b")
points(ds, gamma((ds-1)/2)/ (sqrt(2) * gamma(ds/2)), type = "b", col = "blue")



