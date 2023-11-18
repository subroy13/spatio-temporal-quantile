#####################
# Idea: X ~ F, then F(X) ~ U, 
# For quantile, we know Q(F(X)) = X, so, Q(U) = X ~ F, where U ~ Unif(0, 1)
############################################

##################
# simulation of model for H01
# x ~ N(0, 1)
# Q(tau, x) = (-0.1 + x + 0.1x^2 + tau * abs(x))(1, .... , 1)
##################
n <- 100
p <- 3
sigma.err <- 1e-1  # increase the error variance to see the effect
quant_level <- 0.75
set.seed(1234)

X <- rnorm(n)  # x-covariate
U <- runif(n)
Y.true <- matrix(
    (-0.1 + X + 0.1 * X^2 + (2 * U - 1) * abs(X)),
    nrow = n, ncol = p, byrow = FALSE
)
err <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- Y.true + sigma.err * err

test.type1(Y, as.matrix(X), quant = quant_level)


##################
# simulation of model for H02
# Q(tau, t) = t * (tau, -tau, tau^2 - tau)
##################
n <- 100
p <- 3
sigma.err <- 1e-1  # increase the error variance to see the effect
quant_level <- 0.75
set.seed(1234)

X <- (1:n)/n  # x-covariate
U <- runif(n)
Y.true <- matrix(NA, nrow = n, ncol = p)
for (i in 1:n) {
    tau <- 2 * U[i] - 1
    Y.true[i, ] <- X[i] * c(tau, -tau, tau^2 - tau)
}
err <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- Y.true + sigma.err * err

test.type2(Y, quant = quant_level)


