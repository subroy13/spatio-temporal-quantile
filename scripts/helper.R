##########################
# Some helper functions
#########################

# function to compute inner norm of a vector
inner.norm <- function(x) sqrt(sum(x^2))

# function to compute inner product of two vectors
inner.prod <- function(x, y) sum(x*y)

# matrix square root function of symmetric matrices, or inverse square root
matsq <- function(X, inverse = FALSE) {
    s <- svd(X, symmetric = TRUE)
    evals.sq <- ifelse(s$d > 0, ifelse(inverse, s$d^(-0.5), s$d^0.5), 0)
    return(tcrossprod( t(t(s$u) * evals.sq), s$v ))
}

# folded normal distribution related normalizing function C_t(z)
folded.crit <- function(t, z) {
    if (t <= 1) {
        stop("Invalid value of t")
    }
    return((z - (log(log(t)) + log(4 * pi))/2) / sqrt(2 * log(t)))
}

# inverse of folded normal critical value find z such that C_t(z) = X
inv.folded.crit <- function(t, x) {
    if (t <= 1) {
        stop("Invalid value of t")
    }
    return(x * sqrt(2 * log(t)) + (log(log(t)) + log(4*pi))/2)
}

