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
p <- 10

sigmas <- c(0, 0.1, 0.5, 1, 1.5, 2, 5, 10)  # increase the error variance to see the effect
quant_levels <- c(0.5, 0.75, 0.9, 0.95, 0.99)
df <- expand.grid(sigma = sigmas, quant = quant_levels)
df$pval <- NULL

for (i in 1:nrow(df)) {
    set.seed(1234)
    cat('Iteration: ', i, '\n')
    sigma.err <- df$sigma[i]
    quant_level <- df$quant[i]

    X <- 10 * rnorm(n)  # x-covariate
    U <- runif(n)
    Y.true <- matrix(
        (-0.1 + X + 0.1 * X^2 + (2 * U - 1) * abs(X)),
        nrow = n, ncol = p, byrow = FALSE
    )

    err <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- Y.true + sigma.err * err

    res <- test.type1(Y, as.matrix(X), quant = quant_level)
    df$pval[i] <- res$pvalue
}

write.csv(df, './output/hypothesis-sim1.csv')


##################
# simulation of model for H02
# Q(tau, t) = t * (tau, -tau, tau^2 - tau)
##################
n <- 100
p <- 10

sigmas <- c(0, 0.5, by = 0.1)  # increase the error variance to see the effect
quant_levels <- c(0.5, 0.75, 0.9, 0.95, 0.99)
df <- expand.grid(sigma = sigmas, quant = quant_levels)
df$pval <- NULL

for (i in 1:nrow(df)) {
    set.seed(1234)
    sigma.err <- df$sigma[i]  # increase the error variance to see the effect
    quant_level <- df$quant[i]

    X <- (1:n)/n  # x-covariate
    U <- runif(n)
    Y.true <- matrix(NA, nrow = n, ncol = p)
    for (j in 1:n) {
        tau <- 2 * U[j] - 1
        Y.true[j, ] <- X[j] * c(tau, -tau, tau^2 - tau)
    }
    err <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- Y.true + sigma.err * err

    df$pval[i] <- test.type2(Y, quant = quant_level)
}

write.csv(df, './output/hypothesis-sim2.csv')



############
# Now apply for the dataset

library(readr)
library(dplyr)
library(tidyr)

df <- read_csv('./data/EPM Data/my_epm1_processed.csv')
colSums(table(df$Household, df$Cluster) > 1)

clustNum <- 12   # change the cluster number here
quant_level <- 0.75  # change the quantile level here

Xfull <- df %>% 
  filter(Cluster == clustNum) %>%
  select(WeekT, DayN, HH, DayW, hh, hw, wh, ww) %>% distinct() %>%
  as.matrix()

# using only HH as covariate
Xfull <- df %>% 
  filter(Cluster == clustNum) %>%
  select(HH) %>% distinct() %>%
  as.matrix()

dim(Xfull)

Yfull <- df %>% 
  filter(Cluster == clustNum) %>%
  pivot_wider(names_from = "Household", values_from = "Demand") %>%
  select(-c(DateTime, WeekT, DayN, HH, DayW, Season, DayC, hh, hw, wh, ww, Cluster, Train)) %>%
  as.matrix()
dim(Yfull)

# do the hypothesis testing at specified quantile level
test.type1(Yfull, Xfull, quant = quant_level, verbose = TRUE, tgap = 10)
