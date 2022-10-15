#############################
# Main Simulation Loop

library(matrixStats)
library(parallel)
library(pbapply)
library(dplyr)
library(readr)
library(tidyr)
source('./scripts/estimator.R')

##################################
# Model Generating Function

library(matrixStats)
library(mlVAR)

# Model-1
# Xt ~ VAR(1) process
sim.model1 <- function(ntime, nloc, B = 1000, seed = 1) {
    set.seed(seed)
    
    # generate the parameters
    params.VAR <- matrix(c(0.2, 0.1, -0.3, 0.4), nrow = 2, ncol = 2)
    params.slope <- matrix(runif( 3*nloc, min = -1, max = 1), nrow = 3, ncol = nloc)
    params.slope[-1, ] <- params.slope[-1, ] * 20
    
    params.intercept <- rnorm(nloc)
    Xt <- as.matrix(cbind(1:ntime, simulateVAR(params.VAR, Nt = ntime )))
    
    Ydata <- list()
    pb <- txtProgressBar(max = B, style = 3)
    for (b in 1:B) {
        Y <- matrix(NA, nrow = ntime, ncol = nloc)
        for (t in 1:ntime) {
            mu <- params.intercept + Xt[t, 1] * params.slope[1, ] + sin(2 * pi * Xt[t, 2]) * params.slope[2] + cos(2 * pi * Xt[t, 3]) * params.slope[3]
            
            # reduce error variance so that the signal is prominent
            Sigma <- exp(-0.1 * outer(1:nloc, 1:nloc, FUN = "-")^2) * mean(abs(mu)) * 0.1
            Y[t, ] <- MASS::mvrnorm(mu = mu, Sigma = Sigma)
        }
        Ydata[[b]] <- Y
        
        setTxtProgressBar(pb, value = b)
    }
    close(pb)
    
    return(list("Ydata" = Ydata, "Xt" = Xt, "mu" = mu,
                "params" = list("var" = params.VAR, "slope" = params.slope, 
                                "intercept" = params.intercept )))
}



###############################################

# Simulate data from Model-1 (Exp-1)
dat <- sim.model1(ntime = 300, nloc = 10, B = 1000, seed = 1234)
saveRDS(dat, './scripts/output/model1-exp1-data.Rds')

simulation <- function(dat, params, nsim = 100, seed = 1234, holdout = 0) {
    # Get the true parameters first
    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl, c("true.quantile", "dat", "params", 
                        "rowMedians", "inner.norm", "inner.prod"))
    truetheta <- pblapply(1:nrow(params), FUN = function(i) {
        true.quantile(dat$Ydata, dat$Xt, params$alpha[[i]], params$x[i, ])
    }, cl = cl)
    closeAllConnections()
    
    # Do nsim simulations
    set.seed(seed)
    index <- sample(1:length(dat$Ydata), size = nsim)
    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl, c("estimate.quantile", "dat", "params", "rowMedians", "colMedians",
                        "inner.norm", "inner.prod"))
    
    n <- nrow(dat$Ydata[[1]])
    holdout_start <- floor(n * (1 - holdout))
    
    esttheta <- pblapply(1:nsim, FUN = function(i){
        id <- index[i]
        x <- lapply(1:nrow(params), FUN = function(k) {
            estimate.quantile(dat$Ydata[[id]][1:holdout_start,], 
                              dat$Xt[1:holdout_start,], 
                              params$alpha[[k]], 
                              params$x[k, ])
        })
        return(x)
    }, cl = cl)
    closeAllConnections()
    
    return(list("truth" = truetheta, "estimate" = esttheta, "params" = params))
}


calc.RMSE <- function(res) {
    nparams <- length(res$truth)
    nsim <- length(res$estimate)
    mse <- matrix(NA, nrow = nparams, ncol = nsim)
    pe <- matrix(NA, nrow = nparams, ncol = nsim)
    mae <- matrix(NA, nrow = nparams, ncol = nsim)
    
    for (b in 1:nsim) {
        for (l in 1:nparams) {
            mse[l, b] <- mean(
                (res$truth[[l]] - res$estimate[[b]][[l]])^2
            )
            mae[l, b] <- mean(
                abs(res$truth[[l]] - res$estimate[[b]][[l]])
            )
            pe[l, b] <- inner.norm(
                (res$truth[[l]] - res$estimate[[b]][[l]])
            ) / inner.norm(res$estimate[[b]][[l]]) * 100
        }
    }
    
    df <- res$params
    df$RMSE <- sqrt(rowMeans(mse))
    df$MAE <- rowMeans(mae)
    df$MAPE <- rowMeans(pe)
    return(df)
}



{
    # Let's find the true quantile for 
    #   alpha = 0.5, 0.75, 0.8, 0.85, 0.9, 0.95
    #   t = 10, 50, 150, 250, 290
    dat <- readRDS('./scripts/output/model1-exp1-data.Rds')
    params <- expand.grid(alpha = c(0.5, 0.75, 0.8, 0.85, 0.9, 0.95), 
                          time = c(10, 50, 150, 250, 290) )
    params$x <- dat$Xt[params$time, ]
    
    res <- simulation(dat, params, 100)
    saveRDS(res, './scripts/output/model1-exp1-sim.Rds')
}

# Calculate RMSE
res <- readRDS('./scripts/output/model1-exp1-sim.Rds')
calc.RMSE(res) %>%
    pivot_wider(id_cols = "alpha", names_from = "time", values_from = "MAPE") %>%
    xtable::xtable()




{
    # Let's find the true and predicted quantile for 
    #   alpha = 0.5, 0.75, 0.8, 0.85, 0.9, 0.95
    #   t = 296,... 300
    dat <- readRDS('./scripts/output/model1-exp1-data.Rds')
    params <- expand.grid(alpha = c(0.5, 0.75, 0.8, 0.85, 0.9, 0.95), time = 296:300 )
    params$x <- dat$Xt[params$time, ]
    
    res <- simulation(dat, params, 100, holdout = (5/300))
    saveRDS(res, './scripts/output/model1-exp1-simpred.Rds')
}


# Calculate RMSE
res <- readRDS('./scripts/output/model1-exp1-simpred.Rds')
calc.RMSE(res) %>%
    mutate(OUT = paste0(round(RMSE,3), " (", round(MAPE, 3), ")")) %>%
    pivot_wider(id_cols = "alpha", names_from = "time", values_from = "OUT") %>%
    xtable::xtable()






























