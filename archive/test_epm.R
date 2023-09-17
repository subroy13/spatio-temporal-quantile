library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(lubridate)
library(dtw)

theme_custom <- function(...) {
    theme_bw() +
        theme(
            axis.title = element_text(size = 14),
            axis.text.y = element_text(size = 11),
            axis.text.x = element_text(size = 11, angle = 0),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 12),
            legend.spacing = unit(5.0, 'cm'),
            plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm"),
            plot.title = element_text(size = 16, hjust = 0.5, vjust = 0.5)
        )
}


###############################
# Read and Explore the Data

rawdf <- read_csv('./data/my_epm1.csv')
df <- rawdf %>% 
    mutate(DateTime = dmy_hm(DateTime), DayCType = 1, DayC_new = DayC) %>%
    spread(DayC_new, DayCType, fill = 0) %>%
    pivot_longer(paste0("F", 1:226), names_to = "Household", values_to = "Demand") %>%
    mutate(Household = as.numeric(stringr::str_remove(Household, "F")))
    
summary(df$DateTime)

# Simple GGPlots as in the book
# Q) What is the total demand across different households?
df %>% 
    group_by(Household) %>% summarise(TotalDemand = sum(Demand)) %>%
    ggplot(aes(x = Household, y = TotalDemand)) + 
    geom_col() 

# Q) What is the max demand across different households?
df %>% 
    group_by(Household) %>% summarise(MaxDemand = max(Demand)) %>%
    ggplot(aes(x = Household, y = MaxDemand)) + 
    geom_col() 

# Q) Boxplot of demands across different households.
df %>% ggplot(aes(x = as.factor(Household), y = Demand)) + geom_boxplot() 

# Q) Average consumption pattern
# Mean pattern <- smoothed
# Median pattern <- much spiky
df %>% group_by(HH) %>% summarise(AvgDemand = mean(Demand), SdDemand = sd(Demand)) %>%
    mutate(Time =  ymd_hm("1970-01-01 00:15") + (HH - 1) * 1800) %>%
    ggplot(aes(x = Time, y = AvgDemand)) +
    geom_ribbon(aes(ymin = pmax(AvgDemand - SdDemand, 0), ymax = AvgDemand + SdDemand), fill = "blue", alpha = 0.1) +
    geom_line() + geom_point()

df %>% group_by(HH) %>% summarise(AvgDemand = median(Demand), SdDemand = mad(Demand)) %>%
    mutate(Time =  ymd_hm("1970-01-01 00:15") + (HH - 1) * 1800) %>%
    ggplot(aes(x = Time, y = AvgDemand)) +
    geom_ribbon(aes(ymin = pmax(AvgDemand - SdDemand, 0), ymax = AvgDemand + SdDemand), fill = "blue", alpha = 0.1) +
    geom_line() + geom_point()


df %>% group_by(DayN) %>% summarise(AvgDemand = mean(Demand)) %>%
    mutate(Time =  DayN) %>%
    ggplot(aes(x = Time, y = AvgDemand)) +
    geom_line() + geom_point()

df %>% group_by(DayW) %>% summarise(AvgDemand = mean(Demand)) %>%
    mutate(Time =  DayW) %>%
    ggplot(aes(x = Time, y = AvgDemand)) +
    geom_line() + geom_point()

df %>% group_by(HH, DayC) %>% 
    summarise(AvgDemand = mean(Demand), SdDemand = sd(Demand)) %>%
    mutate(Time =  ymd_hm("1970-01-01 00:15") + (HH - 1) * 1800) %>%
    ggplot(aes(x = Time, y = AvgDemand, color = DayC)) +
    geom_line() + geom_point()

# Maybe first we need to cluster Households according to common consumption pattern
df %>% 
    group_by(HH, Household) %>%
    summarise(AvgDemand = mean(Demand)) %>%
    mutate(Time =  ymd_hm("1970-01-01 00:15") + (HH - 1) * 1800) %>%
    ggplot(aes(x = Time, y = AvgDemand, color = as.factor(Household))) + 
    geom_line() + theme(legend.position = "none")




# Q) Number of Days
# There are 3168 observations
# Train Time = (July 19 - Sept 21)
df2 <- df %>% filter(DateTime < as.POSIXct("2015-09-21 05:30:00")) # tz-adjustment 
summary(df2$DateTime)

traindf <- df2 %>% filter(DateTime < as.POSIXct("2015-09-14 05:30:00"))  # first 8 week for training
testdf <- df2 %>% filter(DateTime >= as.POSIXct("2015-09-14 05:30:00"))  # ninth week for prediction

# calculate dynamic time warping distances
Xmat <- t(as.matrix(traindf %>% 
                        group_by(Household, HH) %>%
                        summarise(AvgDemand = mean(Demand)) %>%
                        pivot_wider(names_from = "Household", values_from = "AvgDemand") %>%
                        arrange(HH) %>%
                        select(as.character(1:226))
))
d <- dtwDist(Xmat)

# apply Hclustering
hc <- hclust(as.dist(d))
plot(hc)

gplots::heatmap.2(log(1 + d),dendrogram='none', Rowv=TRUE, Colv=TRUE,trace='none')


sqrt(nrow(traindf) / 226)
table(cutree(hc, k = 12))   # ?? Good criterion?

# apply grouping tables on households
grouplabels <- cutree(hc, k = 12)
df2 <- df2 %>% mutate(Cluster = grouplabels[Household], 
                      Train = (DateTime < as.POSIXct("2015-09-14 05:30:00") ))
write_csv(df2, './data/my_epm1_processed.csv')

#################################
# Now, estimate the quantile levels for different clusters
# load the estimator
source('./scripts/estimator.R')

df <- read_csv('./data/my_epm1_processed.csv')
traindf <- df %>% filter(Train == TRUE)
testdf <- df %>% filter(Train == FALSE)

Xtrain <- as.matrix(traindf %>% select(WeekT, DayN, HH, DayW, hh, hw, wh, ww) %>% distinct())
Xtest <- as.matrix(testdf %>% select(WeekT, DayN, HH, DayW, hh, hw, wh, ww) %>% distinct())
Xfull <- rbind(Xtrain, Xtest)

colSums(table(testdf$Household, testdf$Cluster) > 1)

{
    clustNum <- 6
    quants <- c(0.5, 0.75, 0.8, 0.9, 0.95, 0.99)
    dflist <- list()
    
    for (i in 1:length(quants)) {
        quant <- quants[i]
        message(paste("Processing quantile level = ", quant))
        
        # Subsetting the training and testing part
        Ytrain <- traindf %>% 
            filter(Cluster == clustNum) %>%
            pivot_wider(names_from = "Household", values_from = "Demand") %>%
            select(-c(DateTime, WeekT, DayN, HH, DayW, Season, DayC, hh, hw, wh, ww, Cluster, Train)) %>%
            as.matrix()
        Ytest <- testdf %>% 
            filter(Cluster == clustNum) %>%
            pivot_wider(names_from = "Household", values_from = "Demand") %>%
            select(-c(DateTime, WeekT, DayN, HH, DayW, Season, DayC, hh, hw, wh, ww, Cluster, Train)) %>%
            as.matrix()
        
        Yfull <- df %>% 
            filter(Cluster == clustNum) %>%
            pivot_wider(names_from = "Household", values_from = "Demand") %>%
            select(-c(DateTime, WeekT, DayN, HH, DayW, Season, DayC, hh, hw, wh, ww, Cluster, Train)) %>%
            as.matrix()
        
        
        EstMat <- matrix(NA, nrow(Xfull), ncol(Yfull))
        colnames(EstMat) <- paste("Household",unique(traindf$Household[traindf$Cluster == clustNum]))
        
        pb <- txtProgressBar(max = nrow(Xfull), style = 3)
        for (j in 1:nrow(Xfull)) {
            EstMat[j, ] <- estimate2.quantile(Ytrain, Xtrain, quant = quant, x = Xfull[j, ])
            setTxtProgressBar(pb, value = j)
        }
        close(pb)
        
        dflist[[i]] <- as_tibble(EstMat)
        dflist[[i]]$Quantile = quant
        dflist[[i]]$Datetime <- sort(unique(df$DateTime))
    }
    
    write_csv(bind_rows(dflist), paste0('./scripts/output/epm/quants-cluster', clustNum,'.csv'))
}

###################################
# Do some plotting

df <- read_csv('./data/my_epm1_processed.csv')
traindf <- df %>% filter(Train == TRUE)
testdf <- df %>% filter(Train == FALSE)
times <- sort(unique(testdf$DateTime))

clust <- 1
qdf <- read_csv(paste0('./scripts/output/epm/quants-cluster',clust,'.csv'))

qdf %>% 
    mutate(Time = rep(times, 6), Quantile = as.factor(Quantile)) %>%
    pivot_longer(-c(Quantile, Time), names_to = "HH", values_to = "Demand") %>%
    group_by(Quantile, Time) %>%
    summarise(Demand = mean(Demand)) %>%
    ggplot(aes(x = Time, y = Demand, color = Quantile)) +
    geom_line(size = 1) + 
    theme_custom() +
    ggtitle(paste("Estimated Quantile Levels for Electricity Demand for Cluster", clust))

ggsave(paste0('./scripts/figure/quantlevel-',clust,'.png'))





























