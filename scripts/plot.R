library(dplyr)
library(readr)
library(ggplot2)

# Custom ggplot theme
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


df <- read_csv('./output/hypothesis-sim1.csv')
df$quant[df$quant == 0.501] <- 0.5
df$logp <- ifelse(df$pval > 0, log(df$pval), log(1e-10))

df %>%
    filter(sigma <= 3) %>%
    ggplot(aes(x = sigma, y = logp, color = as.factor(quant))) +
    geom_line(linewidth = 1) +
    xlab(expression(sigma)) +
    ylab('Log p-value') +
    guides(color = guide_legend(title = 'Quantile')) +
    theme_custom()

ggsave(paste0('./figures/hypothesis-sim1.png'), width = 6, height = 3.5)



## Same thing for 2nd simulation
df <- read_csv('./output/hypothesis-sim2.csv')
df$quant[df$quant == 0.501] <- 0.5
df$logp <- ifelse(df$pval > 0, log(df$pval), log(1e-10))

df %>%
    ggplot(aes(x = sigma, y = logp, color = as.factor(quant))) +
    geom_line(linewidth = 1) +
    xlab(expression(sigma)) +
    ylab('Log p-value') +
    guides(color = guide_legend(title = 'Quantile')) +
    theme_custom()

ggsave(paste0('./figures/hypothesis-sim2.png'), width = 6, height = 3.5)