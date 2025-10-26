# Load required packages and dataset ----
pacman::p_load(
  dplyr, tibble, tidyverse, broom, broom.mixed,
  ape, arm, brms, broom.mixed, cmdstanr, emmeans, glmmTMB, MASS, phytools, rstan, TreeTools,
  DHARMa, loo, MuMIn, parallel, ggeffects,
  bayesplot, ggplot2, patchwork, tidybayes,
  gt, here, kableExtra, knitr
)

load(here("repro_package", "all_data.RData"))
dat <- all_data$dat_BFB %>%subset(REDUCTION == 0)

dat<- dat%>%
  dplyr::select(-TIME, -HATCHING_DATE,-REDUCTION,-RING) %>%
  mutate(SMI = as.numeric(SMI),# Scaled mass index
         lnSMI=log(SMI),
         NEST = as.factor(NEST), 
         WORKYEAR = as.factor(WORKYEAR),
         RANK = factor(RANK, levels = c("1", "2")))

# Run model ----
####  
#'* Note – the models will take some time to finish running. However, you can access the results by downloading and loading the .rds files of the models that have already been run. *
####

## Location-scale model ----
m3 <- bf(log(SMI) ~ 1 + RANK+ (1|q|NEST) + (1|WORKYEAR),
         sigma~ 1 + RANK + (1|q|NEST) + (1|WORKYEAR))
prior3<-default_prior(m3,data = dat,family = gaussian())

fit3 <- brm(
  m3,
  prior= prior3,
  data = dat,
  family = gaussian(),
  iter = 6000,     
  warmup = 1000,   
  chains = 4,  
  cores=4,
  backend = "cmdstanr",
  control = list(
    adapt_delta = 0.99,  
    max_treedepth = 15  
  ),
  seed = 123,      
  refresh = 500)

# fit3 <- readRDS(here("Rdata", "mod3RED.rds"))
summary(fit3)

## Plot result ----
m3_draws<-fit3 %>%
  tidybayes::epred_draws(newdata = expand.grid(
    RANK = unique(dat$RANK)
  ), re_formula = NA, dpar =TRUE) 


plot_lm3 <- ggplot(m3_draws, aes(x = factor(RANK), y = .epred,fill = factor(RANK))) +
  geom_violin(alpha = 0.5, show.legend = FALSE)+
  stat_pointinterval(show.legend = FALSE) +
  labs(
    title = "Location",
    x = "Hatching order", 
    y = "log(SMI)"
  ) +
  theme_classic() + scale_y_continuous(breaks = seq(7.2,7.5,0.05),limits = c(7.2, 7.5))


plot_sm3 <- ggplot(m3_draws, aes(x = factor(RANK), y = sigma,fill = factor(RANK))) +
  geom_violin(alpha = 0.5, show.legend = FALSE)+
  stat_pointinterval(show.legend = FALSE) +
  labs(
    title = "Scale",
    x = "Hatching order", 
    y = "log(Standard Deviation)"
  ) +
  theme_classic() + scale_y_continuous(breaks = seq(0.04,0.12,0.01),limits = c(0.04, 0.12))


correlation_draws <- fit3 %>%
  spread_draws(cor_NEST__Intercept__sigma_Intercept)





corr_plot<-ggplot(data=correlation_draws,aes(x = cor_NEST__Intercept__sigma_Intercept)) +
  stat_halfeye() + 
  labs(
    title = "Correlation between location and the scale intercepts in the NEST group",
    x = "cor(Intercept, sigma_Intercept)",
    y = "Density"
  ) +
  theme_classic()



combinded_plot<-((plot_lm3 + plot_sm3)/corr_plot)& 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
combinded_plot
