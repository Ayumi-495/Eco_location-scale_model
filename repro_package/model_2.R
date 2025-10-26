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

# Visualise dataset ----
ggplot(dat, aes(x = RANK, y = SMI, fill = RANK, color = RANK)) +
  geom_violin(aes(fill = RANK),
              color = "#8B8B83",
              width = 0.8, 
              alpha = 0.3,
              position = position_dodge(width = 0.7)) +
  geom_jitter(aes(color = RANK),
              size = 3,
              alpha = 0.4,
              shape = 1,
              position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.15)) + 
  stat_summary(fun = mean,              
               geom = "crossbar",       
               width = 0.1,             
               color = "black",         
               linewidth = 0.5,         
               position = position_dodge(width = 0.7))+
  labs(
    title = "Body Condition (SMI) by Hatching Order",
    x = "Hatching Order",
    y = "Scaled Mass Index (g)"
  ) +
  scale_fill_manual(
    values = c("1" = "#1F78B4", "2" = "#E31A1C") # 
  ) +
  scale_color_manual(
    values = c("1" = "#1F78B4", "2" = "#E31A1C") # 
  ) +
  scale_x_discrete(
    labels = c("1" = "First-Hatched", "2" = "Second-Hatched"),
    expand = expansion(add = 0.5)
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(color = "#6E7B8B", size = 14),
    axis.title = element_text(color = "#6E7B8B", size = 14),
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )


# Run models ----
## Location-only model ----

model2_1<-glmmTMB(lnSMI ~ 1 + RANK + (1|NEST)+(1|WORKYEAR),
                  family = gaussian(), data=dat)

summary(model2_1)
confint(model2_1)

simulationOutput <- simulateResiduals(fittedModel = model2_1, plot = F)
plot(simulationOutput)

## Location-scale model ----

model2_2 <- glmmTMB(
  lnSMI ~ 1 + RANK + (1|NEST)+(1|WORKYEAR),
  dispformula = ~ 1 + RANK,     
  data = dat, 
  family = gaussian
)

summary(model2_2)
confint(model2_2)

## Model comparison ----
model.sel(model2_1, model2_2)

## Bayesian model ----
####  
#'* Note – the models will take some time to finish running. However, you can access the results by downloading and loading the .rds files of the models that have already been run. *
####

### Location-only model ----
m1 <- bf(log(SMI) ~ 1 + RANK + (1|NEST) + (1|WORKYEAR))
prior1<-default_prior(m1, data = dat, family = gaussian())

fit1 <- brm(
  m1,
  prior = prior1,
  data = dat,
  family = gaussian(),
  iter = 6000,     
  warmup = 1000,   
  chains = 4,  cores=4,
  backend = "cmdstanr",
  control = list(
    adapt_delta = 0.99,  # Keep high if you have divergent transitions
    max_treedepth = 15   # Keep high if hitting max_treedepth warnings
  ),
  seed = 123,      
  refresh = 500    # Less frequent progress updates (reduces overhead)
)

# fit1 <- readRDS(here("Rdata","mod1RED.rds"))
summary(fit1)

### Location-scale model ----

m2 <- bf(log(SMI) ~ 1 + RANK + (1|NEST) + (1|WORKYEAR),
         sigma~ 1 + RANK)
prior2 <- default_prior(m2,data = dat,family = gaussian())

fit2 <- brm(
  m2,
  prior= prior2,
  data = dat,
  family = gaussian(),
  iter = 6000,     
  warmup = 1000,   
  chains = 4,  cores=4,
  backend = "cmdstanr",
  control = list(
    adapt_delta = 0.99,  # Keep high if you have divergent transitions
    max_treedepth = 15   # Keep high if hitting max_treedepth warnings
  ),
  seed = 123,      
  refresh = 500    # Less frequent progress updates (reduces overhead)
)

# fit2 <- readRDS(here("Rdata","mod2RED.rds"))
summary(fit2)

### Model comparison ----
f1loo <- loo::loo(fit1)
f2loo <- loo::loo(fit2)

fc <- loo::loo_compare(f1loo, f2loo)
fc

## Plot result ----
m2_draws<-fit2 %>%
  tidybayes::epred_draws(newdata = expand.grid(
    RANK = unique(dat$RANK)
  ), re_formula = NA, dpar =TRUE) 


plot_lm2 <- ggplot(m2_draws, aes(x = factor(RANK), y = .epred,fill = factor(RANK))) +
  geom_violin(alpha = 0.5, show.legend = FALSE)+
  stat_pointinterval(show.legend = FALSE) +
  labs(
    title = "Location",
    x = "Hatching order", 
    y = "log(SMI)"
  ) +
  theme_classic() + scale_y_continuous(breaks = seq(7.2,7.5,0.05),limits = c(7.2, 7.5))


plot_sm2 <- ggplot(m2_draws, aes(x = factor(RANK), y = sigma,fill = factor(RANK))) +
  geom_violin(alpha = 0.5, show.legend = FALSE)+
  stat_pointinterval(show.legend = FALSE) +
  labs(
    title = "Scale",
    x = "Hatching order", 
    y = "log(Standard Deviation)"
  ) +
  theme_classic() + scale_y_continuous(breaks = seq(0.07,0.12,0.01),limits = c(0.07, 0.12))

combinded_plot<-(plot_lm2 + plot_sm2)& 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
combinded_plot