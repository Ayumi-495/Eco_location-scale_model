# Load required packages and dataset ----
pacman::p_load(
  dplyr, tibble, tidyverse, broom, broom.mixed,
  ape, arm, brms, broom.mixed, cmdstanr, emmeans, glmmTMB, MASS, phytools, rstan, TreeTools,
  DHARMa, loo, MuMIn, parallel, ggeffects,
  bayesplot, ggplot2, patchwork, tidybayes,
  gt, here, kableExtra, knitr
  )

load(here("repro_package", "all_data.RData"))
dat_tarsus <- all_data$dat_tarsus

# Visualise dataset ----
ggplot(dat_tarsus, aes(x = Treatment, y = log(AdTarsus), fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, position = position_dodge(width = 0.8)) +
  geom_jitter(
    aes(color = Sex),
    size = 2, alpha = 0.7,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  ) +
  scale_fill_manual(values = c("Male" = "#1f78b4", "Female" = "#e31a1c")) + 
  scale_color_manual(values = c("Male" = "#1f78b4", "Female" = "#e31a1c")) +
  labs(title = "Adult tarsus length by treatment and sex",
       x = "Treatment", y = "Log-transformed tarsus length") +
  theme_classic() +
  theme(legend.position = "right")


# Run models ----

## Location-only model ----
model_0 <- glmmTMB(
  log(AdTarsus) ~ 1 + Sex + Treatment + Sex:Treatment,
  data = dat_tarsus, 
  family = gaussian)

summary(model_0)
confint(model_0) # check 95%CI

# plot a q-q plot of residuals to visually assess the normality assumption
res <- residuals(model_0)

qqnorm(res) # visual check for normality of residuals
qqline(res) # reference line for normal distribution

## Location-scale model ----
model_1 <- glmmTMB(
  log(AdTarsus) ~ 1 + Sex + Treatment + Sex:Treatment, # location part
  dispformula = ~ 1 + Sex + Treatment + Sex:Treatment, # scale part
  data = dat_tarsus, 
  family = gaussian
)

summary(model_1)
confint(model_1)

## Model comparison ----
## we can use the anova() function to compare the two models of AIC 
anova(model_0, model_1)

## model.sel() from the MuMIn package can be used to compare AICc values.
model.sel(model_0, model_1)


## Bayesian model ----
####  
#'* Note – the model will take some time to finish running. However, you can access the results by downloading and loading the .rds files of the models that have already been run. *
####

# specify the model using bf()
formula1 <- bf(
  log(AdTarsus) ~  1 + Sex + Treatment + Sex:Treatment, 
  sigma = ~ 1 + Sex + Treatment + Sex:Treatment
)

# generate default priors based on the formula and data
default_priors <- default_prior(
  formula1,
  data = dat_tarsus,                             
  family = gaussian() # default link function for gaussian family                                 
)

# fit the model - you can change N of iter, warmup, thin, and also chains.
# adapt_delta = 0.95 helps to reduce divergent transitions
system.time(
  brms_g1 <- brm(formula1,
                 data = dat_tarsus,           
                 family = gaussian(),                   
                 prior = default_priors,                
                 iter = 2000,                          
                 warmup = 1000, 
                 thin = 1,                              
                 chains = 2,                            
                 control = list(adapt_delta = 0.95) 
  )
)

# brms_g1 <- readRDS(here("Rdata", "brms_SN1.rds"))
summary(brms_g1)
