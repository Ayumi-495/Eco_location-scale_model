# Load required packages ----
pacman::p_load(
  dplyr, tibble, tidyverse, broom, broom.mixed,
  ape, arm, brms, broom.mixed, cmdstanr, emmeans, glmmTMB, MASS, phytools, rstan, TreeTools,
  DHARMa, loo, MuMIn, parallel, ggeffects,
  bayesplot, ggplot2, patchwork, tidybayes,
  gt, here, kableExtra, knitr
  )

load(here("repro_package", "all_data.RData"))

####
# Model 1 ----
####

dat_tarsus <- all_data$dat_tarsus

## Visualise dataset ----
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


## Run models ----

### Location-only model ----
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

### Location-scale model ----
model_1 <- glmmTMB(
  log(AdTarsus) ~ 1 + Sex + Treatment + Sex:Treatment, # location part
  dispformula = ~ 1 + Sex + Treatment + Sex:Treatment, # scale part
  data = dat_tarsus, 
  family = gaussian
)

summary(model_1)
confint(model_1)

### Model comparison ----
## we can use the anova() function to compare the two models of AIC 
anova(model_0, model_1)

## model.sel() from the MuMIn package can be used to compare AICc values.
model.sel(model_0, model_1)


### Bayesian model ----
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

####
# Model 2 ----
####

dat <- all_data$dat_BFB %>%subset(REDUCTION == 0)

dat<- dat%>%
  dplyr::select(-TIME, -HATCHING_DATE,-REDUCTION,-RING) %>%
  mutate(SMI = as.numeric(SMI),# Scaled mass index
         lnSMI=log(SMI),
         NEST = as.factor(NEST), 
         WORKYEAR = as.factor(WORKYEAR),
         RANK = factor(RANK, levels = c("1", "2")))

## Visualise dataset ----
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


## Run models ----
### Location-only model ----

model2_1<-glmmTMB(lnSMI ~ 1 + RANK + (1|NEST)+(1|WORKYEAR),
                  family = gaussian(), data=dat)

summary(model2_1)
confint(model2_1)

simulationOutput <- simulateResiduals(fittedModel = model2_1, plot = F)
plot(simulationOutput)

### Location-scale model ----

model2_2 <- glmmTMB(
  lnSMI ~ 1 + RANK + (1|NEST)+(1|WORKYEAR),
  dispformula = ~ 1 + RANK,     
  data = dat, 
  family = gaussian
)

summary(model2_2)
confint(model2_2)

### Model comparison ----
model.sel(model2_1, model2_2)

### Bayesian model ----
####  
#'* Note – the models will take some time to finish running. However, you can access the results by downloading and loading the .rds files of the models that have already been run. *
####

#### Location-only model ----
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

#### Location-scale model ----

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

#### Model comparison ----
f1loo <- loo::loo(fit1)
f2loo <- loo::loo(fit2)

fc <- loo::loo_compare(f1loo, f2loo)
fc

#### Plot result ----
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


####
# Model 3 ----
####

## Run model ----
####  
#'* Note – the models will take some time to finish running. However, you can access the results by downloading and loading the .rds files of the models that have already been run. *
####

### Location-scale model ----
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

### Plot result ----
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

####
# Beyond Gaussian 1 ----
####

dat_pref <- all_data$dat_pref

dat_pref <- dat_pref %>%
  dplyr::select(-stripe, -full, -subset) %>%
  rename(frequency = dot) %>%
  mutate(species = phylo, across(c(condition, sex), as.factor))

## Visualise dataset ----
ggplot(dat_pref, aes(x = condition,
                     y = frequency)) +
  geom_violin(color = "#8B8B83", fill = "white",
              width = 1.2, alpha = 0.3) +
  geom_jitter(aes(shape = sex, color = sex),
              size = 3, alpha = 0.8,
              width = 0.15, height = 0) +
  labs(
    title = "Total frequency of gazes towards dot patterns",
    x = "Condition",
    y = "Total frequency of gazes (1hr)"
  ) +
  scale_shape_manual(values = c("M" = 17, "F" = 16),
                     labels = c("M" = "Male", "F" = "Female")) +
  scale_color_manual(values = c("M" = "#009ACD", "F" = "#FF4D4D"),
                     labels = c("M" = "Male", "F" = "Female")) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(color = "#6E7B8B", size = 14),
    axis.title = element_text(color = "#6E7B8B", size = 14),
    legend.title = element_text(color = "#6E7B8B"),
    legend.text = element_text(color = "#6E7B8B"),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(dat_pref, aes(x = condition, y = frequency)) +
  geom_point(alpha = 0.5) +
  labs(
    title = "Total frequency of gazes towards dot patterns by species",
    x = "Condition",
    y = "Total frequency of gazes (1hr)"
  ) +
  stat_summary(fun = mean, geom = "line", aes(group = species, color = species)) +
  theme_classic() +
  facet_wrap(~ species)

## Run models ----
### Location-only model ----
model_0 <- glmmTMB(
  frequency ~ 1 +  condition + sex + (1 | species) + (1 | id), # location part (mean)
  data = dat_pref,
  family = nbinom2(link = "log")
)

summary(model_0)
confint(model_0)

model0_res <- simulateResiduals(model_0, plot = TRUE, seed = 42)
testDispersion(model0_res) 

### Location-scale model ----
# NB model
model_1 <- glmmTMB(
  frequency ~ 1 + condition + sex + (1 | species) + (1 | id),
  dispformula = ~ condition + sex,     
  data = dat_pref,                          
  family = nbinom2(link = "log")
)
summary(model_1)
confint(model_1)

# CMP model
model_2 <- glmmTMB(
  frequency ~ 1 + condition + sex + (1 | species) + (1 | id),
  dispformula = ~ condition + sex, 
  data = dat_pref,
  family = compois(link = "log")
)
summary(model_2)
confint(model_2)

### Model comparison ----
model.sel(model_0, model_1)
model.sel(model_0, model_1, model_2)

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


####
# beyond Gaussian 2 ----
####

dat <- all_data$dat_Cougar

dat <- dat %>%
  dplyr::select(Site, Pool, if_kill, cover) %>%
  mutate(Site=as.factor(Site),
         if_kill=as.factor(if_kill),
         cover=as.numeric(cover)
  )

## Visualise dataset ----

ggplot(dat, aes(x = if_kill, y = cover, fill = if_kill)) +
  geom_violin(
    aes(fill = if_kill), # Fill violins based on 'if_kill'
    color = "#8B8B83", # Outline color for violins
    width = 0.8,
    alpha = 0.3,
    position = position_dodge(width = 0.7)
  ) +
  geom_jitter(
    aes(color = if_kill), # Color jittered points based on 'if_kill'
    size = 3,
    alpha = 0.4,
    shape = 1, # Open circles for jittered points
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.15)
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.1,
    color = "black", # Black crossbar for mean
    linewidth = 0.5,
    position = position_dodge(width = 0.7)
  ) +
  labs(
    title = "Donkey Trampling by Cougar Kill Presence",
    x = "Cougar Kill Presence",
    y = "Proportion Trampled Bare Ground"
  ) +
  scale_fill_manual(
    values = c("burro kills absent" = "cornflowerblue", "burro kills present" = "firebrick")
  ) +
  scale_color_manual( # Add scale_color_manual for jitter points
    values = c("burro kills absent" = "cornflowerblue", "burro kills present" = "firebrick")
  ) +
  scale_x_discrete(
    labels = c("burro kills absent" = "No", "burro kills present" = "Yes"),
    expand = expansion(add = 0.5)
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(color = "#6E7B8B", size = 14),
    axis.title = element_text(color = "#6E7B8B", size = 14),
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

## Run model ----
####
#'* Note – the models will take some time to finish running. However, you can access the results by downloading and loading the .rds files of the models that have already been run. *
####

### Location-only model ----

m0<-bf(cover ~ if_kill + (1|Pool),
       zoi~  if_kill, 
       coi~  if_kill) 
prior1<-default_prior(m0, family=zero_one_inflated_beta(), data=dat)

# Since this model is time-consuming, reload without running it:
fit0 <- brm(
  m0,
  data = dat,
  family = zero_one_inflated_beta(),
  prior = prior1,
  iter = 6000,
  warmup = 1000,
  chains = 2,  cores=2,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  seed = 123,
  refresh = 500
)

# fit0 <- readRDS(here("Rdata", "fit0_BETA_Burros.rds"))
summary(fit0)

### Location-scale model ----
m1<-bf(cover ~ if_kill + (1|Pool),
       zoi ~ if_kill, 
       coi ~ if_kill,
       phi ~ if_kill) 
prior2<-default_prior(m1, family=zero_one_inflated_beta(), data=dat)

fit1 <- brm(
  m1,
  data = dat,
  family = zero_one_inflated_beta(),
  prior = prior2,
  iter = 6000,
  warmup = 1000,
  chains = 2,  cores=2,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  seed = 123,      #
  refresh = 500    #
)

# fit1 <- readRDS(here("Rdata", "fit1_BETA_Burros.rds"))
summary(fit1)

### Model comparison ----
f0loo <- loo::loo(fit0)
f1loo <- loo::loo(fit1)

fc<-loo::loo_compare(f0loo, f1loo)
fc

### Plot results ----
burro_draws<-fit1 %>%
  epred_draws(newdata = expand.grid(if_kill=unique(dat$if_kill)
  ),re_formula = NA, dpar=TRUE) 

# Plot for location
plot_mu <- ggplot(burro_draws, aes(x = factor(if_kill), y = mu, fill = factor(if_kill))) +
  geom_violin(alpha = 0.5, show.legend = FALSE) +
  stat_pointinterval(show.legend = FALSE) +
  labs(
    title = "Location",
    subtitle = "Trampled Bare Ground",
    x = "Cougar Kill Presence",
    y = "Proportion "
  ) +
  theme_classic()+scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0, 1))

# Plot for Zero-Inflation (zoi)
plot_zoi <- ggplot(burro_draws, aes(x = factor(if_kill), y = zoi, fill = factor(if_kill))) +
  geom_violin(alpha = 0.5, show.legend = FALSE) +
  stat_pointinterval(show.legend = FALSE) +
  labs(
    title = "Zero-Inflation (zoi)",
    x = "Cougar Kill Presence",
    y = "Probability of Zero"
  ) +
  theme_classic()+scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0, 1))

# Plot for One-Inflation (coi)
plot_coi <- ggplot(burro_draws, aes(x = factor(if_kill), y = coi, fill = factor(if_kill))) +
  geom_violin(alpha = 0.5, show.legend = FALSE) +
  stat_pointinterval(show.legend = FALSE) +
  labs(
    title = "One-Inflation (coi)",
    x = "Cougar Kill Presence",
    y = "Probability of One"
  ) +
  theme_classic() +scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0, 1))

# Plot for Precision (phi)
plot_phi <- ggplot(burro_draws, aes(x = factor(if_kill), y = phi, fill = factor(if_kill))) +
  geom_violin(alpha = 0.5, show.legend = FALSE) +
  stat_pointinterval(show.legend = FALSE) +
  labs(
    title = "Scale",
    x = "Cougar Kill Presence",
    y = "Precision"
  ) +
  theme_classic()+scale_y_continuous(breaks = seq(0,25,5),limits = c(0, 25))

# --- 4. Combine Plots into a Single Figure ---
# Arrange the four plots in a 2x2 grid
combined_plot <- (plot_mu + plot_zoi) / (plot_coi + plot_phi) &
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Display the final combined plot
combined_plot


####
# Model selection 1 ----
####

dat <- all_data$dat_BFB %>%subset(REDUCTION == 0)

dat<- dat%>%
  dplyr::select(-TIME, -HATCHING_DATE,-REDUCTION,-RING) %>%
  mutate(SMI = as.numeric(SMI),# Scaled mass index
         lnSMI=log(SMI),
         NEST = as.factor(NEST), 
         WORKYEAR = as.factor(WORKYEAR),
         RANK = factor(RANK, levels = c("1", "2")))

## Run models ----
model2_1<-glmmTMB(lnSMI ~ 1 + RANK + (1|NEST)+(1|WORKYEAR),
                  family = gaussian(),data=dat)

simulationOutput <- simulateResiduals(fittedModel = model2_1, plot = F)
plot(simulationOutput)

dat<- dat%>%
  dplyr::mutate(ST_DATE=scale(ST_DATE,center=T,scale = F))

model2_2 <- glmmTMB(
  lnSMI ~ 1 + RANK + (1|NEST)+(1|WORKYEAR),
  dispformula = ~ 1 + RANK,     
  data = dat, 
  family = gaussian
)

model2_3a <- glmmTMB(
  lnSMI ~ 1 + RANK + ST_DATE + (1|NEST)+(1|WORKYEAR),
  dispformula = ~ 1 + RANK,     
  data = dat, 
  family = gaussian
)

model2_3b<- glmmTMB(
  lnSMI ~ 1 + RANK + ST_DATE + RANK*ST_DATE +(1|NEST)+(1|WORKYEAR),
  dispformula = ~ 1 + RANK,     
  data = dat, 
  family = gaussian
)

model.sel(model2_2, model2_3a)
model.sel(model2_3a, model2_3b)
summary(model2_3b)




st_date_seq <- seq(min(dat$ST_DATE), max(dat$ST_DATE), length.out = 50)


emm_location <- emmeans(model2_3b, ~ST_DATE*RANK,
                        component = "cond", type = "response",
                        at= list(ST_DATE = st_date_seq,RANK=unique(dat$RANK)))

df_location <- as.data.frame(emm_location) %>%
  rename(
    mean_log = emmean, 
    lwr = lower.CL,    
    upr = upper.CL     
  )
emm_dispersion <- emmeans(model2_3b, ~ RANK,
                          component = "disp", type = "response")
df_dispersion <- as.data.frame(emm_dispersion) %>%
  rename(
    log_sigma = response, 
    lwr = lower.CL,     
    upr = upper.CL     
  )
pos <- position_dodge(width = 0.4)


plot_location_emmeans <- ggplot(df_location, aes(x = ST_DATE, y = mean_log,color=factor(RANK),group=factor(RANK))) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = factor(RANK)), alpha = 0.2, linetype = 0) +
  geom_line(linewidth = 1)+
  labs(
    title = "Location",
    x = "Hatching date (days)",
    y = "ln(SMI)",
    color="Hatching order",
    fill= "Hatching order"
  ) +
  theme_classic()+ 
  scale_y_continuous(breaks = seq(6.5,8,0.2),limits = c(6.5,8))+
  scale_x_continuous(breaks= seq(0,160,30),limits = c(0, 160)) +theme(legend.position = 'bottom')

plot_dispersion_emmeans <- ggplot(df_dispersion, aes(x = RANK, y = log_sigma,color=factor(RANK))) +
  geom_point(position = pos, size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), position = pos, width = 0.15) +
  labs(
    title = "Scale",
    x = "Hatching order",
    y = "Standard Deviation",
    color="Hatching order"
  ) +
  theme_classic() + scale_y_continuous(breaks = seq(0.07,0.12,0.01),limits = c(0.07, 0.12))+
  theme(legend.position = "none")

plot_location_emmeans + plot_dispersion_emmeans


####
# Model selection 2 ----
####

####  
#'* Note – the models will take some time to finish running. However, you can access the results by downloading and loading the .rds files of the models that have already been run. *
####

# model_eg2.0 <- readRDS(here("Rdata", "model_eg2.0.rds"))
# model_eg2.1 <- readRDS(here("Rdata", "model_eg2.1.rds"))
# model_eg2.2 <- readRDS(here("Rdata", "model_eg2.2.rds"))
# model_eg2.3 <- readRDS(here("Rdata", "model_eg2.3.rds"))

tree <- all_data$pref_tree
tree　<-force.ultrametric(tree, method = "extend")

## Run models ----

model_0 <- glmmTMB(
  frequency ~ 1 +  condition + sex + (1 | species) + (1 | id), # location part (mean)
  data = dat_pref,
  family = nbinom2(link = "log")
)
model0_res <- simulateResiduals(model_0, plot = TRUE, seed = 42)

# formal test for over/underdispersion
testDispersion(model0_res) 

## Bayesian model (location-only) ---- 
formula_eg2.0 <- bf(
  frequency ~ 1 + condition + sex + (1 | species) + (1 | id)
)

prior_eg2.0 <- default_prior(formula_eg2.0, 
                             data = dat_pref, 
                             family = negbinomial(link = "log")
)

model_eg2.0 <- brm(formula_eg2.0, 
                   data = dat_pref, 
                   prior = prior_eg2.0,
                   chains = 2, 
                   iter = 5000, 
                   warmup = 3000,
                   thin = 1,
                   family = negbinomial(link = "log"),
                   save_pars = save_pars(all = TRUE)
                   # control = list(adapt_delta = 0.95)
)


# Create phylogenetic correlation matrix
A <- ape::vcv.phylo(tree, corr = TRUE)

# Specify the model formula (location-only)
formula_eg2.1 <- bf(
  frequency ~ 1 + condition + sex + 
    (1 | a | gr(phylo, cov = A)) +  # phylogenetic random effect
    (1 | species) +             # non-phylogenetic species-level random effect (ecological factors)
    (1 | id)                    # individual-level random effect
)

prior_eg2.1 <- brms::get_prior(
  formula = formula_eg2.1, 
  data = dat_pref, 
  data2 = list(A = A),
  family = negbinomial(link = "log")
)

model_eg2.1 <- brm(
  formula = formula_eg2.1, 
  data = dat_pref, 
  data2 = list(A = A),
  chains = 2, 
  iter = 12000, 
  warmup = 10000,
  thin = 1,
  family = negbinomial(link = "log"),
  prior = prior_eg2.1,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE)
)


### Location-scale models ----
# Specify the model formula (location-scale) - the scale part does not include random effects
formula_eg2.2 <- bf(
  frequency ~ 1 + condition + sex + 
    (1 | a | gr(phylo, cov = A)) +  
    (1 | species) +          
    (1 | id),              
  shape ~ 1 + condition + sex
)

prior_eg2.2 <- brms::get_prior(
  formula = formula_eg2.2, 
  data = dat_pref, 
  data2 = list(A = A),
  family = negbinomial(link = "log", link_shape = "log")
)

model_eg2.2 <- brm(
  formula = formula_eg2.2, 
  data = dat_pref, 
  data2 = list(A = A),
  chains = 2, 
  iter = 12000, 
  warmup = 10000,
  thin = 1,
  family = negbinomial(link = "log", link_shape = "log"),
  prior = prior_eg2.2,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE)
)


formula_eg2.3 <- bf(
  frequency ~ 1 + condition + sex + 
    (1 |a| gr(phylo, cov = A)) +  # phylogenetic random effect
    (1 | species) +            # non-phylogenetic species-level random effect (ecological factors)
    (1 | id),              # individual-level random effect
  shape ~ 1 + condition + sex +
    (1 |a| gr(phylo, cov = A)) +  
    (1 | species) +  
    (1 | id)          
)

prior_eg2.3 <- brms::get_prior(
  formula = formula_eg2.3, 
  data = dat_pref, 
  data2 = list(A = A),
  family = negbinomial(link = "log", link_shape = "log")
)

model_eg2.3 <- brm(
  formula = formula_eg2.3, 
  data = dat_pref, 
  data2 = list(A = A),
  chains = 2, 
  iter = 12000, 
  warmup = 10000,
  thin = 1,
  family = negbinomial(link = "log", link_shape = "log"),
  prior = prior_eg2.3,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE)
)

summary(model_eg2.3)


options(future.globals.maxSize = 2 * 1024^3)
loo_eg2.0 <- loo::loo(model_eg2.0, moment_match = TRUE,
                      reloo = TRUE, cores = 2)
loo_eg2.1 <- loo::loo(model_eg2.1, moment_match = TRUE, 
                      reloo = TRUE, cores = 2)
loo_eg2.2 <- loo::loo(model_eg2.2, moment_match = TRUE,
                      reloo = TRUE, cores = 2)
loo_eg2.3 <- loo::loo(model_eg2.3, moment_match = TRUE,
                      reloo = TRUE, cores = 2)

fc_eg2 <- loo::loo_compare(loo_eg2.0, loo_eg2.1, loo_eg2.2, loo_eg2.3)

print(fc_eg2)


# Create a grid of predictor values for posterior distributions
nd <- tidyr::expand_grid(
  condition = levels(dat_pref$condition),
  sex  = levels(dat_pref$sex)
)

# Get posterior draws for mu and shape and reshape the data
draws_link <- tidybayes::linpred_draws(
  model_eg2.3, 
  newdata = nd, 
  re_formula = NA, 
  dpar = TRUE
) %>%
  transmute(
    sex, condition, .draw,
    log_mu = mu,
    log_theta = shape
  ) %>%
  tidyr::pivot_longer(
    cols = c(log_mu, log_theta), 
    names_to = "parameter", 
    values_to = "value"
  )


# Get the posterior median of the linear predictor for each observation
pred_obs <- add_linpred_draws(
  newdata = dat_pref,    # Use the original data here
  object = model_eg2.3,
  dpar   = "mu",
  re_formula = NULL      # Include group-level effects for specific predictions
) %>%
  group_by(.row, condition, sex) %>% # .row is added by add_linpred_draws
  summarise(
    predicted_logmu = median(.linpred), # Calculate the posterior median
    .groups = "drop"
  )

# Create the plot for the Location parameter (mu) with predicted points
p_loc <- draws_link %>%
  filter(parameter == "log_mu") %>%
  ggplot(aes(x = sex, y = value, fill = condition)) +
  geom_violin(
    width = 0.9, trim = FALSE, alpha = 0.6,
    color = NA, position = position_dodge(0.6)
  ) +
  stat_pointinterval(
    aes(color = condition),
    position = position_dodge(width = 0.6),
    .width = 0.95, size = 0.8
  ) +
  # Add the new layer for predicted data points
  geom_point(
    data = pred_obs,
    aes(y = predicted_logmu, color = condition), # y-aesthetic is the predicted value
    position = position_jitterdodge(
      dodge.width = 0.6,     # Must match the dodge width above
      jitter.width = 0.1
    ),
    alpha = 0.4,
    size = 1.5
  ) +
  labs(
    title = "Location", 
    x = "Sex", 
    y = "log(mean frequency)"
  ) +
  scale_y_continuous(limits = c(1.0, 7.0), breaks = seq(1.0, 7.0, by = 0.5)) +
  theme_classic() +
  theme(legend.position = "right")

# Create the plot for the Scale parameter (shape)
p_scl <- draws_link %>%
  filter(parameter == "log_theta") %>%
  ggplot(aes(x = sex, y = value, fill = condition)) +
  geom_violin(
    width = 0.9, trim = FALSE, alpha = 0.6,
    color = NA, position = position_dodge(0.6)
  ) +
  stat_pointinterval(
    aes(color = condition),
    position = position_dodge(width = 0.6),
    .width = 0.95, size = 0.8
  ) +
  labs(
    title = "Scale", 
    x = "Sex", 
    y = "log(theta)"
  ) +
  scale_y_continuous(limits = c(-2.0, 4.5), breaks = seq(-2.0, 4.0, by = 0.5)) +
  theme_classic() +
  theme(legend.position = "right")

p_loc + p_scl
