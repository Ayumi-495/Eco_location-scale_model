# Load required packages and dataset ----
pacman::p_load(
  dplyr, tibble, tidyverse, broom, broom.mixed,
  ape, arm, brms, broom.mixed, cmdstanr, emmeans, glmmTMB, MASS, phytools, rstan, TreeTools,
  DHARMa, loo, MuMIn, parallel, ggeffects,
  bayesplot, ggplot2, patchwork, tidybayes,
  gt, here, kableExtra, knitr
)

load(here("repro_package", "all_data.RData"))
dat <- all_data$dat_Cougar

dat <- dat %>%
  dplyr::select(Site, Pool, if_kill, cover) %>%
  mutate(Site=as.factor(Site),
         if_kill=as.factor(if_kill),
         cover=as.numeric(cover)
  )

# Visualise dataset ----

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

# Run model ----
####
#'* Note – the models will take some time to finish running. However, you can access the results by downloading and loading the .rds files of the models that have already been run. *
####

## Location-only model ----

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

## Location-scale model ----
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

## Model comparison ----
f0loo <- loo::loo(fit0)
f1loo <- loo::loo(fit1)

fc<-loo::loo_compare(f0loo, f1loo)
fc

## Plot results ----
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