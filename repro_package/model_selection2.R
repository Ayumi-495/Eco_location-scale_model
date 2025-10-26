####  
#'* Note – the models will take some time to finish running. However, you can access the results by downloading and loading the .rds files of the models that have already been run. *
####

# model_eg2.0 <- readRDS(here("Rdata", "model_eg2.0.rds"))
# model_eg2.1 <- readRDS(here("Rdata", "model_eg2.1.rds"))
# model_eg2.2 <- readRDS(here("Rdata", "model_eg2.2.rds"))
# model_eg2.3 <- readRDS(here("Rdata", "model_eg2.3.rds"))

# Load required packages and dataset ----
pacman::p_load(
  dplyr, tibble, tidyverse, broom, broom.mixed,
  ape, arm, brms, broom.mixed, cmdstanr, emmeans, glmmTMB, MASS, phytools, rstan, TreeTools,
  DHARMa, loo, MuMIn, parallel, ggeffects,
  bayesplot, ggplot2, patchwork, tidybayes,
  gt, here, kableExtra, knitr
)

load(here("repro_package", "all_data.RData"))
dat_pref <- all_data$dat_pref
tree <- all_data$pref_tree
tree　<-force.ultrametric(tree, method = "extend")

# Run models ----

model_0 <- glmmTMB(
  frequency ~ 1 +  condition + sex + (1 | species) + (1 | id), # location part (mean)
  data = dat_pref,
  family = nbinom2(link = "log")
)
model0_res <- simulateResiduals(model_0, plot = TRUE, seed = 42)

# formal test for over/underdispersion
testDispersion(model0_res) 

## Baysian model (location-only) ---- 
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


## Location-scale models ----
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
