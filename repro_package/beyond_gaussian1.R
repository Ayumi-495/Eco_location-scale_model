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

dat_pref <- dat_pref %>%
  dplyr::select(-stripe, -full, -subset) %>%
  rename(frequency = dot) %>%
  mutate(species = phylo, across(c(condition, sex), as.factor))

# Visualise dataset ----
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

# Run models ----
## Location-only model ----
model_0 <- glmmTMB(
  frequency ~ 1 +  condition + sex + (1 | species) + (1 | id), # location part (mean)
  data = dat_pref,
  family = nbinom2(link = "log")
)

summary(model_0)
confint(model_0)

model0_res <- simulateResiduals(model_0, plot = TRUE, seed = 42)
testDispersion(model0_res) 

## Location-scale model ----
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

## Model comparison ----
model.sel(model_0, model_1)
model.sel(model_0, model_1, model_2)
