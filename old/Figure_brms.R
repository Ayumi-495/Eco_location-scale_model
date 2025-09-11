
library(dplyr)
library(ggplot2)
library(tidybayes)
library(patchwork)

# model 1 ----
nd <- dat_tarsus %>%
  distinct(Sex, Treatment) %>%
  arrange(Sex, Treatment)

draws_mu <- epred_draws(brms_g1, newdata = nd, re_formula = NA) %>%
  median_qi(.epred, .width = 0.95)

draws_sig <- epred_draws(brms_g1, newdata = nd, re_formula = NA, dpar = TRUE) %>%
  transmute(Sex, Treatment, log_sigma = log(sigma)) %>%
  median_qi(log_sigma, .width = 0.95)

dat_points_loc <- dat_tarsus %>%
  mutate(y_log = log(AdTarsus))

lp_sigma <- posterior_linpred(
  brms_g1, newdata = dat_tarsus, re_formula = NA,
  dpar = "sigma", transform = FALSE
)
dat_points_scl <- dat_tarsus %>%
  mutate(log_sigma_hat = apply(lp_sigma, 2, median))

pos <- position_dodge(width = 0.5)

p_loc <- ggplot(draws_mu, aes(x = Sex, y = .epred, color = Treatment)) +
  geom_point(position = pos, size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), position = pos, width = 0.15) +
  geom_point(
    data = dat_points_loc,
    aes(x = Sex, y = y_log, color = Treatment),
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.08),
    size = 1.6, alpha = 0.35, inherit.aes = FALSE
  ) +
  labs(title = "Location",
       x = "Sex", y = "log(AdTarsus)"
       ) +
  scale_y_continuous(limits = c(2.5, 3.0), breaks = seq(2.5, 3.0, by = 0.1) ) +
  theme_classic()

p_scl <- ggplot(draws_sig, aes(x = Sex, y = log_sigma, color = Treatment)) +
  geom_point(position = pos, size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), position = pos, width = 0.15) +
  geom_point(
    data = dat_points_scl,
    aes(x = Sex, y = log_sigma_hat, color = Treatment),
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.08),
    size = 1.6, alpha = 0.35, inherit.aes = FALSE
  ) +
  labs(title = "Scale",
       x = "Sex", y = "log(sigma)"
       ) +
  scale_y_continuous(limits = c(-4.5, -2.5), breaks = seq(-4.5, -2.5, by = 0.1) ) + 
  theme_classic()

p_brms_m1 <- p_loc + p_scl


# beyond gausian 1 ----
library(tidybayes) # add_linpred_draws() for posterior draws on the link scale

dat_pref <- dat_pref %>%
  mutate(
    condition = factor(condition),
    sex = factor(sex)
  )

# Newdata grid for all Sex and Condition
nd <- tidyr::expand_grid(
  condition = levels(dat_pref$condition),
  sex  = levels(dat_pref$sex)
)

# Location(mu): log(mean) on the model scale
# add_linpred_draws(..., dpar = "mu") returns draws of the linear predictor (log link) for mu
draws_mu <- add_linpred_draws(
  nd, 
  object = model_nb_brms, 
  dpar   = "mu", # distributional parameter 'mu'
  re_formula = NA # population-level (exclude group-level REs) for EMM-like estimates
)

# Summarise posterior draws to median and 95% CrI on the log scale
sum_mu <- draws_mu %>%
  group_by(condition, sex) %>%
  summarise(
    mean_log = median(.linpred), # posterior median of log(mu)
    lwr      = quantile(.linpred, 0.025), # 95% CrI lower (log scale)
    upr      = quantile(.linpred, 0.975), # 95% CrI upper (log scale)
    .groups = "drop"
  )

# Scale (shape theta): log(theta) on the model scale
# add_linpred_draws(..., dpar = "shape") returns draws for log(theta)
draws_shape <- add_linpred_draws(
  nd,
  object = model_nb_brms, 
  dpar   = "shape",
  re_formula = NA
)

sum_shape <- draws_shape %>%
  group_by(condition, sex) %>%
  summarise(
    log_theta = median(shape), # posterior median of log(theta)
    lwr       = quantile(shape, 0.025),
    upr       = quantile(shape, 0.975),
    .groups = "drop"
  )

# Plot: Location (log(mu))
pos <- position_dodge(width = 0.4)

p_loc <- ggplot(sum_mu, aes(x = sex, y = mean_log, color = condition)) +
  geom_point(position = pos, size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), position = pos, width = 0.15) +
  labs(
    title = "Location",
    x = "Sex", y = "log(mean frequency)"
  ) +
  scale_y_continuous(limits = c(3.0, 6.0), breaks = seq(3.0, 6.0, by = 0.5) ) +
  theme_classic()

p_mean_link_raw <- p_loc +
  geom_point(
    data = dat_pref,
    aes(x = sex, y = log(frequency), color = condition),
    position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.08),
    alpha = 0.35, size = 1.6, inherit.aes = FALSE
  )


# Plot: Scale (log(theta)) 
p_shape <- ggplot(sum_shape, aes(x = sex, y = log_theta, color = condition)) +
  geom_point(position = pos, size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), position = pos, width = 0.15) +
  labs(
    title = "Scale",
    x = "Sex", y = "log(theta)"
  ) +
  scale_y_continuous(limits = c(-0.5, 2.0), breaks = seq(-0.5, 2.0, by = 0.5) ) + 
  theme_classic()

p_brms_nb <- p_mean_link_raw + p_shape


# All ----

## model 1 ----
(p_glmmTMB_m1 / p_brms_m1) +
  plot_annotation(title = "top: glmmTMB / bottom: brms")

## nb model ----
(p_glmmTMB_nb / p_brms_nb) +
plot_annotation(
  title = "top: glmmTMB / bottom: brms"
)

### ### ###
## assume m is a fitted glmmTMB model with family nbinom2 or nbinom1
## and data frame d used to fit it

## 1) fixed-effect part of the dispersion linear predictor
X_disp   <- model.matrix(model_pref, component = "disp")        # design matrix for disp
b_disp   <- fixef(model_pref)$disp                               # log-dispersion coefs
eta_fix  <- as.vector(X_disp %*% b_disp)                # fixed-effects contribution

## 2) (optional) random-effect contribution in the dispersion model
eta_re <- 0
re_list <- try(ranef(model_pref)$disp, silent = TRUE)            # BLUPs for disp REs
if (!inherits(re_list, "try-error") && length(re_list)) {
  ## common cases: (1|g) and (x|g) in the dispformula
  for (gf in names(re_list)) {
    RE <- re_list[[gf]]                                  # data.frame of RE per level
    ii <- match(d[[gf]], rownames(RE))                   # map obs -> group row
    if ("(Intercept)" %in% names(RE)) {
      eta_re <- eta_re + RE[ii, "(Intercept)"]
    }
    ## add random slopes present in the disp RE structure
    slope_names <- setdiff(colnames(RE), "(Intercept)")
    for (sn in slope_names) {
      eta_re <- eta_re + RE[ii, sn] * d[[sn]]
    }
  }
}

## 3) observation-level dispersion on the natural scale
eta_disp <- eta_fix + eta_re
disp_par <- exp(eta_disp)
hist(disp_par)
## 4) interpret according to the family
fam <- family(model_pref)$family
if (identical(fam, "nbinom2")) {
  theta_hat <- disp_par                     # θ_i
} else if (identical(fam, "nbinom1")) {
  alpha_hat <- disp_par                     # α_i
} else stop("This code is for NB1/NB2 only.")

