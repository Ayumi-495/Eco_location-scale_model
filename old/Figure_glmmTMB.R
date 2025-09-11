# library(ggeffects)
library(emmeans)

# model 1----
# Extract estimated marginal means for location part (conditional model)
# returns predicted log(AdTarsus) by Sex*Treatment
emm_mu_link <- emmeans(model_1, ~ Sex * Treatment,
                       component = "cond", type = "link")
df_mu_link <- as.data.frame(emm_mu_link) %>%
  rename(
    mean_log = emmean, # estimated log(mean AdTarsus)
    lwr = lower.CL, # lower 95% CI (log scale)
    upr = upper.CL # upper 95% CI (log scale)
  )

# Extract estimated marginal means for scale part (dispersion model)
emm_sigma_link <- emmeans(model_1, ~ Sex * Treatment,
                          component = "disp", type = "link")
df_sigma_link <- as.data.frame(emm_sigma_link) %>%
  rename(
    log_sigma = emmean, # estimated log(sigma)
    lwr = lower.CL, # lower 95% CI (log scale)
    upr = upper.CL # upper 95% CI (log scale)
  )

pos <- position_dodge(width = 0.4)

# Visualisation (location part)
# Plot predicted log(AdTarsus) ± 95% CI for Sex*Treatment
p_mean_link <- ggplot(df_mu_link, aes(x = Sex, y = mean_log, color = Treatment)) +
  geom_point(position = pos, size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), position = pos, width = 0.15) +
  labs(
    title = "Location",
    x = "Sex",
    y = " log(AdTarsus)"
  ) +
  scale_y_continuous(limits = c(2.5, 3.0), breaks = seq(2.5, 3.0, by = 0.1) ) + 
  theme_classic()


# Visualisation (scale part)
# Plot predicted log(sigma) ± 95% CI for Sex*Treatment
p_sigma_link <- ggplot(df_sigma_link, aes(x = Sex, y = log_sigma, color = Treatment)) +
  geom_point(position = pos, size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), position = pos, width = 0.15) +
  labs(
    title = "Scale",
    x = "Sex",
    y = "log(sigma)"
  ) +
  scale_y_continuous(limits = c(-4.5, -2.5), breaks = seq(-4.5, -2.5, by = 0.1) ) + 
  theme_classic()

# scale_data <- log(sqrt(residuals(model_1)^2)) # this is the residual  - we can add the figure, but it does not come from the model
     
p_glmmTMB_m1 <- p_mean_link + p_sigma_link
model_1
# beyond gaussian1 ----

# Location part (conditional model): log(mean frequency)
emm_mu_link <- emmeans(model_pref, ~ condition * sex,
                       component = "cond", type = "link")  # return log(mu)

head(as.data.frame(emm_mu_link)) 
# condition sex   emmean        SE  df asymp.LCL asymp.UCL
# deprived  F   4.873695 0.2041941 Inf  4.473482  5.273908
# supplied  F   4.028497 0.2231513 Inf  3.591129  4.465866
# deprived  M   4.769738 0.2079726 Inf  4.362120  5.177357
# supplied  M   3.924540 0.2228790 Inf  3.487705  4.361375
# 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 

df_mu_link <- as.data.frame(emm_mu_link) %>%
  rename(
    mean_log = emmean, # estimated log(mean frequency)
    lwr = asymp.LCL, # lower 95% CI on log scale
    upr = asymp.UCL # upper 95% CI on log scale
  )

# Scale part (dispersion model): log(theta)
emm_disp_link <- emmeans(model_pref, ~ condition * sex,
                         component = "disp", type = "link")  # return log(θ)
df_disp_link <- as.data.frame(emm_disp_link) %>%
  rename(
    log_theta = emmean, # estimated log(dispersion parameter theta)
    lwr = asymp.LCL, # lower 95% CI on log scale
    upr = asymp.UCL # upper 95% CI on log scale
  )

pos <- position_dodge(width = 0.4)

# Plot location (mean, on log scale)
p_mu_link <- ggplot(df_mu_link, aes(x = sex, y = mean_log, color = condition)) +
  geom_point(position = pos, size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), position = pos, width = 0.15) +
  labs(
    title = "Location",
    x = "Sex", y = "log(mean frequency)" # log scale
  ) +
  scale_y_continuous(limits = c(3.0, 6.0), breaks = seq(3.0, 6.0, by = 0.1) ) + 
  theme_classic()

p_mu_link_raw <- p_mu_link +
  geom_point(
    data = dat_pref,
    aes(x = sex, y = log(frequency), color = condition),
    position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.08),
    alpha = 0.35, size = 1.6, inherit.aes = FALSE
  )


# Plot dispersion (theta, on log scale)
p_disp_link <- ggplot(df_disp_link, aes(x = sex, y = log_theta, color = condition)) +
  geom_point(position = pos, size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), position = pos, width = 0.15) +
  labs(
    title = "Scale",
    x = "Sex", y = "log(theta)" # log scale
  ) +
  scale_y_continuous(limits = c(-0.5, 2.0), breaks = seq(-0.5, 2.0, by = 0.5) ) + 
  theme_classic()

p_glmmTMB_nb <- p_mu_link_raw + p_disp_link

