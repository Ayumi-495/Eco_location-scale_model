dat <- dat %>%
  mutate(
    if_kill = factor(if_kill,
                     levels = c("burro kills absent", "burro kills present")),
    Pool    = as.factor(Pool) 
  )

nd <- tidyr::expand_grid(
  if_kill = levels(dat$if_kill),
  Pool = levels(dat$Pool)
) %>%
  dplyr::mutate(
    if_kill = factor(if_kill, levels = levels(dat$if_kill)),
    Pool    = factor(Pool, levels = levels(dat$Pool))
  )

# --- Get linear predictors (link scale) for all four submodels at once ---
# .linpred = linear predictor for mu (here logit(mu))
# phi, zoi, and coi are also returned on the link scale:
#   phi → log(phi), zoi → logit(zoi), coi → logit(coi)
lp_draws <- linpred_draws(
  fit1, newdata = nd, re_formula = NULL, dpar = TRUE, transform = FALSE
)

# Reshape and summarise posterior draws (with enforced facet order) ----
df_lp <- lp_draws %>%
  dplyr::select(if_kill, .draw, .linpred, phi, zoi, coi) %>%
  dplyr::rename(mu = .linpred) %>%
  tidyr::pivot_longer(mu:coi, names_to = "param", values_to = "value") %>%
  dplyr::mutate(
    param = factor(param, levels = c("mu", "phi", "zoi", "coi"))
  ) %>%
  dplyr::group_by(if_kill, param) %>%
  mean_qi(value, .width = 0.95)

df_lp
# A tibble: 8 × 8
# if_kill             param  value   .lower .upper .width .point .interval
# <fct>               <fct>  <dbl>    <dbl>  <dbl>  <dbl> <chr>  <chr>    
# 1 burro kills absent  mu     0.915 -0.509   2.37    0.95 mean   qi       
# 2 burro kills absent  phi    1.73   0.723   2.61    0.95 mean   qi       
# 3 burro kills absent  zoi    0.724  0.0984  1.39    0.95 mean   qi       
# 4 burro kills absent  coi    6.28   2.61   14.1     0.95 mean   qi       
# 5 burro kills present mu    -0.303 -1.61    0.922   0.95 mean   qi       
# 6 burro kills present phi    0.664  0.313   1.00    0.95 mean   qi       
# 7 burro kills present zoi   -0.985 -1.53   -0.485   0.95 mean   qi       
# 8 burro kills present coi    2.08   0.845   3.65    0.95 mean   qi       


# More descriptive y-axis labels for each parameter ----
lab_param <- c(
  mu  = "Mean of Beta component (logit(mu))",
  phi = "Precision of Beta component (log(phi))",
  zoi = "Zero-inflation probability (logit(zoi))",
  coi = "One-inflation probability (logit(coi))"
  )

# Plot linear predictors with 95% credible intervals ----
p_params <- ggplot(df_lp, aes(x = if_kill, 
                              y = value,
                              ymin = .lower, 
                              ymax = .upper, 
                              color = if_kill)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.15) +
  scale_color_manual(values = c("steelblue4", "orange3")) +
  labs(x = "Burro kill present?") +
  facet_wrap(~ param, scales = "free_y",
             labeller = ggplot2::labeller(param = lab_param)) +
  theme_classic() +
  theme(legend.position = "none")

# Predicted mean on the response scale (overall expected cover) ----
# df_ep <- newdat_pop %>% add_epred_draws(fit1, re_formula = NA) %>%
#   dplyr::group_by(if_kill) %>% median_qi(.epred)
# 
# p_ep <- ggplot(df_ep, aes(if_kill, 
#                           .epred, 
#                           ymin = .lower, 
#                           ymax = .upper, 
#                           color = if_kill)) +
#   geom_point(size = 2) +
#   geom_errorbar(width = 0.15) +
#   scale_color_manual(values = c("steelblue4", "orange3")) +
#   labs(x = "Burro kill present?", y = "Expected cover (0–1)") +
#   theme_classic() +
#   theme(legend.position = "none")

# Combined ----
p_params
# / p_ep

### predicted values ----

pre <- as.data.frame(predict(fit1))
predict(fit1)

pp_check(fit1, ndraws = 30)

dat <- dat %>%
  mutate(
    if_kill = factor(if_kill,
                     levels = c("burro kills absent", "burro kills present")),
    Pool    = as.factor(Pool) 
  )

nd_pool <- tidyr::expand_grid(
  if_kill = levels(dat$if_kill),
  Pool = levels(dat$Pool)
) %>%
  dplyr::mutate(
    if_kill = factor(if_kill, levels = levels(dat$if_kill)),
    Pool    = factor(Pool, levels = levels(dat$Pool))
  )

df_ppc_RE <- nd_pool %>%
  tidybayes::add_predicted_draws(fit1, re_formula = NULL) %>%   
  dplyr::group_by(.draw, if_kill) %>%                           
  dplyr::summarise(pred_mean = mean(.prediction), .groups = "drop") %>%  
  tidybayes::mean_qi(pred_mean, .width = 0.95)

p_ppc_RE <- ggplot(df_ppc_RE,
                   aes(x = if_kill, y = pred_mean, ymin = .lower, ymax = .upper, color = if_kill)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.15) +
  scale_color_manual(values = c("steelblue4","orange3")) +
  labs(x = "Burro kill present?", y = "Posterior predictive (0–1)\n(RE included, pooled over Pools)") +
  theme_classic() +
  theme(legend.position = "none")

