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

# Run models ----
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