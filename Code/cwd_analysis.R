## Header ----
## Script name: 
##
## Purpose of script:
##
## Author: Nick Marzolf
## Date Created: 2021-12-07
## Date Modified: 
## Email: nmarzol@ncsu.edu
##
## load packages:  
{library(tidyverse)
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(car)
  library(drc)
  library(nlme)
  library(nlstools)
  library(qpcR)
  library(ggeffects)
}
##
## clear the environment if needed
rm(list = ls())
##
## set the ggplot theme
source("C:/Users/Nick Marzolf/Desktop/Research/R code/theme_nick.R")
theme_set(theme_nick())

# stream chemistry ----
# import
chem <- read_excel('Data/LT_CWD_datasheet.xlsx',
                   sheet = 'Chemistry') %>%
  rename(site = 'Site',
         srp = `SRP (ug/L)`,
         no3_n = `NO3-N (ug/L)`,
         nh4_n = `NH4-N (ug/L)`)

# mean of all measurmements
chem_sum <- chem %>% 
  group_by(site) %>%
  summarise(across(.cols = 2:7,
                   .fns = mean, na.rm = TRUE))

# calculate SD from stream chemistry data
chem_sd <- chem %>% 
  group_by(site) %>%
  summarise(across(.cols = 2:7,
                   .fns = sd, na.rm = TRUE))

# create object that sorts sites by decreasing mean conductivity
sites <- chem_sum %>%
  group_by(site) %>%
  summarise(mean_cond = mean(Cond, na.rm = TRUE)) %>%
  arrange(desc(mean_cond)) 

sites <- as.character(sites$site)

sites_long <- c(`Arb` = 'Arboleda 30',
                `Sur30` = 'Sura 30',
                `Tito60` = 'Saltito 60',
                `Piper` = 'Piper',
                `Tac` = 'Taconazo 30')


chem_sd$site <- fct_relevel(chem_sd$site,
                            sites)
chem_sd <- arrange(chem_sd, desc(Cond))

# fill in NA from Arb and Tac with long-term means
# the monthly Arb and Tac data from those sites isn't available yet
chem_sum[1,2] = 201.8
chem_sum[1,3] = 200.0
chem_sum[1,4] = 31.3

chem_sum[4,2] = 3.92
chem_sum[4,3] = 175.4
chem_sum[4,4] = 37.9

# calculate DIN and N:P ratio
chem_sum <- chem_sum %>%
  mutate(din = no3_n + nh4_n,
         n_p = (din/14.0067)/(srp/30.973762))


# relevel the site factor
chem_sum$site <- fct_relevel(chem_sum$site,
                             sites)

chem_sum <- arrange(chem_sum, desc(Cond))


# wood decomp rates ----

# load data
cwd <- read_excel('Data/LT_CWD_datasheet.xlsx',
                  sheet = 'Sheet1') %>% 
  filter(Flag == 0) %>%
  dplyr::select(site = Site, 
                month = `Collection Month`,
                rep = Rep,
                init_mass = `initial CWD mass (g)`,
                dry_mass = `CWD Pack Dry Mass (g)`,
                init_den = `init wood density (g/cm3)`,
                fin_den = `final wood density (g/cm3)`) %>%
  filter(month < 24)

cwd$site <- forcats::fct_relevel(cwd$site, sites)


cwd_calc <- cwd %>%
  dplyr::filter(site != 'Sac') %>%
  dplyr::mutate(percent_mass = (dry_mass/init_mass)*100) # same as for density

# final dry mass to export for bug analysis
final_dry_mass <- cwd_calc %>% 
  dplyr::select(site, month, rep, dry_mass) %>% 
  dplyr::filter(month != 0)

readr::write_csv(final_dry_mass,
                 'data/final_dry_mass.csv')

ggplot(cwd_calc, 
       aes(x = month, y = percent_mass,
           color = rep))+
  geom_point()+
  geom_smooth(method = 'lm')+
  facet_wrap(site ~ ., ncol = 5)

# # summarise by site and date
# cwd_sum <- cwd_calc %>%
#   group_by(site, month) %>%
#   summarise(mean_mass = mean(dry_mass, na.rm = TRUE),
#             sd_mass = sd(dry_mass, na.rm = TRUE),
#             mean_per_mass = mean(percent_mass, na.rm = TRUE),
#             sd_per_mass = sd(percent_mass, na.rm = TRUE),
#             mean_den = mean(fin_den, na.rm = TRUE),
#             sd_den = sd(fin_den, na.rm = TRUE)) 
# 
# # plot- percent change in mass (same as density)
# # yes, the Piper increases but it's better than the other options!
# ggplot(cwd_sum,
#        aes(x = month, 
#            y = mean_per_mass,
#            color = site))+
#   geom_point(size = 5)+
#   geom_smooth(method = 'lm',
#               alpha = 0.3)+
#   geom_errorbar(aes(ymin = mean_per_mass - sd_per_mass,
#                     ymax = mean_per_mass + sd_per_mass))+
#   facet_grid(. ~ site, 
#              #ncol = 5,
#              labeller = as_labeller(sites_long))+
#   scale_y_log10()+
#   labs(x = 'Months in Stream',
#        y = '% Dry Mass Remaining')+
#   scale_color_viridis_d()+
#   theme(legend.position = 'none')
# 
# 
# ancova <- aov(data = cwd_sum, 
#               mean_per_mass ~ month + site)
# summary(ancova)
# Anova(ancova, type = 'III')
# 
# tukey <- HSD.test(ancova, 'site')
# tukey


# calculate decay rates

k_cwd_add <- lm(data = cwd_calc,
                log(percent_mass) ~ month + site)
k_cwd_int <- lm(data = cwd_calc,
                log(percent_mass) ~ month * site)

anova(k_cwd_add,
      k_cwd_int) # use the interaction model (p = 0.02)

summary(k_cwd_int)
anova(k_cwd_int)
Anova(k_cwd_int, type = 'III')

# Because there is a significant interaction term, we can run a lm for each site
# to extract the slope and intercept
decay_rates <- cwd_calc %>%
  filter(site != 'Sac') %>%
  group_by(site) %>%
  summarise(
    int = round((coef(lm(log(percent_mass) ~ month))[1]), 3),
    k_yr = (round((coef(lm(log(percent_mass) ~ month))[2])*12, 3))*-1,
    error = summary(lm(log(percent_mass) ~ month))$coefficient[3],
    #df = summary(lm(log(percent_mass) ~ month))$fstatistic,
    r2 = summary(lm(log(percent_mass) ~ month))$r.squared,
    p = anova(lm(log(percent_mass) ~ month))$'Pr(>F'[1]
  )
decay_rates

ggplot(data = cwd_calc,
       aes(y = log(percent_mass),
           x = month))+
  geom_point(aes(color = site))+
  geom_smooth(aes(color = site),
              method = 'lm',se = FALSE)+
  geom_label(data = decay_rates %>% 
               dplyr::select(site, k_yr),
             aes(x = 2.5, y = 3.6,
                 label = paste0('k = ',k_yr)))+
  facet_grid(. ~ site, 
             labeller = as_labeller(sites_long))+
  labs(x = 'Months',
       y = 'ln(% DM Remaining)')+
  scale_color_viridis_d()+
  scale_x_continuous(breaks = c(0, 1, 3, 6, 12))+
  theme(legend.position = 'none')

# merge chemistry with decay rates ----

# extract k from decay_rates
merged <- decay_rates %>%
  dplyr::select(site, k_yr, error) %>% 
  right_join(chem_sum, 'site') %>%
  dplyr::select(-Temp) %>%
  pivot_longer(srp:n_p)

names_long <- c(`Cond` = 'Conductivity (µS/cm)',
                `din` = 'DIN (µg/L)',
                `n_p` = 'DIN:SRP',
                `nh4_n` = 'NH4-N (µg/L)',
                `no3_n` = 'NO3-N (µg/L)',
                `pH` = 'pH',
                `srp` = 'SRP (µg/L)')


# decay rates as a function of stream chemistry ----

vars <- unique(merged$name)

# linear: Y ~ a + bX
linear <- data.frame(mod = character(),
                     var = character(),
                     a = numeric(),
                     a_2.5 = numeric(),
                     a_97.5 = numeric(),
                     b = numeric(),
                     b_2.5 = numeric(),
                     b_97.5 = numeric(),
                     AIC = numeric(),
                     p = numeric(),
                     rse = numeric(),
                     df = numeric()
)

for(i in 1:length(vars)){
  var <- vars[i]
  lm <- lm(data = merged %>%
             filter(name == var),
           k_yr ~ value)
  
  print(plot(ggpredict(lm)))
  
  linear <- linear %>% 
    add_row(
      mod = 'Linear',
      var = var,
      a = coef(lm)[1],
      a_2.5 = confint(lm)[1,1],
      a_97.5 = confint(lm)[1,2],
      b = coef(lm)[2],
      b_2.5 = confint(lm)[2,1],
      b_97.5 = confint(lm)[2,2],
      AIC = AIC(lm),
      p = summary(lm)$coefficients[2,4],
      rse = summary(lm)$sigma,
      df = summary(lm)$df[2]
    )
}

linear %>% 
  filter(p < 0.05)



# non-linear
## Michaelis-Menten: Y ~ aX/b+X
micmen <- data.frame(mod = character(),
                     fit = character(),
                     var = character(),
                     a = numeric(),
                     a_2.5 = numeric(),
                     a_97.5 = numeric(),
                     b = numeric(),
                     b_2.5 = numeric(),
                     b_97.5 = numeric(),
                     AIC = numeric(),
                     p = numeric()
)

for(i in 1:length(vars)) {
  var <- vars[i]
  mm_mod <- try(nls(data = merged %>% 
                      filter(name == var),
                    k_yr ~ SSmicmen(value, a, b))
  )
  
  if(inherits(mm_mod, 'try-error')){
    micmen <- micmen %>% 
      add_row(fit = 'error',
              var = var,
              mod = 'Michaelis-Menten')
    next
  }
  
  micmen <- micmen %>% 
    add_row(
      mod = 'Michaelis-Menten',
      fit = 'success',
      var = var,
      a = coef(mm_mod)[1],
      a_2.5 = confint2(mm_mod, level = 0.95)['a',1],
      a_97.5 = confint2(mm_mod, level = 0.95)['a',2],
      b = coef(mm_mod)[2],
      b_2.5 = confint2(mm_mod, level = 0.95)['b',1],
      b_97.5 = confint2(mm_mod, level = 0.95)['b',2],
      AIC = AIC(mm_mod),
      p = summary(mm_mod)$coefficients[8]
    )
}

micmen %>% 
  filter(p < 0.05)

## Logarithmic: Y ~ a + b*log(X)
logarithmic <- data.frame(mod = character(),
                          var = character(),
                          a = numeric(),
                          a_2.5 = numeric(),
                          a_97.5 = numeric(),
                          b = numeric(),
                          b_2.5 = numeric(),
                          b_97.5 = numeric(),
                          AIC = numeric(),
                          R2 = numeric(),
                          p = numeric(),
                          rse = numeric(),
                          df = numeric()
)

for(i in 1:length(vars)) {
  var <- vars[i]
  
  loga <- lm(data = merged %>% 
               filter(name == var),
             k_yr ~ log10(value)) 
  
  logarithmic <- logarithmic %>% 
    add_row(
      mod = 'Logarithmic',
      var = var,
      a = coef(loga)[1],
      a_2.5 = confint(loga)[1],
      a_97.5 = confint(loga)[3],
      b = coef(loga)[2],
      b_2.5 = confint(loga)[2],
      b_97.5 = confint(loga)[4],
      AIC = AIC(loga),
      R2 = summary(loga)$r.squared,
      p = summary(loga)$coefficients[8],
      rse = summary(loga)$sigma,
      df = summary(loga)$df[2]
    )
}

logarithmic %>% 
  filter(p < 0.05)


## Logistic: Y ~ 1/1+exp(X)
logistic <- data.frame(mod = character(),
                       var = character(),
                       b = numeric(),
                       b_2.5 = numeric(),
                       b_97.5 = numeric(),
                       d = numeric(),
                       d_2.5 = numeric(),
                       d_97.5 = numeric(),
                       e = numeric(),
                       e_2.5 = numeric(),
                       e_97.5 = numeric(),
                       AIC = numeric(),
                       rse = numeric(),
                       df = numeric())

for(i in 1:length(vars)){
  var <- vars[i]
  
  logi <- drm(k_yr ~ value, 
              data = merged %>% 
                filter(name == var),
              fct = L.3())
  
  logistic <- logistic %>% 
    add_row(
      mod = 'Logistic',
      var = var,
      b = coef(logi)[1],
      b_2.5 = confint2(logi)['b:(Intercept)',1],
      b_97.5 = confint2(logi)['b:(Intercept)',2],
      d = coef(logi)[2],
      d_2.5 = confint2(logi)['d:(Intercept)',1],
      d_97.5 = confint2(logi)['d:(Intercept)',2],
      e = coef(logi)[3],
      e_2.5 = confint2(logi)['e:(Intercept)',1],
      e_97.5 = confint2(logi)['b:(Intercept)',2],
      AIC = AIC(logi),
      rse = summary(logi)$rseMat[1],
      df = summary(logi)$rseMat[2])
}


# combine model outputs

aic_all <- rbind(logistic %>% 
                   dplyr::select(mod, var, AIC) %>% 
                   mutate(p = NA),
                 logarithmic %>% 
                   dplyr::select(mod, var, AIC, p),
                 linear %>% 
                   dplyr::select(mod, var, AIC, p),
                 micmen %>% 
                   dplyr::select(mod, var, AIC, p))

best_fits <- aic_all %>% 
  arrange(var,
          AIC) %>% 
  group_by(var) %>% 
  dplyr::slice_head(n=1)


aic_wts <- data.frame()
for(i in 1:length(vars)) {
  use <- vars[i]
  
  df <- aic_all %>% 
    filter(!is.na(AIC),
           var %in% use) %>% 
    mutate(delAIC = akaike.weights(AIC)$deltaAIC,
           weights = akaike.weights(AIC)$weights)
  
  aic_wts <- rbind(aic_wts, df)
}


aic_wts %>% 
  dplyr::select(var, mod, weights) %>% 
  pivot_wider(names_from = var, 
              values_from = weights) %>% 
  knitr::kable(format = 'rst')


# plot decay rates as a function of chemistry
# add lines for SRP (linear), NH4 (M-M), N:P(M-M), and cond (linear)

merged_preds <- rbind(merged %>% 
                        filter(name == 'srp') %>% 
                        mutate(pred_k = -0.12456 + (-0.002889*value)),
                      merged %>% 
                        filter(name == 'Cond') %>% 
                        mutate(pred_k = -0.075 + (-0.0022*value)) ,
                      merged %>% 
                        filter(name == 'nh4_n') %>% 
                        mutate(pred_k = (-0.0768*value)/(-28.13+value)) ,
                      merged %>% 
                        filter(name == 'n_p') %>% 
                        mutate(pred_k = (-0.179*value)/(-1.944+value)),
                      merged %>% 
                        filter(name %in% c('no3_n', 'din', 'pH')) %>% 
                        mutate(pred_k = NA))

merged_preds <- data.frame()
for(i in 1:length(best_fits$var)) {
  fit <- best_fits[i,]
  
  type <- pull(fit[1])
  var <- pull(fit[2])
  
  if(type == 'Logistic'){
    logis <- logistic %>% 
      filter(var == !!var)
    
    df <- merged %>% 
      filter(name == var) %>% 
      mutate(pred_k = logis$d/(1+exp(logis$b*(value - logis$e))),
             best_fit = type)
  
  } # end logistic if statement
  
  if(type == 'Linear') {
    line <- linear %>% 
      filter(var == !!var)
    
    if(line$p <= 0.05){
    df <- merged %>% 
      filter(name == var) %>% 
      mutate(pred_k = line$a + (value*line$b),
             best_fit = type)
    } else {
      df <- merged %>% 
        filter(name == var) %>% 
        mutate(pred_k = NA,
               best_fit = NA)
    }
  } # end linear if statement
  
  if(type == 'Michaelis-Menten'){
    mm <- micmen %>% 
      filter(var == !!var)
    
    if(mm$p <= 0.05){
    df <- merged %>% 
      filter(name == var) %>% 
      mutate(pred_k = (mm$a*value)/(value + mm$b),
             best_fit = type)
    } else {
      df <- merged %>% 
        filter(name == var) %>% 
        mutate(pred_k = NA,
               best_fit = NA)
    }
  } # end M-M for loop
  
  if(type == 'Logarithmic'){
    logar <- logarithmic %>% 
      filter(var == !!var)
    
    if(logar$p <= 0.05) {
      df <- merged %>% 
        filter(name == var) %>% 
        mutate(pred_k = logar$a + (logar$b*value),
               best_fit = type)
    } else {
      df <- merged %>% 
        filter(name == var) %>% 
        mutate(pred_k = NA,
               best_fit = NA)
    }
  } # end logarithmic if statement
  
  merged_preds <- rbind(merged_preds, df)
} # end for loop


fig3 <- ggplot(merged_preds)+
  geom_point(aes(x = value,
                 y = k_yr,
                 color = site),
             size = 4)+
  geom_line(aes(x = value, 
                y = pred_k),
            size = 1,
            linetype = 'dashed')+
  ylab(expression(paste(k[CWD], ' (', y^-1,')')))+
  facet_wrap(. ~ name,
             nrow = 2,
             scales = 'free',
             labeller = as_labeller(names_long),
             drop = FALSE)+
  scale_color_viridis_d(name = 'Site',
                        labels = as_labeller(sites_long))+
  theme(axis.title.x = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

shift_legend2 <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"
  
  # now we just need a simple call to reposition the legend
  reposition_legend(p, 'center', panel=names)
}

shift_legend2(fig3)


# old plots ----
# plot- raw change in density (VALUES INCREASE, DON'T USE)
ggplot(cwd_sum,
       aes(x = month, 
           y = mean_den,
           color = site))+
  geom_point(size = 5)+
  geom_smooth(method = 'lm',
              alpha = 0.3)+
  geom_errorbar(aes(ymin = mean_den - sd_den,
                    ymax = mean_den + sd_den))+
  facet_wrap(site ~ ., 
             ncol = 5,
             labeller = as_labeller(sites_long))+
  labs(x = 'Months in Stream')+
  scale_color_viridis_d()+
  theme(legend.position = 'none')


# plot- raw dry mass (VALUES INCREASE, DON'T USE)
ggplot(cwd_sum,
       aes(x = month, 
           y = mean_mass,
           color = site))+
  geom_point(size = 5)+
  geom_smooth(method = 'lm',
              alpha = 0.3)+
  geom_errorbar(aes(ymin = mean_mass - sd_mass,
                    ymax = mean_mass + sd_mass))+
  facet_wrap(site ~ ., 
             ncol = 5,
             labeller = as_labeller(sites_long))+
  labs(x = 'Months in Stream')+
  scale_color_viridis_d()+
  theme(legend.position = 'none')