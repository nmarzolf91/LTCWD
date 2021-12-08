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
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readxl)
##
## clear the environment if needed
rm(list = ls())
##
## set the ggplot theme
source("C:/Users/nmarz/Desktop/Research/R code/theme_nick.R")
theme_set(theme_nick())

# wood decomp rates ----

# load data
dat <- read_excel('Data/LT_CWD_datasheet.xlsx',
           sheet = 'Sheet1') %>% 
  select(site = Site, 
         month = `Collection Month`,
         rep = Rep,
         init_mass = `initial CWD mass (g)`,
         dry_mass = `Dry Mass (g)`,
         ash_mass = `Sub Ash Mass (g)`,
         ash_frac = `Ash Fraction`,
         org_frac = `Organic Fraction`,
         pack_afdm = `Pack AFDM (g)`,
         afdm_remain = `AFDM Remaining`)

# summarise by site and date
dat_sum <- dat %>%
  group_by(site, month) %>%
  summarise(mean_afdm_remain = mean(afdm_remain, na.rm = TRUE),
            sd_afdm_remain = sd(afdm_remain, na.rm = TRUE),
            mean_afdm = mean(pack_afdm, na.rm = TRUE),
            sd_afdm = mean(pack_afdm, na.rm = TRUE),
            mean_dry = mean(dry_mass, na.rm = TRUE)) %>%
  filter(site != 'Sac')

# plot the dry mass in each bag at each collection point
# there are problems with the AFDM foil tins and are problematic
ggplot(dat_sum,
       aes(x = month, y = log(mean_dry),
           color = site))+
  geom_point()+
  geom_smooth(method = 'lm',
              alpha = 0.3)
  # geom_errorbar(aes(ymin = mean_afdm - sd_afdm,
  #                   ymax = mean_afdm + sd_afdm))+
  # facet_wrap(site ~ .,
  #            scales = 'free')


# calculate decay rates
decay_rates <- dat %>%
  filter(site != 'Sac') %>%
  group_by(site) %>%
  summarise(
    k = coef(lm(log(dry_mass) ~ month))[2],
    r2 = summary(lm(log(dry_mass) ~ month))$r.squared
    )

# stream chemistry ----
chem <- read_excel('Data/LT_CWD_datasheet.xlsx',
                   sheet = 'Chemistry') %>%
  rename(site = 'Site',
         srp = `SRP (ug/L)`)

chem_sum <- chem %>% 
  group_by(site) %>%
  summarise(across(.cols = 2:7,
                   .fns = mean, na.rm = TRUE))

# fill in NA from Arb and Tac with long-term means
# the monthly Arb and Tac data from those sites isn't available yet
chem_sum[1,2] = 201.8
chem_sum[1,3] = 200.0
chem_sum[1,4] = 31.3

chem_sum[4,2] = 3.92
chem_sum[4,3] = 175.4
chem_sum[4,4] = 37.9


merged <- left_join(decay_rates, chem_sum, 'site')

summary(lm(data = merged, 
           k ~ srp))

merged %>%
  ggplot(., aes(x = srp, y = -k*12,
                color = site))+
  geom_point(size = 5)+
  xlab(expression(paste('SRP (µg', " ", L^-1,")")))+
  ylab(expression(paste('Decay rate ( ', yr^-1,')')))


