## Header ----
## Script name: 
##
## Purpose of script:
##
## Author: Nick Marzolf
## Date Created: 2023-03-01
## Date Modified: 
## Email: nicholas.marzolf@duke.edu
##
## load packages:  
library(tidyverse)
library(dplyr)
library(ggplot2)
library(vegan)
library(ggrepel)
##
## clear the environment if needed
rm(list = ls())
##
## set the ggplot theme
source("C:/Users/Nick Marzolf/Desktop/Research/R code/theme_nick.R")
theme_set(theme_nick())



# 1) import data ----
data <- readxl::read_excel('Data/Samples_Nick.xlsx')
taxonomy <- readxl::read_excel('Data/Samples_Nick.xlsx',
                               sheet = 'taxonomy')

final_dry_mass <- readr::read_csv('data/final_dry_mass.csv')

# raw spreadsheet from Ana had a site 'Saltito 30'.
# Mar 1 2023: Change to Saltito 60, assuming this was a typo

unique(data$Stream) # should have 5 sites

# put in order of increasing conductivity
data$Stream <- fct_relevel(data$Stream,
                           c('Arboleda 30', 'Sura 30', 'Saltito 60',
                             'Piper', 'Taconazo 30'))

final_dry_mass$site <- recode_factor(final_dry_mass$site,
                                     Arb = 'Arboleda 30',
                                     Sur30 = 'Sura 30',
                                     Tito60 = 'Saltito 60',
                                     Piper = 'Piper',
                                     Tac = 'Taconazo 30')

# 2) add full taxonomy to each bug ID ----
tax_data <- list()
for(i in 1:length(unique(data$Taxa))){
  rank <- unique(data$Taxa)[i]
  
  tax_info <- taxonomy %>% 
    filter(Taxa %in% rank)
  
  tax_data[[i]] <- data %>% 
    data.frame() %>% 
    dplyr::filter(Taxa == rank) %>% 
    mutate(rank = tax_info$Rank,
           family = tax_info$Family,
           order = tax_info$Order,
           class = tax_info$Class,
           phylum = tax_info$Phylum,
           month = as.numeric(gsub(".*?([0-9]+).*", "\\1", Sample)),
           rep = gsub("\\d+", "", Sample))
}

# create df
tax_data <- reduce(tax_data, rbind)

# clean the data
tax_data_clean <- tax_data %>% 
  filter(QAQC == 0,
         month != 24) 


# summary stats
sum(tax_data_clean$Total)

tax_data_clean %>% 
  group_by(Stream) %>% 
  summarise(sum = sum(Total))

tax_data_clean %>% 
  group_by(order) %>% 
  summarise(sum = sum(Total))


# 3) NMDS's ----

# 3.1) All dates & sites, at Order level

# which orders are present
unique_orders <- na.omit(unique(tax_data_clean$order))

length(unique_orders)

tax_data_clean_order <- tax_data_clean %>% 
  dplyr::filter(!is.na(order)#,
                #order != 'Neoophora'
  ) %>%                               # remove NAs from Order
  dplyr::select(Stream, order, month, rep, Total) %>%    # get the necessary columns
  dplyr::group_by(Stream, month, order) %>%         # and do the grouping
  dplyr::summarise(mean_order = mean(Total, na.rm = TRUE)) %>%   # sum by order in each possible grouping 
  dplyr::mutate(log_total = log10(1 + mean_order)) %>%          # log10 + 1 transform data
  dplyr::select(-mean_order) %>% 
  pivot_wider(names_from = order,                         # pivot data
              values_from = log_total,
              values_fill = 0) 


tax_data_matrix <- tax_data_clean_order[,-c(1:2)]
tax_data_meta <- tax_data_clean_order[,c(1:2)]

nmds_order <- vegan::metaMDS(tax_data_matrix,
                             distance = 'bray',
                             k = 3,
                             autotransform = FALSE)
stress <- round(nmds_order$stress, 3)

nmds_order_out <- data.frame(x = nmds_order$points[,1],
                             y = nmds_order$points[,2])

nmds_order_out <- cbind(tax_data_meta,
                        nmds_order_out)

fit <- (vegan::envfit(nmds_order,
                      tax_data_matrix,
                      perm = 9999))

scrs <- data.frame(vegan::scores(fit, 'vectors'))

scrs$pvals <- fit$vectors$pvals

scrs_sig <- subset(scrs, pvals <= 0.05)

scrs_sig$env.variables <- row.names(scrs_sig)



plot_nmds_order_scrs <- ggplot(nmds_order_out,
                               aes(x = x, y = y))+
  geom_point(size = 2)+
  geom_segment(data = scrs_sig,
               aes(x = 0, xend = NMDS1,
                   y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, 'cm')),
               color = 'black')+
  ggrepel::geom_label_repel(data = scrs_sig,
                            aes(NMDS1, NMDS2,
                                label = env.variables),
                            size = 3)+
  labs(x = 'NMDS 1',
       y = 'NMDS 2')+
  lims(x = c(-1,1),
       y = c(-0.75, 0.75))+
  geom_label(label = paste('3D Stress = ', stress),
             x = 0.55,
             y = 0.75)
plot_nmds_order_scrs


source('Code/veganCovEllipse.R')
streams <- data.frame()
for(i in unique(nmds_order_out$Stream)){
  streams <- rbind(streams,
                   cbind(
                     as.data.frame(
                       with(nmds_order_out[nmds_order_out$Stream == i,],
                            veganCovEllipse(cov.wt(cbind(x, y),
                                                   wt = rep(1/length(x),
                                                            length(x)))$cov,
                                            center = c(mean(x),
                                                       mean(y)
                                            )
                            )
                       )
                     ), 
                     Stream = i)
  )
} # end for loop


plot_nmds_streams <- ggplot(data = nmds_order_out,
                            aes(x = x, y = y))+
  geom_point(aes(color = Stream),
             size = 2)+
  geom_path(data = streams,
            linewidth = 1,
            aes(x = x, y = y, color = Stream))+
  labs(x = "NMDS 1", 
       y = "NMDS 2")+
  lims(x = c(-1,1),
       y = c(-0.75,0.75))+
  theme(legend.background = element_blank())+
  scale_color_viridis_d(name = element_blank())
plot_nmds_streams

ggpubr::ggarrange(plot_nmds_order_scrs,
                  plot_nmds_streams,
                  ncol = 2,
                  align = 'h', widths = c(1,1.4),
                  labels = 'auto',
                  label.x = c(0.18,0.13))

# Month NMDS
months <- data.frame()
for(i in unique(nmds_order_out$month)){
  months <- rbind(months,
                  cbind(
                    as.data.frame(
                      with(nmds_order_out[nmds_order_out$month == i,],
                           veganCovEllipse(cov.wt(cbind(x, y),
                                                  wt = rep(1/length(x),
                                                           length(x)))$cov,
                                           center = c(mean(x),
                                                      mean(y)
                                           )
                           )
                      )
                    ), 
                    month = i)
  )
} # end for loop

plot_nmds_months <- ggplot(data = nmds_order_out,
                           aes(x = x, 
                               y = y))+
  geom_point(aes(color = factor(month)),
             size = 2)+
  geom_path(data = months,
            linewidth = 1,
            aes(x = x, 
                y = y, 
                color = factor(month)))+
  labs(x = "NMDS 1", 
       y = "NMDS 2")+
  lims(x = c(-1,1),
       y = c(-0.75,0.75))+
  theme(legend.background = element_blank())+
  scale_color_viridis_d(name = element_blank())
plot_nmds_months


# 3.2) NMDS at family level
unique_fams <- na.omit(unique(tax_data_clean$family))

length(unique_fams)

tax_data_clean_family <- tax_data_clean %>% 
  dplyr::filter(!is.na(family)) %>%                               # remove NAs from Order
  dplyr::select(Stream, family, month, rep, Total) %>%    # get the necessary columns
  dplyr::group_by(Stream, month, family) %>%         # and do the grouping
  dplyr::summarise(mean_family = mean(Total, na.rm = TRUE)) %>%   # sum by order in each possible grouping 
  dplyr::mutate(log_total = log10(1 + mean_family)) %>%          # log10 + 1 transform data
  dplyr::select(-mean_family) %>% 
  pivot_wider(names_from = family,                         # pivot data
              values_from = log_total,
              values_fill = 0)


tax_data_matrix_fam <- tax_data_clean_family[,-c(1:2)]
tax_data_meta_fam <- tax_data_clean_family[,c(1:2)]

nmds_family <- vegan::metaMDS(tax_data_matrix_fam,
                              distance = 'bray',
                              k = 3,
                              autotransform = FALSE)
stress_fam <- round(nmds_family$stress, 3)

nmds_family_out <- data.frame(x = nmds_family$points[,1],
                              y = nmds_family$points[,2])

nmds_family_out <- cbind(tax_data_meta_fam,
                         nmds_family_out)

fit_fam <- (vegan::envfit(nmds_family,
                      tax_data_matrix_fam,
                      perm = 9999))

scrs_fam <- data.frame(vegan::scores(fit_fam, 'vectors'))

scrs_fam$pvals <- fit_fam$vectors$pvals

scrs_fam_sig <- subset(scrs_fam, pvals <= 0.05)

scrs_fam_sig$env.variables <- row.names(scrs_fam_sig)


plot_nmds_family_scrs <- ggplot(nmds_family_out,
                               aes(x = x, y = y))+
  geom_point(size = 2)+
  geom_segment(data = scrs_fam_sig,
               aes(x = 0, xend = NMDS1,
                   y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, 'cm')),
               color = 'black')+
  ggrepel::geom_label_repel(data = scrs_fam_sig,
                            aes(NMDS1, NMDS2,
                                label = env.variables),
                            size = 3)+
  labs(x = 'NMDS 1',
       y = 'NMDS 2')+
  # lims(x = c(-1,1),
  #      y = c(-0.75, 0.75))+
  geom_label(label = paste('3D Stress = ', stress),
             x = -0.65,
             y = -0.75)
plot_nmds_family_scrs

# 4) bug data over time ----
macro_density <- tax_data_clean %>% 
  dplyr::filter(!is.na(order)) %>%                               # remove NAs from Order
  dplyr::select(site = Stream, order, month, rep, Total) %>%     # get the necessary columns
  dplyr::group_by(site, month, rep) %>%                          # and do the grouping
  dplyr::summarise(total_order = sum(Total)) %>% 
  dplyr::left_join(final_dry_mass,
            by = c('site', 'month', 'rep')) %>% 
  dplyr::mutate(abund_per_dm = total_order/dry_mass)

macro_density_sum <- macro_density %>% 
  group_by(site, month) %>% 
  summarise(mean_den = mean(abund_per_dm, na.rm = TRUE),
            se_den = sd(abund_per_dm, na.rm = TRUE)/length(mean_den))

ggplot(macro_density_sum,
       aes(x = month,
           y = mean_den,
           color = site,
           group = site))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin = mean_den - se_den,
                    ymax = mean_den + se_den),
                width = 0)+
  scale_x_continuous(breaks = c(0, 1, 3, 6, 12))+
  labs(y = expression(paste('Mean macroinvertebrate density (# g   ',DM^-1,')')),
       x = 'Month')+
  scale_color_viridis_d(name = 'Site')+
  theme(axis.title = element_text(size = 12))



macro_lm <- lm(data = macro_density,
               log10(abund_per_dm) ~ factor(month) * site)

summary(macro_lm)
anova(macro_lm)
car::Anova(macro_lm, type = 'III')
hist(macro_lm$residuals)
(agricolae::HSD.test(macro_lm, 'site'))
