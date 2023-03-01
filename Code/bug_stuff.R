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

# raw spreadsheet from Ana had a site 'Saltito 30'.
# Mar 1 2023: Change to Saltito 60, assuming this was a typo

unique(data$Stream) # should have 5 sites

# put in order of increasing conductivity
data$Stream <- fct_relevel(data$Stream,
                           c('Taconazo 30', 'Piper', 'Saltito 60',
                             'Sura 30', 'Arboleda 30'))

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
  filter(QAQC == 0) 


# 3) NMDS's ----

# 3.1) All dates & sites, at Order level

# which orders are present
unique_orders <- na.omit(unique(tax_data_clean$order))

length(unique_orders)

tax_data_clean_order <- tax_data_clean %>% 
  filter(!is.na(order)#,
         #order != 'Neoophora'
         ) %>%                               # remove NAs from Order
  select(Stream, order, month, rep, Total) %>%    # get the necessary columns
  group_by(Stream, month, order) %>%         # and do the grouping
  summarise(mean_order = mean(Total, na.rm = TRUE)) %>%   # sum by order in each possible grouping 
  mutate(log_total = log10(1 + mean_order)) %>%          # log10 + 1 transform data
  select(-mean_order) %>% 
  pivot_wider(names_from = order,                         # pivot data
              values_from = log_total,
              values_fill = 0) 


tax_data_matrix <- tax_data_clean_order[,-c(1:2)]
tax_data_meta <- tax_data_clean_order[,c(1:2)]

nmds_order <- metaMDS(tax_data_matrix,
                      distance = 'bray',
                      k = 3,
                      autotransform = FALSE)
stress <- round(nmds_order$stress, 3)

nmds_order_out <- data.frame(x = nmds_order$points[,1],
                             y = nmds_order$points[,2])

nmds_order_out <- cbind(tax_data_meta,
                        nmds_order_out)

fit <- (envfit(nmds_order,
               tax_data_matrix,
               perm = 9999))

scrs <- data.frame(scores(fit, 'vectors'))

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
  geom_label_repel(data = scrs_sig,
                   aes(NMDS1, NMDS2,
                       label = env.variables),
                   size = 3)+
  labs(x = 'NMDS 1',
       y = 'NMDS 2')+
  lims(x = c(-1,1),
       y = c(-1,1))+
  geom_label(label = paste('3D Stress = ',stress),
           x = -0.55,
           y = 0.9)
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
  labs(x = "NMDS axis 1", 
       y = "NMDS axis 2")+
  lims(x = c(-1,1),
       y = c(-1,1))+
  theme(legend.background = element_blank())+
  scale_color_viridis_d(name = element_blank())
plot_nmds_streams

ggpubr::ggarrange(plot_nmds_order_scrs,
                  plot_nmds_streams,
                  ncol = 2,
                  align = 'h', widths = c(1,1.4))
