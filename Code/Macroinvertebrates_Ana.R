####NMDS####

install.packages("vegan")
install.packages("permute")
install.packages("lattice")

library(permute)
library(lattice)
library(vegan)


data = read.csv("Data/from Ana/NMDS.csv", header= TRUE)

#make community matrix - extract columns with abundance information
com = data[,3:ncol(data)]

#turn abundance data frame into a matrix
m_com = as.matrix(com)

set.seed(123)
nmds = metaMDS(m_com, K=3, distance = "bray")
nmds 

plot(nmds)

#you need to obtain the coordinates for your NMDS1 and NMDS2 axes and put them in a new data frame:
data.scores = as.data.frame(scores(nmds)$sites)

#The newest version of the vegan package (>2.6-2) changed the format of the scores(nmds) object to a list, and so the above code will throw an error saying arguments imply differing number of rows. If you get this error try the following code instead to extract your site scores:

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))

#Next, you can add columns from your original data (data) to your new NMDS coordinates data frame. This will come in handy when you plot your data and want to differentiate groups or treatments:

#add columns to data frame 
data.scores$Stream = data$Stream

head(data.scores)

#Now we can plot our NMDS in ggplot2

library(ggplot2)

NMDS_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Stream, colour = Stream))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Stream", y = "NMDS2", shape = "Stream")  + 
  scale_colour_manual(values = c("deepskyblue4","azure4", "darkslategray4", "gold3","coral" )) 

NMDS_plot


####PERMANOVA####

#PERMANOVA
#Permutational multivariate analysis of variance (PERMANOVA) tests are used to determine whether the clustering patterns in NMDS plots represent consistent and significant differences in community composition)

install.packages(c("permute","lattice", "vegan"))

library(permute)
library(lattice)
library(vegan)


#Species matrix

Spp<-read.csv("Data/from Ana/PERMANOVA_Family.csv", header = TRUE)

#Treatment matrix (independent variable)

Treatments<-read.csv("Data/from Ana/PERMANOVA_Treatments.csv", header = TRUE)

adonis2(Spp ~ Stream, data = Treatments, permutations = 9999, method="bray", by= NULL)

#ANOSIM
#Analysis of similarities (ANOSIM) sI used to evaluate a dissimilarity matrix 
#ANOSIM is used after performing the NMDS to see if there are differences between the groups.
#Together, the dimension reduction and visualization capacities of NMDS and the hypothesis testing offered by ANOSIM are complementary approaches in evaluating nonparametric multivariate data.
# ANOSIM tests whether distances between groups are greater than within groups PERMANOVA tests whether distances differ between groups.

#Species matrix
Spp<-read.csv("Data/from Ana/PERMANOVA_Spp_Example.csv", header = TRUE)

#Treatment matrix (independent variable)
Treatments<-read.csv("Data/from Ana/PERMANOVA_Treatments_Example.csv", header = TRUE)

dune.dist <- vegdist(Spp)
attach(Treatments)
dune.ano <- anosim(dune.dist, Habitat, permutations = 9999, distance = "bray")
summary(dune.ano)
plot(dune.ano, xlab="Habitat",
     ylab= "Disimilarity", las=3, cex.lab=0)
boxplot(dune.ano, xlab="Period",
        ylab= "Disimilarity", las=3, cex.lab=0.8)


library(ggplot2)
autoplot(dune.ano, notch = FALSE)



####SIMPER#####
#Similarity Percentages (SIMPER) identify variables (e.g. species) that are likely to be the major contributors to any difference between groups detected by methods such as PERMANOVA.

library(permute)
library(lattice)
library(vegan)

Data_Spp=read.csv("Data/from Ana/PERMANOVA_Family.csv", header= TRUE)#Species matrix
Data_Treatments=read.csv("Data/from Ana/PERMANOVA_Treatments.csv", header= TRUE)#Treatment matrix (independent variable)

(sim <- with(Data_Treatments, simper(Data_Spp, Stream)))
summary(sim)

