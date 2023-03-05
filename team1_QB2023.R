#QB 2023 Team 1 Project#
#Lauren and Jonathan#

#clear environment of other files and set working directory to Team 1 folder in GitHub#
rm(list = ls())
getwd()

#install and load packages#
library(plyr)
library(vegan)
library(readr)
library(tidyverse)
library(codyn)
library(ade4)
library(indicspecies)
library(GGally)

####DATA WRANGLE####
#import data set as zoopdata, look at structure
zoopdata <- read_csv("AllYears_Summary_111020.csv")
str(zoopdata)


#removing rounds with 'Marta' and replacing zeros with NA
zoopdata2 <- subset(zoopdata, Round != "Marta")
zoopdata2[zoopdata2 == 0] <- NA   
str(zoopdata2)
#formatting date, creating a days column that counts days since first date entry (2009-08-09)
startdate <- as.Date("2009-08-09","%Y-%m-%d")
zoopdata2$days <- as.numeric(difftime(zoopdata2$date, startdate, unit = "days"))

#get rid of NAs
zoopdata3 <- na.omit(zoopdata2)
#creating year column
zoopdata3$year <- format(as.Date(zoopdata3$date, format="%Y/%m/%d"),"%Y")


#also want to subset into site by species matrix and environmental matrix
zoopspecies <- zoopdata3[, 6:9]
zoopenv     <- zoopdata3[, c(1,3,5,10:25)]


####ALPHA DIVERSITY####
#function for simpson diversity
SimpD <- function(x =""){
  D = 0
  N = sum(x)
  for (n_i in x){
    D = D + (n_i^2)/(N^2)
  }
  return(D)
}


#calculating some alpha diversity indices
zoopdata3$simpson     <- SimpD(zoopdata3[,6:9])
zoopdata3$shannon     <- diversity(zoopdata3[,6:9], index = "shannon")

#Hill Numbers 
zoopdata3$expshannon  <- exp(diversity(zoopdata3[,6:9], index = "shannon")) #q = 1
zoopdata3$richness    <- specnumber(zoopdata3[,6:9]) #q = 0
zoopdata3$invsimp     <- 1/SimpD(zoopdata3[,6:9]) #q = 2

#plotting over time for each lake
ggplot(zoopdata3, aes(x = year, y = shannon, group = Lake_Name, color = Lake_Name))+
  geom_jitter(aes(x = year, y = shannon), 
              size=3, width = 0.1, alpha = .5)+
  geom_line(aes(x = year, y = shannon, color = Lake_Name, group = Lake_Name))+
  labs (x = "year", y = "shannon diversity")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        legend.title = element_blank())+
  facet_wrap(~Lake_Name)+
  theme(legend.position = "none")

ggplot(zoopdata3, aes(x = year, y = expshannon, group = Lake_Name, color = Lake_Name))+
  geom_jitter(aes(x = year, y = expshannon), 
              size=3, width = 0.1, alpha = .5)+
  geom_line(aes(x = year, y = expshannon, color = Lake_Name, group = Lake_Name))+
  labs (x = "year", y = "exponential shannon entropy")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        legend.title = element_blank())+
  facet_wrap(~Lake_Name)+
  theme(legend.position = "none")
  
ggplot(zoopdata3, aes(x = year, y = simpson, group = Lake_Name, color = Lake_Name))+
  geom_jitter(aes(x = year, y = simpson), 
              size=3, width = 0.1, alpha = .5)+
  geom_line(aes(x = year, y = simpson, color = Lake_Name, group = Lake_Name))+
  labs (x = "year", y = "simpson diversity")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        legend.title = element_blank())+
  facet_wrap(~Lake_Name)+
  theme(legend.position = "none")

ggplot(zoopdata3, aes(x = year, y = richness, group = Lake_Name, color = Lake_Name))+
  geom_jitter(aes(x = year, y = richness), 
              size=3, width = 0.1, alpha = .5)+
  labs (x = "year", y = "species richness")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        legend.title = element_blank())+
  facet_wrap(~Lake_Name)+
  theme(legend.position = "none")
  
ggplot(zoopdata3, aes(x = year, y = invsimp, group = Lake_Name, color = Lake_Name))+
  geom_jitter(aes(x = year, y = invsimp), 
              size=3, width = 0.1, alpha = .5)+
  geom_line(aes(x = year, y = invsimp, color = Lake_Name, group = Lake_Name))+
  labs (x = "year", y = "inverse simpson diversity")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        legend.title = element_blank())+
  facet_wrap(~Lake_Name)+
  theme(legend.position = "none")

ggplot(zoopdata3, aes(x = richness, y = shannon, group = Lake_Name, color = Lake_Name))+
  geom_jitter(aes(x = richness, y = shannon), 
              size=3, width = 0.1, alpha = .5)+
  labs (x = "species richness", y = "shannon diversity")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        legend.title = element_blank())+
  facet_wrap(~Lake_Name)+
  theme(legend.position = "none")


ggplot(zoopdata3, aes(x = richness, y = expshannon, group = Lake_Name, color = Lake_Name))+
  geom_jitter(aes(x = richness, y = expshannon), 
              size=3, width = 0.1, alpha = .5)+
  labs (x = "species richness", y = "exponential shannon diversity")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        legend.title = element_blank())+
  facet_wrap(~Lake_Name)+
  theme(legend.position = "none")

####VISUALIZING ENVIRONMENTAL PARAMETERS####
#calculating N:P ratio
zoopenv$NP_ratio <- zoopenv$N_ug/zoopenv$P_ug

#see if any have correlation
ggcorr(zoopenv[c(3:20)],method=c("pairwise","pearson"))

#creating a long version of the environmental matrix
zoopenvlong <- zoopenv %>%
  pivot_longer(cols = Edible_Chl:NP_ratio, names_to = "env", values_to = "values")

#make long form of environmental matrix
zoopenvlong$year <- format(as.Date(zoopenvlong$date, format="%Y/%m/%d"),"%Y")
  
#Plot through time
#scale function to get z scores of environmental parameters   
zoopenvlong$scalevalues <- scale(zoopenvlong$values)
#making plot over time for each lake
ggplot(zoopenvlong, aes(x = year, y = scalevalues, group = env, color = env))+
  ylim(c(-5.0,10))+
  geom_jitter(aes(x = year, y = scalevalues), 
              position = "jitter")+
  labs (x = "year", y = "scaled environmental values")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        legend.title = element_blank())+
  facet_wrap(~Lake_Name)

####BETA DIVERSITY####
#constructing a resemblance matrix using Bray-Curtis 
dist      <- as.matrix(vegdist(zoopdata3[,6:9], "bray"))
interval  <- diff(zoopdata3$days, lag = 1)

#dbRDA
zoop.dbrda <- dbrda(dist ~ ., as.data.frame(zoopenv[,3:20]))
ordiplot(zoop.dbrda)

#model only the intercept
zoop.dbrda.mod0 <- dbrda(dist ~ 1, as.data.frame(zoopenv[,3:20]))
ordiplot(zoop.dbrda.mod0)

#model full model, with all explanatory variables
zoop.dbrda.mod1 <- dbrda(dist ~ ., as.data.frame(zoopenv[,3:20]))

#all combinations of explanatory variables in model-- function returns model with lowest AIC value
zoop.dbrda <- ordiR2step(zoop.dbrda.mod0, zoop.dbrda.mod1, perm.max = 200)

zoop.dbrda$call
zoop.dbrda$anova
ordiplot(zoop.dbrda)

#permutation tests to evaluate significance
permutest(zoop.dbrda, permutations = 999)
envfit(zoop.dbrda, zoopenv, perm = 999)

#calculate explained variation
dbrda.explainvar1 <- round(zoop.dbrda$CCA$eig[1] / sum(c(zoop.dbrda$CCA$eig, zoop.dbrda$CA$eig)), 3) * 100
dbrda.explainvar2 <- round(zoop.dbrda$CCA$eig[2] / sum(c(zoop.dbrda$CCA$eig, zoop.dbrda$CA$eig)), 3) * 100

#define plot parameters
par(mar = c(5, 5, 4, 4) + 0.1)

#initiate plot
plot(scores(zoop.dbrda, display = "wa"), xlim = c(-1.8, 2.5),
     ylim = c(-3.5, 3.0), xlab = paste("dbRDA 1 (", dbrda.explainvar1, "%)",
                                       sep = ""), ylab = paste("dbRDA 2 (", dbrda.explainvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

#add axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

#add points and labels
points(scores(zoop.dbrda, display = "wa"),
       pch = 19, cex = 3, bg = "gray", col = "gray")
text(scores(zoop.dbrda, display = "wa"),
     labels = row.names(scores(zoop.dbrda, display = "wa")))

#add environmental vectors
vectors <- scores(zoop.dbrda, display = "bp")

arrows(0, 0, vectors[,1], vectors[,2],
       lwd = 2, lty = 1, length = 0.2, col = "red")
text(vectors[,1], vectors[,2], pos = 3,
     labels = row.names(vectors))
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, col = "red", lwd = 2.2,
     at = pretty(range(vectors[,1]))* 2, labels = pretty(range(vectors[, 1])))
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, col = "red", lwd = 2.2,
     at = pretty(range(vectors[, 2])) * 2, labels = pretty(range(vectors[, 2])))


#creating vectors of environmental parameters that were "significant" 
#first find median, min, max to make arbitrary limits to each vector
#nitrogen
median(zoopenv$N_ug)
min(zoopenv$N_ug)
max(zoopenv$N_ug)

#mean temperature
median(zoopenv$mT)
min(zoopenv$mT)
max(zoopenv$mT)

#conductivity
median(zoopenv$mspc)
min(zoopenv$mspc)
max(zoopenv$mspc)

#mean dissolved oxygen
median(zoopenv$mDO)
min(zoopenv$mDO)
max(zoopenv$mDO)

#new columns of vectors 
zoopenv.2 <- zoopenv%>%
  mutate(Nug   = case_when(N_ug < 364 ~ 'low',
                           N_ug < 500 ~ 'med',
                           N_ug > 501 ~ 'high'))%>%
  mutate(meanT = case_when(mT < 10 ~ 'low',
                           mT < 19.1186 ~ 'med',
                           mT > 19.1186 ~ 'high'))%>%
  mutate(cond  = case_when(mspc  < 1000.00~ 'low',
                           mspc  < 1750.00 ~ 'med',
                           mspc  > 1751.00 ~ 'high'))%>%
  mutate(DO = case_when(mDO  < 4.0 ~ 'low',
                         mDO < 7.90685 ~ 'med',
                         mDO >  7.90686 ~ 'high'))

Nug <- zoopenv.2[,21]
table(Nug)
Nug <- c(rep("low", 70), rep("med", 41), rep("high", 31))

meanT <- zoopenv.2[,22]
table(meanT)
meanT <- c(rep("high", 71), rep("med", 65), rep("low", 6))

cond <- zoopenv.2[,23]
table(cond)
cond <- c(rep("high", 31), rep("med", 87), rep("low", 24))

DO <- zoopenv.2[,24]
table(DO)
DO <- c(rep("high", 71), rep("med", 69), rep("low", 2))


#species preference to the selected environmental parameters#
zoop.rel <- decostand(zoopspecies, method = "total")
phi2 <- multipatt(zoop.rel, cluster = Nug, func = "r.g", control = how(nperm = 999))
summary(phi2)

phi3 <- multipatt(zoop.rel, cluster = meanT, func = "r.g", control = how(nperm = 999))
summary(phi3)
#nothing significant with temperature

phi4 <-  multipatt(zoop.rel, cluster = cond, func = "r.g", control = how(nperm = 999))
summary(phi4)

phi5 <-  multipatt(zoop.rel, cluster = DO, func = "r.g", control = how(nperm = 999))
summary(phi5)

#creating a long version of the species matrix
zooplong <- zoopdata3 %>%
  pivot_longer(cols = dent:Bosmina, names_to = "species", values_to = "density")

#plotting density of each species with N_ug
ggplot(zooplong, aes(x = N_ug, y = log(density)))+
  geom_point()+
  geom_smooth(method = NULL, colour = "#CD9B1D")+
  labs(x = "Nitrogen (ug)", y = "log(density)")+
  theme_bw()+
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 24, color = "black"),
        legend.title = element_blank())+
  facet_wrap(~species)
#plotting density of each species with conductivity 
ggplot(zooplong, aes(x = mspc, y = log(density)))+
  geom_point()+
  geom_smooth(method = NULL, colour = "#CD9B1D")+
  labs(x = "Conductivity", y = "log(density)")+
  theme_bw()+
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 24, color = "black"),
        legend.title = element_blank())+
  facet_wrap(~species)


#NMDS and scores
nmds2 <- metaMDS(zoopspecies, distance = "bray")
scores2 <- as.data.frame(scores(nmds2)$sites)

#environmental vectors
en <- envfit(nmds2, zoopenv, perm = 999)
env_vectors <- as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

scores2$Lake_Name = zoopdata3$Lake_Name

ggplot(data = scores2, aes(x = NMDS1, y = NMDS2))+
  geom_point(data = scores2, aes(colour = Lake_Name), size = 3, alpha = 0.5)+
  coord_fixed()+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = env_vectors, colour = "#CD9B1D",
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_text(data = env_vectors, aes(x = NMDS1, y = NMDS2), colour = "black",
            label = row.names(env_vectors), size = 2.0)+
  theme_bw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        legend.title = element_blank())+
  theme(legend.position = "bottom")
  