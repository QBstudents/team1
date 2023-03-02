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


#removing rounds with 'Marta' 
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
#calculate some diversity indices
#function for eveness
Evar <- function(x){
  x <- as.vector(x[x > 0])
  1 - (2/pi) * atan(var(log(x)))
}

#calculating some alpha diversity indices
zoopdata3$shannon     <- diversity(zoopdata3[,6:9], index = "shannon")
zoopdata3$richness    <- specnumber(zoopdata3[,6:9])
zoopdata3$eveness     <- Evar(zoopdata3[,6:9])

#plotting over time for each lake
ggplot(zoopdata3, aes(x = year, y = shannon, group = Lake_Name, color = Lake_Name))+
  geom_jitter(aes(x = year, y = shannon), 
              size=3, width = 0.1)+
  labs (x = "year", y = "shannon diversity")+
  theme_bw()+
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 24, color = "black"),
        legend.title = element_blank())

ggplot(zoopdata3, aes(x = year, y = eveness, group = Lake_Name, color = Lake_Name))+
  geom_jitter(aes(x = year, y = eveness), 
              size=3, width = 0.1)+
  labs (x = "year", y = "eveness")+
  theme_bw()+
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 24, color = "black"),
        legend.title = element_blank())

####VISUALIZING ENV CHARACTERISTICS####
#environmental characteristics#
#calculating N:P ratio
zoopenv$NP_ratio <- zoopenv$N_ug/zoopenv$P_ug

#see which characteristics have correlation
ggcorr(zoopenv[c(3:20)],method=c("pairwise","pearson"))

#creating a long version of the environmental matrix
zoopenvlong <- zoopenv %>%
  pivot_longer(cols = Edible_Chl:NP_ratio, names_to = "env", values_to = "values")

zoopenvlong$year <- format(as.Date(zoopenvlong$date, format="%Y/%m/%d"),"%Y")
  
#Plot through time
ggplot(zoopenvlong, aes(x = year, y = log(values), group = Lake_Name, color = Lake_Name))+
  geom_jitter(aes(x = year, y = log(values)), 
              size=3, width = 0.1)+
  labs (x = "year", y = "metric")+
  facet_wrap(~env)


####BETA DIVERSITY####
#constructing a resemblance matrix 
dist      <- as.matrix(vegdist(zoopdata3[,6:9], "bray"))
interval  <- diff(zoopdata3$days, lag = 1)

#zoop pcoa to visualize beta diversity 
zoop.pcoa <- cmdscale(dist, eig = TRUE, k = 3)

explainvar1 <- round(zoop.pcoa$eig[1] / sum(zoop.pcoa$eig), 3) * 100
explainvar2 <- round(zoop.pcoa$eig[2] / sum(zoop.pcoa$eig), 3) * 100
explainvar3 <- round(zoop.pcoa$eig[3] / sum(zoop.pcoa$eig), 3) * 100 
sum.eig     <- sum(explainvar1, explainvar2, explainvar3)

par(mar = c(5,5,1,2) + 0.1)

plot(zoop.pcoa$point[ ,1], zoop.pcoa$points[ ,2], ylim = c(-0.2, 0.7),
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5,
     cex.axis = 1.2, axes = FALSE
)

axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

points(zoop.pcoa$points[,1], zoop.pcoa$points[,2],
       pch = 19, cex = 3, bg = "gray", col = "gray")
text(zoop.pcoa$points[ ,1], zoop.pcoa$points[ ,2],
     labels = row.names(zoop.pcoa$points))

#perform dbRDA
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

ggplot(zoopenv, aes(x = date, y = N_ug))+
  geom_jitter(aes(x = date, y = N_ug), 
              size=3, width = 0.1)+
  labs (x = "year", y = "metric")

#creating vectors of environmental characteristics that were significant
median(zoopenv$N_ug)
min(zoopenv$N_ug)
max(zoopenv$N_ug)

median(zoopenv$mT)
min(zoopenv$mT)
max(zoopenv$mT)

median(zoopenv$mspc)
min(zoopenv$mspc)
max(zoopenv$mspc)

zoopenv.2 <- zoopenv%>%
  mutate(Nug   = case_when(N_ug < 364 ~ 'low',
                           N_ug < 500 ~ 'med',
                           N_ug > 501 ~ 'high'))%>%
  mutate(meanT = case_when(mT < 10 ~ 'low',
                           mT < 19.1186 ~ 'med',
                           mT > 19.1186 ~ 'high'))%>%
  mutate(mspc  = case_when(mspc  < 1000.00~ 'low',
                           mspc  < 1750.00 ~ 'med',
                           mspc  > 1751.00 ~ 'high'))

Nug <- zoopenv.2[,21]
table(Nug)
Nug <- c(rep("low", 70), rep("med", 41), rep("high", 31))

meanT <- zoopenv.2[,21]
table(meanT)
meanT <- c(rep("high", 71), rep("med", 65), rep("low", 6))

mspc <- zoopenv.2[,21]
table(mspc)
mspc <- c(rep("high", 31), rep("med", 41), rep("low", 70))

#species associations#
zoop.rel <- decostand(zoopspecies, method = "total")
phi2 <- multipatt(zoop.rel, cluster = Nug, func = "r.g", control = how(nperm = 999))
summary(phi2)

phi3 <- multipatt(zoop.rel, cluster = meanT, func = "r.g", control = how(nperm = 999))
summary(phi3)
#nothing significant with temperature

phi4 <-  multipatt(zoop.rel, cluster = mspc, func = "r.g", control = how(nperm = 999))
summary(phi4)

#plotting density of each species with N_ug
ggplot(zoopdataLong, aes(x = N_ug, y = log(density)))+
  geom_point()+
  geom_smooth(method = NULL, colour = "#CD9B1D")+
  labs(x = "Nitrogen (ug)", y = "log(density)")+
  theme_bw()+
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 24, color = "black"),
        legend.title = element_blank())+
  facet_wrap(~species)

ggplot(zoopdataLong, aes(x = mspc, y = log(density)))+
  geom_point()+
  geom_smooth(method = NULL, colour = "#CD9B1D")+
  labs(x = "Conductivity", y = "log(density)")+
  theme_bw()+
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 24, color = "black"),
        legend.title = element_blank())+
  facet_wrap(~species)

#NMDS and scores
nmds1 <-metaMDS(dist, distance = "bray")
scores1 <- as.data.frame(scores(nmds1))

#subsetting nmds so can call as axis in plot
zoopdata3$axis1 <- scores1$NMDS1 
zoopdata3$axis2 <- scores1$NMDS2

#plotting NMDS
ggplot(zoopdata3, aes(x = days, y = axis1))+
  geom_point()+
  geom_smooth(method = "loess", span = 1, colour = "#CD9B1D")+
  theme_bw()+
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 24, color = "black"),
        legend.title = element_blank())
