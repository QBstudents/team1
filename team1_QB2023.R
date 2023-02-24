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

####Data wrangle####
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

#filling NAs with 0 
zoopdata3 <- na.omit(zoopdata2)

#calculating shannon diversity and species number 
zoopdata3$shannon    <- diversity(zoopdata3[,6:9], index = "shannon")
zoopdata3$specnum    <- specnumber(zoopdata3[,6:9])


#plotting shannon diversity over time 
ggplot(zoopdata3, aes(x = days, y = shannon))+
  geom_point()

ggplot(zoopdata3, aes(x = days, y = shannon))+
  geom_point()

#constructing a resemblance matrix 
dist      <- as.matrix(vegdist(zoopdata3[,6:9], "bray"))
interval  <- diff(zoopdata3$days, lag = 1)

levelplot(dist)


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


#also want to subset into site by species matrix and environmental matrix
zoopspecies <- zoopdata3[, 6:9]
zoopenv <- zoopdata3[, c(10:25)]


#perform dbRDA
zoop.dbrda <- dbrda(dist ~ ., as.data.frame(zoopenv))
ordiplot(zoop.dbrda)

#model only the intercept
zoop.dbrda.mod0 <- dbrda(dist ~ 1, as.data.frame(zoopenv))
ordiplot(zoop.dbrda.mod0)

#model full model, with all explanatory variables
zoop.dbrda.mod1 <- dbrda(dist ~ ., as.data.frame(zoopenv))

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
plot(scores(zoop.dbrda, display = "wa"), xlim = c(-1.3, 1.1),
     ylim = c(-1.1, 2.7), xlab = paste("dbRDA 1 (", dbrda.explainvar1, "%)",
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
#row.names(vectors) <- rownames(vectors)

arrows(0, 0, vectors[,1], vectors[,2],
       lwd = 2, lty = 1, length = 0.2, col = "red")
text(vectors[,1], vectors[,2], pos = 3,
     labels = row.names(vectors))
axis(side = 3, lwd.ticks = 2, cex.axis = 1.2, las = 1, col = "red", lwd = 2.2,
     at = pretty(range(vectors[,1]))* 2, labels = pretty(range(vectors[, 1])))
axis(side = 4, lwd.ticks = 2, cex.axis = 1.2, las = 1, col = "red", lwd = 2.2,
     at = pretty(range(vectors[, 2])) * 2, labels = pretty(range(vectors[, 2])))


#NMDS and scores

nmds1 <-metaMDS(dist, distance = "bray")
scores1 <- as.data.frame(scores(nmds1))

#subsetting nmds so can call as axis in plot
zoopdata3$axis1 <- scores1$NMDS1 
zoopdata3$axis2 <- scores1$NMDS2

#plotting NMDS
ggplot(zoopdata3, aes(x = days, y = axis1))+
  geom_point()+
  geom_smooth(method = "loess", span = 1)

ggplot(zoopdata3, aes(x = axis1, y = axis2, color = ))+
  geom_point()


#transposing data into long form (new row of species in own columns)
zoopdataLong <- zoopdata3 %>%
  pivot_longer(cols = dent:cyc, names_to = "species", values_to = "density") %>%
  filter(!Lake_Name == "Hale")

