#QB 2023 Team 1 Project#
#Lauren and Jonathan#


#clear environment of other files and set working directory to Team 1 folder in GitHub#
rm(list = ls())
getwd()

#install and load packages#
install.packages("plyr")
install.packages("vegan")
install.packages("readr")
library(plyr)
library(vegan)
library(readr)

####Data wrangle####
#import data set as zoopdata, look at structure
zoopdata <- read_csv("AllYears_Summary_111020.csv")
str(zoopdata)

#removing rounds with 'Marta'
zoopdata <- subset(zoopdata, Round != "Marta")

#make site by species matrix
#first changing lake names to numbered sites
zoopdata$Lake_Name <- revalue(zoopdata$Lake_Name, 
                  c("Airline" ="1", "Beaver Dam" ="2", "Beaver dam" ="2", "Benefiel" ="3", "Canvasback" = "4", "Dogwood" = "5", 
                    "Downing" = "6", "Gambill" = "7", "Goodman" = "8", "Goose" = "9", "Hale" = "10", "Island" = "11", "Long" = "12",
                    "Mayfield" = "13", "Midland" = "14", "Pump" = "15", "Scott" = "16", "T-Lake" = "17","T-lake"="17", "Todd" = "18", 
                    "University" ="19", "university" ="19", "Willow" = "20", "Wampler" ="21", "Walnut" = "22", "Trout" = "23", 
                    "Tree" = "24", "Sycamore" = "25", "Star" = "26", "Spencer" = "27", "South Lake" = "28", "Shouel West" = "29", 
                    "Shouel East" = "30", "Shop" = "31", "Shake1" = "32", "Shake 1" = "32", "Shake 2" ="33", "Redbud" = "34", "Narrow" = "35", 
                    "Lonnie" = "36", "Long-Hillenbrand" = "37", "Horseshoe" = "38", "Hackberry"="39", "Giluore" ="40", "Front"="41","Frank"="42",
                    "Crystal"="43","Corky"="44","Clear"="45","Chapel"="46"))
str(zoopdata)

#subset into site by species matrix
zoopspecies <- zoopdata[, 6:18]
#subset into environmental matrix
zoopenv <- zoopdata[, c(1:3, 5, 19:35)]










