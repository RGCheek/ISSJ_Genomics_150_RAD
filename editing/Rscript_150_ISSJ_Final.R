#############################################
# ISSJ RDA Analysis Final version       #####
# Author: Rebecca Cheek                 #####
# Nov 2019                             #####
#############################################

library(data.table)
library(adegenet)
library(ape)
library("pegas")
library("seqinr")
library("ggplot2")
library("hierfstat")
library(LEA)
library(StAMPP)
library(dplyr)
library(RColorBrewer)
library(mapplots)
library(fields)
library(stringr)
library(pcadapt)
library(qvalue)
library(ggthemes)
library(vegan)
library(qqman)
library(dartR)
library(sf)
library(lme4)
library(psych)
library(sp)
library(dartR)
library(parallel)




###############################################################
# RDA with more environmental predictors 

#############################################################

##Set working directory to

setwd("C:/Users/Rebecca/Dropbox/ISSJ/issj_genomics/r_issj_genomics")
#

#Convert genomic data. Using the .raw file from Plink 

gen <- read.PLINK("C:/Users/Rebecca/Dropbox/ISSJ/issj_genomics/r_issj_genomics/ISSJ_clean.raw", parallel = F)
dim(gen)

#Both LFMM and RDA require complete data frames (i.e., no missing genetic data). 
#For this example, we'll use a simple approach to imputing missing genotype values: 
#we will impute using the most common genotype at each SNP across all individuals.

gen <- as.matrix(gen)

sum(is.na(gen)) 

gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), 
                                             as.numeric(names(which.max(table(x))))))

sum(is.na(gen.imp)) # No NAs



genfile <- as.matrix(gen.imp [,-1], row.names = FALSE, col.names=FALSE)

mode(genfile) <- 'numeric'  #needs to be numeric for some reason and matrix was character


#Read in environmental data
env <- read.csv("../ISSJ_env2.csv")

pred <- env %>% 
  dplyr::select(-c("Latitude", "Longitude")) %>% 
  filter(tissue_number != "F11-093") %>% 
  filter(tissue_number != "F11-156") %>%
  filter(tissue_number != "F09-11") %>%
  filter(tissue_number != "F09-45") %>%
  filter(tissue_number != "F10-44") %>%
  filter(tissue_number != "F10-154") %>%
  filter(tissue_number != "F09-27") %>% 
  filter(tissue_number != "F10-162") %>% 
  filter(tissue_number != "F09-141") %>%
  filter(tissue_number != "F11-170") 


pred$tissue_number <- as.character(pred$tissue_number) # Make individual names characters (not factors)


# Confirm that genotypes and environmental data are in the same order
identical(rownames(gen.imp), pred[,1])

#filter tissue names out of env file since we're not testing those
pred <- pred%>% 
  dplyr::select(-c("tissue_number","POP")) 

#remove some of the predictors that are correlated (everything is fekin correlated)
pred <- pred %>%
  dplyr::select(-c("NFFD","difference_coldest_and_warmest_month_measure_continentality_deg_C", "NFFD",
                   "mean_summer_precip_mm", "degree.days_above_5_deg_C", "degree.days_above_18_deg_C","extreme_min_temp_over_30_years",
                   "annual_heat_moisture_index","Hargreave_climatic_moisture_index","Hargreave_reference_evaporation",
                   "summer_heat_moisture_index","distance","degree.days_below_18_deg_C", "mean_temp__warmest_month_deg_C", "mean_temp_coldest_month_deg_C"))

#check correlations of predictors
pairs.panels(pred[,1:5], scale=T)


#temp and precip are correlated, but since we're dealing with such a small scale I'm not freaking surprised. 

#First, will test with temp and percent pine and percent oak
pred <- pred %>%
  dplyr::select(-c("mean_annual_precip_mm", "percent_other"))

#check correlations of predictors
pairs.panels(pred[,1:3], scale=T)

#Run the RDA
issj.rda <- vegan::rda(genfile~., data=pred, scale=T)
issj.rda

vegan::RsquareAdj(issj.rda) #r2= 0.027, adusted r2= 0.005

summary(eigenvals(issj.rda, model="constrained"))

#We can visualize this information using a screeplot of the canonical eigenvalues by calling `screeplot`:

screeplot(issj.rda)


signif.full <- anova.cca(issj.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

signif.axis <- anova.cca(issj.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis


#You can run a formal test of statistical significance of each constrained axis using:
#anova.cca(issj.rda, by="axis") #We can assess both the full model and each constrained axis using F-statistics
# (Legendre et al, 2010). The null hypothesis is that no linear relationship exists between the SNP data and the
#'# environmental predictors. See `?anova.cca` for more details and options.


#The permutation process to test the signficiance of each axis takes a while
#(up to a few hours on large data sets), so we'll just use the screeplot for a first assessment.
#If we did run the formal test, we would find that the first three constrained axes are significant (p = 0.001)

#Finally, `vegan` has a simple function for checking Variance Inflation Factors for the predictor variables used in the model:

vegan::vif.cca(issj.rda) 

#All values are below 10, co collinarity shouldn't be a problem. 

#quick plot of the RDA output using the default plotting in `vegan`

plot(issj.rda, scaling=3)          # default is axes 1 and 2


########################################
#Plotting with fancy figure         ###
########################################

###### To highlight individuals by population
env <- read.csv("../Pop_assignment3.csv")

env <- env %>%
  filter(tissue_number != "F11-093") %>% 
  filter(tissue_number != "F11-156") %>%
  filter(tissue_number != "F09-11") %>%
  filter(tissue_number != "F09-45") %>%
  filter(tissue_number != "F10-44") %>%
  filter(tissue_number != "F10-154") %>%
  filter(tissue_number != "F09-27") %>% 
  filter(tissue_number != "F10-162") %>% 
  filter(tissue_number != "F09-141") %>%
  filter(tissue_number != "F11-170") 


eco <- env %>% 
  mutate(POP = recode(POP,"E-O"="Eastern Oak","W-O"="Western Oak","C-O"="Central Oak",
                      "E-P"="Eastern Pines","W-P"="Western Pines","C-P"="Central Pines"))

myCol2 <- c("darkorange","deepskyblue1", "yellow", "lightcyan2","orangered2","navy")

#Plot the individuals on axis 1 &2 

plot(issj.rda, type="n", scaling=3)
points(issj.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(issj.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=myCol2[eco$POP]) # the jays
text(issj.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco$POP), bty="n", col="gray32", pch=21, cex=1, pt.bg=myCol2)

# axes 1 & 3
plot(issj.rda, type="n", scaling=3, choices=c(1,3))
points(issj.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(issj.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=myCol2[eco$POP], choices=c(1,3))
text(issj.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topleft", legend=levels(eco$POP), bty="n", col="gray32", pch=21, cex=1, pt.bg=myCol2)

#########################################
### b) Identify RDA candidates        ###
#########################################

#Use the loadings of the SNPs (their location) in the ordination
#space to determine which SNPs are candidates for local adaptation.
#The SNP loadings are stored as `species` in the RDA object. We'll extract the SNP loadings from the first three constrained axes:

load.rda <- scores(issj.rda, choices=c(1:3), display="species")

#If we look at histograms of the loadings on each RDA axis,
#we can see their (relatively normal) distribution. SNPs loading at the center of the distribution are not showing
#a relationship with the environmental predictors; those loading in the tails are, and are more likely to be under 
#selection as a function of those predictors (or some other predictor correlated with them).


hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 


#define the function here as `outliers`, where `x` is the vector of loadings and `z` is the number of standard deviations to use:

#find loadings +/-z sd from mean loading
#locus names of the tails

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)       
  x[x < lims[1] | x > lims[2]]           
}

#Now let's apply it to the first three constrained axes using sd of 3.5 (strong selection)

# cand1 <- outliers(load.rda[,1],3.5) # 11
# cand2 <- outliers(load.rda[,2],3.5) # 1
# cand3 <- outliers(load.rda[,3],3.5) # 9
# 
# ncand <- length(cand1)+length(cand2)+length(cand3)
# ncand  ##22 candidates


#using an sd of 3 

ecand1 <- outliers(load.rda[,1],3) # 35
ecand2 <- outliers(load.rda[,2],3) # 54
ecand3 <- outliers(load.rda[,3],3) # 37



encand <- length(ecand1)+length(ecand2)+length(ecand3)
encand

#We have 102 potential candidate SNPs under moderate selection

#make a single data frame with the axis, SNP name, & loading:

ecand1 <- cbind.data.frame(rep(1,times=length(ecand1)), names(ecand1), unname(ecand1))
ecand2 <- cbind.data.frame(rep(2,times=length(ecand2)), names(ecand2), unname(ecand2))
ecand3 <- cbind.data.frame(rep(3,times=length(ecand3)), names(ecand3), unname(ecand3))

colnames(ecand1) <- colnames(ecand2)<- colnames(ecand3) <- c("axis","snp","loading")


ecand <- rbind(ecand1, ecand2, ecand3)

ecand$snp <- as.character(ecand$snp)

foo <- matrix(nrow=(encand), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("mean_annual_temp_deg_C", "percent_pine", "percent_oak")

for (i in 1:length(ecand$snp)) {
  nam <- ecand[i,2]
  snp.gen <- genfile[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

ecand <- cbind.data.frame(ecand,foo)  
head(ecand)

dim(ecand)

#Some of these may be duplicate detections; let's check:


# duplicate detections:
length(ecand$snp[duplicated(ecand$snp)])
foo <- cbind(ecand$axis,duplicated(ecand$snp))
table(foo[foo[,1]==1,2])
table(foo[foo[,1]==2,2])
table(foo[foo[,1]==3,2])

ecand <- ecand[!duplicated(ecand$snp),] #remove duplicates

dim(ecand)

## 102 unique outlier loci of 5768 SNPs

#which predictor each snp is more correlated with

for (i in 1:length(ecand$snp)) {
  bar <- ecand[i,]
  ecand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable
  ecand[i,8] <- max(abs(bar[4:6]))              # gives the correlation
}

colnames(ecand)[7] <- "predictor"
colnames(ecand)[8] <- "correlation"

table(ecand$predictor)
#mean_annual_temp_deg_C            percent_oak           percent_pine 
#30                                  33                     39                

##To highlight SNPs
sel <- ecand$snp
env <-ecand$predictor

env[env=="mean_annual_temp_deg_C"] <- '#e31a1c'
env[env=="percent_pine"] <- '#81be41'
env[env=="percent_oak"] <- '#ffde00'


# color by predictor:
col.pred <- rownames(issj.rda$CCA$v) # pull the SNP names


for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

colors <- c("#e31a1c","#81be41","#ffde00")

col.pred[!grepl(paste(colors,collapse="|"), col.pred)] <- '#f1eef6' # non-candidate SNPs


empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent

empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c("#e31a1c","#81be41","#ffde00")


# axes 1 & 2
plot(issj.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(issj.rda, display="species", pch=21, cex=1, col="grey32", bg=col.pred, scaling=3)
points(issj.rda, display="species", pch=21, cex=1,col=empty.outline, bg=empty, scaling=3)
text(issj.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomleft",legend=c("Mean Temp","Percent Pine","Percent Oak"),bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

#axes 1 &3 
plot(issj.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(issj.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(issj.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(issj.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomleft",legend=c("Mean Temp","Percent Pine","Percent Oak"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

####################################################
# Running RDA with just pine and oak
###################################################

#Convert genomic data. Using the .raw file from Plink 

gen <- read.PLINK("C:/Users/Rebecca/Dropbox/ISSJ/issj_genomics/r_issj_genomics/ISSJ_clean.raw", parallel = F)
dim(gen)

#Both LFMM and RDA require complete data frames (i.e., no missing genetic data). 
#For this example, we'll use a simple approach to imputing missing genotype values: 
#we will impute using the most common genotype at each SNP across all individuals.

gen <- as.matrix(gen)

sum(is.na(gen)) 

gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), 
                                             as.numeric(names(which.max(table(x))))))

sum(is.na(gen.imp)) # No NAs



genfile <- as.matrix(gen.imp [,-1], row.names = FALSE, col.names=FALSE)

mode(genfile) <- 'numeric'  #needs to be numeric for some reason and matrix was character


#Read in environmental data
env <- read.csv("../ISSJ_env2.csv")

pred <- env %>% 
  dplyr::select(-c("Latitude", "Longitude")) %>% 
  filter(tissue_number != "F11-093") %>% 
  filter(tissue_number != "F11-156") %>%
  filter(tissue_number != "F09-11") %>%
  filter(tissue_number != "F09-45") %>%
  filter(tissue_number != "F10-44") %>%
  filter(tissue_number != "F10-154") %>%
  filter(tissue_number != "F09-27") %>% 
  filter(tissue_number != "F10-162") %>% 
  filter(tissue_number != "F09-141") %>%
  filter(tissue_number != "F11-170") 


pred$tissue_number <- as.character(pred$tissue_number) # Make individual names characters (not factors)


# Confirm that genotypes and environmental data are in the same order
identical(rownames(gen.imp), pred[,1])

#filter tissue names out of env file since we're not testing those
pred <- pred%>% 
  dplyr::select(-c("tissue_number","POP")) 

#remove some of the predictors that are correlated (everything is fekin correlated)
pred <- pred %>%
  dplyr::select(-c("NFFD","difference_coldest_and_warmest_month_measure_continentality_deg_C", "NFFD",
                   "mean_summer_precip_mm", "degree.days_above_5_deg_C", "degree.days_above_18_deg_C","extreme_min_temp_over_30_years",
                   "annual_heat_moisture_index","Hargreave_climatic_moisture_index","Hargreave_reference_evaporation",
                   "summer_heat_moisture_index","distance","degree.days_below_18_deg_C", "mean_temp__warmest_month_deg_C", "mean_temp_coldest_month_deg_C"))

#check correlations of predictors
pairs.panels(pred[,1:5], scale=T)


#temp and precip are correlated, but since we're dealing with such a small scale I'm not freaking surprised. 

#First, will test with percent pine and percent oak
pred <- pred %>%
  dplyr::select(-c("mean_annual_precip_mm", "percent_other", "mean_annual_temp_deg_C"))

#Run the RDA
issj.rda <- vegan::rda(genfile~., data=pred, scale=T)
issj.rda

vegan::RsquareAdj(issj.rda) #r2= 0.0180, adusted r2= 0.0037

summary(eigenvals(issj.rda, model="constrained"))

#We can visualize this information using a screeplot of the canonical eigenvalues by calling `screeplot`:

screeplot(issj.rda)


signif.full <- anova.cca(issj.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

signif.axis <- anova.cca(issj.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis


#You can run a formal test of statistical significance of each constrained axis using:
anova.cca(issj.rda, by="axis") #We can assess both the full model and each constrained axis using F-statistics
# (Legendre et al, 2010). The null hypothesis is that no linear relationship exists between the SNP data and the
#'# environmental predictors. See `?anova.cca` for more details and options.


#The permutation process to test the signficiance of each axis takes a while
#(up to a few hours on large data sets), so we'll just use the screeplot for a first assessment.
#If we did run the formal test, we would find that the first three constrained axes are significant (p = 0.001)

#Finally, `vegan` has a simple function for checking Variance Inflation Factors for the predictor variables used in the model:

vegan::vif.cca(issj.rda) 

#All values are below 10, co collinarity shouldn't be a problem. 

#quick plot of the RDA output using the default plotting in `vegan`

plot(issj.rda, scaling=3)          # default is axes 1 and 2


########################################
#Plotting with fancy figure         ###
########################################

###### To highlight individuals by population
env <- read.csv("../Pop_assignment3.csv")

env <- env %>%
  filter(tissue_number != "F11-093") %>% 
  filter(tissue_number != "F11-156") %>%
  filter(tissue_number != "F09-11") %>%
  filter(tissue_number != "F09-45") %>%
  filter(tissue_number != "F10-44") %>%
  filter(tissue_number != "F10-154") %>%
  filter(tissue_number != "F09-27") %>% 
  filter(tissue_number != "F10-162") %>% 
  filter(tissue_number != "F09-141") %>%
  filter(tissue_number != "F11-170") 


eco <- env %>% 
  mutate(POP = recode(POP,"E-O"="Eastern Oak","W-O"="Western Oak","C-O"="Central Oak",
                      "E-P"="Eastern Pines","W-P"="Western Pines","C-P"="Central Pines"))

myCol2 <- c("darkorange","deepskyblue1", "yellow", "lightcyan2","orangered2","navy")

#Plot the individuals on axis 1 &2 

plot(issj.rda, type="n", scaling=3)
points(issj.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(issj.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=myCol2[eco$POP]) # the jays
text(issj.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(eco$POP), bty="n", col="gray32", pch=21, cex=1, pt.bg=myCol2)

# axes 1 & 3
plot(issj.rda, type="n", scaling=3, choices=c(1,3))
points(issj.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(issj.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=myCol2[eco$POP], choices=c(1,3))
text(issj.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topleft", legend=levels(eco$POP), bty="n", col="gray32", pch=21, cex=1, pt.bg=myCol2)

#########################################
### b) Identify RDA candidates        ###
#########################################

#Use the loadings of the SNPs (their location) in the ordination
#space to determine which SNPs are candidates for local adaptation.
#The SNP loadings are stored as `species` in the RDA object. We'll extract the SNP loadings from the first three constrained axes:

load.rda <- scores(issj.rda, choices=c(1:3), display="species")

#If we look at histograms of the loadings on each RDA axis,
#we can see their (relatively normal) distribution. SNPs loading at the center of the distribution are not showing
#a relationship with the environmental predictors; those loading in the tails are, and are more likely to be under 
#selection as a function of those predictors (or some other predictor correlated with them).


hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 


#define the function here as `outliers`, where `x` is the vector of loadings and `z` is the number of standard deviations to use:

#find loadings +/-z sd from mean loading
#locus names of the tails

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)       
  x[x < lims[1] | x > lims[2]]           
}

#Now let's apply it to the first three constrained axes using sd of 3.5 (strong selection)

# cand1 <- outliers(load.rda[,1],3.5) # 11
# cand2 <- outliers(load.rda[,2],3.5) # 1
# cand3 <- outliers(load.rda[,3],3.5) # 9
# 
# ncand <- length(cand1)+length(cand2)+length(cand3)
# ncand  ##22 candidates


#using an sd of 3 

cand1 <- outliers(load.rda[,1],3) # 35
cand2 <- outliers(load.rda[,2],3) # 54
cand3 <- outliers(load.rda[,3],3) # 37



ncand <- length(cand1)+length(cand2)+length(cand3)
ncand

#We have 164 potential candidate SNPs under moderate selection

#make a single data frame with the axis, SNP name, & loading:

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2)<- colnames(cand3) <- c("axis","snp","loading")


cand <- rbind(cand1, cand2, cand3)

cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=2)  # 3 columns for 3 predictors

colnames(foo) <- c("percent_pine", "percent_oak")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- genfile[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

dim(cand)

#Some of these may be duplicate detections; let's check:


# duplicate detections:
length(cand$snp[duplicated(cand$snp)])
foo <- cbind(cand$axis,duplicated(cand$snp))
table(foo[foo[,1]==1,2])
table(foo[foo[,1]==2,2])
table(foo[foo[,1]==3,2])

cand <- cand[!duplicated(cand$snp),] #remove duplicates

dim(cand)

## 164 unique outlier loci of 5768 SNPs

#which predictor each snp is more correlated with

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,6] <- names(which.max(abs(bar[4:5]))) # gives the variable
  cand[i,7] <- max(abs(bar[4:5]))              # gives the correlation
}

colnames(cand)[6] <- "predictor"
colnames(cand)[7] <- "correlation"

table(cand$predictor)
#           percent_oak           percent_pine 
#                64                   100                

##To highlight SNPs
sel <- cand$snp
env <-cand$predictor

env[env=="percent_pine"] <- '#33a02c'
env[env=="percent_oak"] <- '#e3aa1a'


# color by predictor:
col.pred <- rownames(issj.rda$CCA$v) # pull the SNP names


for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

colors <- c("#33a02c","#e3aa1a")

col.pred[!grepl(paste(colors,collapse="|"), col.pred)] <- '#f1eef6' # non-candidate SNPs


empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent

empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#33a02c','#e3aa1a')


# axes 1 & 2
plot(issj.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(issj.rda, display="species", pch=21, cex=1, col="grey32", bg=col.pred, scaling=3)
points(issj.rda, display="species", pch=21, cex=1,col=empty.outline, bg=empty, scaling=3)
text(issj.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("topright",legend=c("Percent Pine", "Percent Oak"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

#axes 1 &3 
plot(issj.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(issj.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(issj.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(issj.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topright",legend=c("Percent Pine", "Percent Oak"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

#edit the last line
#legend("bottomright", legend=c("Distance From Pine", "Mean Precip","Mean Temp"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)






###############################################################
# RDA with Morphological charactaristics 
# AKA GWAS!

#############################################################
#Convert genomic data. Using the .raw file from Plink 

gen <- read.PLINK("C:/Users/Rebecca/Dropbox/ISSJ/issj_genomics/r_issj_genomics/ISSJ_clean.raw", parallel = F)
dim(gen)

# RDA requires complete data frames (i.e., no missing genetic data). 
#For this example, we'll use a simple approach to imputing missing genotype values: 
#we will impute using the most common genotype at each SNP across all individuals.

gen <- as.matrix(gen)

sum(is.na(gen)) 

gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), 
                                             as.numeric(names(which.max(table(x))))))

sum(is.na(gen.imp)) # No NAs

genfile <- data.frame(gen.imp)
setDT(genfile, keep.rownames = T)


#Read in phenotypic data
env <- read.csv("../Pop_assignment3.csv")

#Filter to the clean genetic dataset   
pheno <- env %>%
  filter(tissue_number != "F11-093") %>% 
  filter(tissue_number != "F11-156") %>%
  filter(tissue_number != "F09-11") %>%
  filter(tissue_number != "F09-45") %>%
  filter(tissue_number != "F10-44") %>%
  filter(tissue_number != "F10-154") %>%
  filter(tissue_number != "F09-27") %>% 
  filter(tissue_number != "F10-162") %>% 
  filter(tissue_number != "F09-141") %>%
  filter(tissue_number != "F11-170") 

#Do any birds have missing bill length data?
pheno[is.na(pheno$bill_len_1),]

#Two: F10-89 (Pelican; East Pine), and F09-10 (Coches; Central Oak). Filter these out of Genotype file 

genfile <- genfile %>% 
  filter(rn != "F10-89") %>% 
  filter(rn !="F09-10") 


#Filter out individuals with missing bill length nares
pheno <-pheno[!is.na(pheno$bill_len_1),]

#Calculate predicted bill length measurements based on body size using formula interface of prcomp. *** Have to use 
#wing measurement because there are missing values of tarsus 

pca_wing <- prcomp(~ +wing_right +tarsus1, data=pheno, center=T, scale= T)

axes_wing <- predict(pca_wing, newdata = pheno)


##glm of length with body size correction
glm_pheno <- cbind.data.frame(pheno, axes_wing)

length_pheno <- lm(bill_len_1~habitat + PC1, data=glm_pheno) 

pheno$predict.bill.length <- predict(length_pheno)

#glm of depth with body size correction
# 
# glm_depth <- cbind.data.frame(pheno, axes_wing)
# 
# depth_pheno <- lm(bill_depth~+PC1, data=glm_depth) 
# 
# pheno$predict.bill.depth <- predict(depth_pheno)

#select for just predicted bill length or depth based on body size
#pheno <- pheno %>%   ###Don't trust this since it treats depth and length as seperate predictors when really they interact
#  dplyr::select(c("predict.bill.length","predict.bill.depth" ))


# Confirm that genotypes and phenotypes are in the same order
pheno <-pheno[order(match(pheno[,1],genfile[,1])),]

colnames(genfile)[1] <- "tissue_number"

#make sure pheno tissue number is a character class like the genotype dataset 
pheno$tissue_number <- as.character(pheno$tissue_number)

identical(genfile$tissue_number, pheno$tissue_number) 


#Filter to just the genotype info
genfile <- as.matrix(genfile [,-1], row.names =FALSE, col.names=FALSE)

mode(genfile) <- 'numeric'  #needs to be numeric for some reason and matrix was character


#filter to just the pheotype info
pheno <- pheno %>% 
  dplyr::select(c("predict.bill.length"))

#OR try checking the corelation fo the raw bill length and depth

# pheno <- pheno %>% 
#   dplyr::select(c("bill_lengt", "bill_depth"))


#check correlations of predictors
#pairs.panels(pheno[,1:2], scale=T)
#temp and precip are 50% correlated... that makes sense. Still less than .7


#Run the RDA (GWAS)
issj.rda <- vegan::rda(genfile~., data=pheno, scale=T)
issj.rda

vegan::RsquareAdj(issj.rda) #r2= 0.008, adusted r2= 0.00085

summary(eigenvals(issj.rda, model="constrained"))

#We can visualize this information using a screeplot of the canonical eigenvalues by calling `screeplot`:

screeplot(issj.rda) #not helpful since its just the one measure


signif.full <- anova.cca(issj.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

signif.axis <- anova.cca(issj.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
#RDA 1 is significant

#You can run a formal test of statistical significance of each constrained axis using:
#anova.cca(issj.rda, by="axis") #We can assess both the full model and each constrained axis using F-statistics
# (Legendre et al, 2010). The null hypothesis is that no linear relationship exists between the SNP data and the
#'# environmental predictors. See `?anova.cca` for more details and options.


#The permutation process to test the signficiance of each axis takes a while
#(up to a few hours on large data sets), so we'll just use the screeplot for a first assessment.
#If we did run the formal test, we would find that the first three constrained axes are significant (p = 0.001)

#Finally, `vegan` has a simple function for checking Variance Inflation Factors for the predictor variables used in the model:

vegan::vif.cca(issj.rda) 

#bill length is 1, but there are no other predictors to worry about coliniarity  

#quick plot of the RDA output using the default plotting in `vegan`

plot(issj.rda, scaling=3)          # default is axes 1 and 2


########################################
#Plotting with fancy figure         ###
########################################

###### To highlight individuals by population
env <- read.csv("../Pop_assignment3.csv")

env <- env %>%
  filter(tissue_number != "F11-093") %>% 
  filter(tissue_number != "F11-156") %>%
  filter(tissue_number != "F09-11") %>%
  filter(tissue_number != "F09-45") %>%
  filter(tissue_number != "F10-44") %>%
  filter(tissue_number != "F10-154") %>%
  filter(tissue_number != "F09-27") %>% 
  filter(tissue_number != "F10-162") %>% 
  filter(tissue_number != "F09-141") %>%
  filter(tissue_number != "F11-170") %>% 
  filter(tissue_number !="F09-11") %>% 
  filter(tissue_number !="F10-12") %>% 
  filter(tissue_number != "F10-89") %>% 
  filter(tissue_number !="F09-10") 


eco <- env %>% 
  dplyr::mutate(POP = dplyr::recode(POP,"E-O"="Eastern Oak","W-O"="Western Oak","C-O"="Central Oak",
                      "E-P"="Eastern Pines","W-P"="Western Pines","C-P"="Central Pines"))

myCol2 <- c("darkorange","deepskyblue1", "yellow", "lightcyan2","orangered2","navy")

#Plot the individuals on axis 1 &2 

plot(issj.rda, type="n", scaling=3)
points(issj.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(issj.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=myCol2[eco$POP]) # the jays
text(issj.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=levels(eco$POP), bty="n", col="gray32", pch=21, cex=1, pt.bg=myCol2)

# axes 1 & 3
plot(issj.rda, type="n", scaling=3, choices=c(1,3))
points(issj.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(issj.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=myCol2[eco$POP], choices=c(1,3))
text(issj.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topleft", legend=levels(eco$POP), bty="n", col="gray32", pch=21, cex=1, pt.bg=myCol2)

#########################################
### b) Identify RDA candidates        ###
#########################################

#Use the loadings of the SNPs (their location) in the ordination
#space to determine which SNPs are candidates for local adaptation.
#The SNP loadings are stored as `species` in the RDA object. We'll extract the SNP loadings from the first three constrained axes:

load.rda <- scores(issj.rda, choices=c(1:3), display="species")

#If we look at histograms of the loadings on each RDA axis,
#we can see their (relatively normal) distribution. SNPs loading at the center of the distribution are not showing
#a relationship with the environmental predictors; those loading in the tails are, and are more likely to be under 
#selection as a function of those predictors (or some other predictor correlated with them).


hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 


#define the function here as `outliers`, where `x` is the vector of loadings and `z` is the number of standard deviations to use:

#find loadings +/-z sd from mean loading
#locus names of the tails

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)       
  x[x < lims[1] | x > lims[2]]           
}

#Now let's apply it to the first three constrained axes using sd of 3.5 (strong selection)

#  cand1 <- outliers(load.rda[,1],3.5) 
#  cand2 <- outliers(load.rda[,2],3.5) 
#  cand3 <- outliers(load.rda[,3],3.5) 
# # 
#  ncand <- length(cand1)+length(cand2)+length(cand3)
#  ncand  ##14 candidates under strong selection


#using an sd of 3.0 since a majority of adaptation is polygenic and therefore under weak to moderate selection (Hendry's 2019 book)
# 
cand1<- outliers(load.rda[,1],3.0) 
cand2 <- outliers(load.rda[,2],3.0) 
cand3 <- outliers(load.rda[,3],3.0) 



ncand <- length(cand1)+length(cand2)+length(cand3)
ncand

#We have 513 potential candidate SNPs under moderate selection

#make a single data frame with the axis, SNP name, & loading:

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2)<- colnames(cand3) <- c("axis","snp","loading")


cand <- rbind(cand1, cand2, cand3)

cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=1)  # 1 columns for 1 predictors

colnames(foo) <- c("bill length")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- genfile[,nam]
  foo[i,] <- apply(pheno,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

dim(cand)

#Some of these may be duplicate detections; let's check:


# duplicate detections:
length(cand$snp[duplicated(cand$snp)])
foo <- cbind(cand$axis,duplicated(cand$snp))
table(foo[foo[,1]==1,2])
table(foo[foo[,1]==2,2])
table(foo[foo[,1]==3,2])

cand <- cand[!duplicated(cand$snp),] #remove duplicates

dim(cand)

## 101 unique outlier loci of 5768 SNPs

#which predictor each snp is more correlated with

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,5] <- names(which.max(abs(bar[4]))) # gives the variable
  cand[i,6] <- max(abs(bar[4]))              # gives the correlation
}

colnames(cand)[5] <- "predictor"
colnames(cand)[6] <- "correlation"

table(cand$predictor)

##To highlight SNPs
sel <- cand$snp
env <-cand$predictor
env[env=="bill length"] <- '#F9E424'




# color by predictor:
col.pred <- rownames(issj.rda$CCA$v) # pull the SNP names


for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

colors <- c("#F9E424")

col.pred[!grepl(paste(colors,collapse="|"), col.pred)] <- '#f1eef6' # non-candidate SNPs

empty <- col.pred

empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")

#Color the length morph SNPs
bg <- c('#F9E424')


# axes 1 & 2
plot(issj.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(issj.rda, display="species", pch=21, cex=1, col="grey32", bg=col.pred, scaling=3)
points(issj.rda, display="species", pch=21, cex=1,col=empty.outline, bg=empty, scaling=3)
text(issj.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("topright",legend=c("Bill Length"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

#axes 1 &3 
plot(issj.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(issj.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(issj.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(issj.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topright",legend=c("Bill Length"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)




###########################
# Overlapping Snps      ###
###########################


#read in results from univaritate GWAS from TASSEL
ISSJgwas <- read.csv("ISSJ_MLM_results.csv")



#Change cand to dataframe 
RDAcand <- as.data.frame(cand)


#Split cand so the snp name matches the output from TASSEL
#make a new df
foo <- data.frame(do.call('rbind', strsplit(as.character(RDAcand$snp),'_',fixed=TRUE)))

#get rid of the X in the snp name 
foo$X1 <- str_remove(foo$X1, "X")

#merge the two snp number collumns 
foo$snp <-paste(foo$X1, foo$X2, sep="_")

#Select just the snp name and replace the snp name from RDA
foo <- foo %>% 
  dplyr::select(c("snp"))

foo$snp <- as.character(foo$snp)

RDAcand$snp <- foo$snp

#make sure the classes match of the two snp datafrmes 

RDAcand$snp <- as.character(RDAcand$snp)
ISSJgwas$Marker <- as.character(ISSJgwas$Marker)

colnames(ISSJgwas)[2]<- "snp"

#filter to just the trait of interest from TASSEL
TASSEL <- ISSJgwas %>% 
  subset(Trait=="nares_predict")

#count how many matches
sharedsnps<- merge(TASSEL, RDAcand, by = "snp", all = FALSE)


sharedsnps <- sharedsnps[!duplicated(sharedsnps$snp),] #remove duplicates

dim(sharedsnps) #90 shared SNPs between tests 

#check to see how many are shared between environment RDA and GWAS dataset 




######################################################
# Overlapping Snps GWAS and environmental RDA      ###
######################################################


#Change ecand to dataframe 
env.RDA <- as.data.frame(ecand)


#Split cand so the snp name matches the output from TASSEL
#make a new df
envfoo <- data.frame(do.call('rbind', strsplit(as.character(env.RDA$snp),'_',fixed=TRUE)))

#get rid of the X in the snp name 
envfoo$X1 <- str_remove(envfoo$X1, "X")

#merge the two snp number collumns 
envfoo$snp <-paste(envfoo$X1, envfoo$X2, sep="_")

#Select just the snp name and replace the snp name from RDA
envfoo <- envfoo %>% 
  dplyr::select(c("snp"))

envfoo$snp <- as.character(envfoo$snp)

env.RDA$snp <- envfoo$snp

#make sure the classes match of the two snp datafrmes 

env.RDA$snp <- as.character(env.RDA$snp)
sharedsnps$snp <- as.character(sharedsnps$snp)


#count how many matches
shared.all<- merge(env.RDA, sharedsnps, by = "snp", all = FALSE)

shared.all <- shared.all[!duplicated(shared.all$snp),] #remove duplicates

dim(shared.all) #1 shared SNPs between env and pheno tests. Most strongly associated with temperature followed by pine 


########################################################
###                     Manhattan Plot               ###
########################################################


#Note that Zang et al. 2019 advocate for a P-value of P = 0.0002,
#or a LOD score of 3.0 in GWAS to void erroneous outlier detection
#But this is for a multi-locus test. Whereas TASSEL or GEMMA are single Locus tests


#use the shared SNPs from the univariate TASSLE GWAS and Multivariate RDA

#Find outliers 
alpha <- 0.05
outliersnares <-sharedsnps %>% 
  filter(sharedsnps$p<alpha)

SNPslist <- paste(outliersnares$snp, collapse = ", ") #get snp list
#SNPslist

#Mark the outlier SNPs seperately from the shared SNPs
InterestingSNPsnares <- c("104266_60", "105207_49", 
                          "108866_30", "110742_8", "113215_92", "139089_94", 
                          "14782_11", "148932_23", "16165_45", "164274_25",
                          "22010_19", "34696_23", "62081_20", "654_56")


##mutate the colnames so they are recognizable by the package of the full GWAS SNP and position
#data from TASSEL
#Omit missing values in the position collumn
ISSJgwas <-ISSJgwas[!is.na(ISSJgwas$Pos),]

ISSJgwas <- ISSJgwas %>%
  mutate(str_remove(ISSJgwas$Chr,"SCAFFOLD_")) %>% 
  dplyr::rename (
    "SNP"="snp",
    "CHR"="str_remove(ISSJgwas$Chr, \"SCAFFOLD_\")",
    "P"= "p"
  )

##Select just the collumns of interest, the SNP, Chromosome (in this case scaffold), and position 
ISSJgwas <- ISSJgwas %>% 
  dplyr::select("SNP", "CHR", "P" ) 

ISSJgwas$BP <- seq.int(nrow(ISSJgwas))

ISSJgwas$CHR <- as.numeric(as.character(ISSJgwas$CHR))
ISSJgwas <- na.omit(ISSJgwas)

ISSJgwas <- ISSJgwas[!duplicated(ISSJgwas$SNP),] #remove duplicates


## Manhattan plot highlighting the SNPS of interest**Add highlight= Interesting SNPs Nares
#something weird going on where multiple SNPS are being highlighted twice. Sometime above the threshold, and sometimes not

manplot <- qqman::manhattan(ISSJgwas, main="P-values by Scaffold for Corrected Bill Length", xlab="Scaffold", 
                     suggestiveline = 1.3,col = c("lightgrey", "black"), genomewideline = F, highlight = InterestingSNPsnares, )


#If you want to change the highlight color you have to edit the source code > trace(qqman::manhattan, edit = T)

