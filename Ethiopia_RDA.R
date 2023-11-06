## R Script: Genotype-environment association analysis with RDA ##
## Written by Guillermo Friis, last edit October 22, 2023; R-4.2.2 ##
## Edited J Borrell, 3/11/23

setwd('C:/Users/jb83kg/Documents/Ethiopia_agrisystems_paper/RDA_run_test')

library(readxl)
library(raster)
library(sdmpredictors)
library(psych)
library(vroom)
library(caret)
library(vegan)
library(RPMG)
library(robust)
library(pcadapt)
library(qvalue)
library(dplyr)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")

# Extracting environmental variables - Coffee example
#---------------------------------------------------------------------------------------------------------------------------------------------
## Sample data and coordinates - Ignoring phenotype and soil data for now
agris.df <- read_excel('Agrisystems_phenotype_data_pasture.xlsx',sheet = 1, col_names = T, na = 'NA')
agris.df[agris.df == 'Na'] <- NA

species.df <- droplevels(agris.df[agris.df$Sp == 'Mai' & agris.df$Concern == 'OK', ])
species.df <- species.df[ , c(7, 20, 19, 11, 12)]
species.df$Pop <- paste(species.df$Transect, species.df$Site, sep = '_')
colnames(species.df)[1:2] <- c('Sample', 'Lon')

table(species.df$Transect, species.df$Site)

selected_variables <- read.csv("~/Ethiopia_agrisystems_paper/RDA_run_test/selected_variables.txt", sep="")


# Genotypes 
##---------------------------------------------------------------------------------------------------------------------------------------------
## Full Dataset - vcf needs to be previously converted to 012 format, easy with vcftools
## There are several inconsistencies between the phenotype table and the vcf file

library(dartR)
library(vcfR)
vcf <- read.vcfR("C:/Users/jb83kg/Documents/Ethiopia_agrisystems_paper/RDA_run_test/maize.postfilt.noExcluded.noIndels.vcf/maize.postfilt.noExcluded.noIndels.vcf", verbose = FALSE )
my_genclone <- vcfR2genlight(vcf)
my_genclone@ind.names
my_genclone@ind.names

my_genclone@gen
gl.alf(my_genclone)

Dat1 <- my_genclone[indNames(my_genclone) != species.df$Sample]
Dat1@pop <- as.factor(species.df$Pop)

?gl2geno
genind2df(Dat1@gen)


snps.dataset <- vroom('coffee_ethio.012', delim = '/t', col_names = F, na = c("-1"), progress = T, num_threads = 4)
snps.ind <- read.table('coffee_ethio.012.indv', sep = '/t', header = F)
snps.dataset <- cbind(snps.ind$V1, snps.dataset[ , -1])
colnames(snps.dataset) <- paste('SNP', seq(0, ncol(snps.dataset)-1), sep = '')


snps.dataset <- merge(coffee.table[ , c(1, 6)], snps.dataset, by.x = 'Sample', by.y = 'SNP0')
snps.dataset$Pop <- as.factor(snps.dataset$Pop)
head(snps.dataset[, 1:10])

## Some informative statistics, they may help to set filters later for each species
snps_stats <- function(x) {

  a <- table(x$Pop)
  b <- table(x$Pop)[which.min(table(x$Pop))]
  c <- (sum(is.na(x))/(nrow(x)*ncol(x)-2))*100
  d <- cbind(x$Sample, (rowSums(is.na(x))/(ncol(x)-2))*100)
  
  aux0.df <- table(snps.dataset$Pop)
  
  for (i in levels(snps.dataset$Pop)){
    
    aux1.df <- snps.dataset[snps.dataset$Pop == i, ]
    aux1.df <- aux1.df[,-c(1, 2)]
    aux2.df <- rowSums(is.na(aux1.df))/(ncol(aux1.df))
    aux0.df[i] <- mean(aux2.df)*100
    
  }
  
  print('Samples per site:')
  print(a)
  print('Site with fewer samples:')
  print(b)
  print('Percentage of total missing data:')
  print(c)
  print('Mean missing data per sampling site:')
  print(aux0.df)  
  print('Missing data per sample:')
  noquote(d)

}

snps_stats(snps.dataset)

## Compute SNP frequencies
# I created a function to filter SNP datasets and compute frequencies, the parameters are:
# x: The SNP dataset AS CREATED ABOVE
# nmin: The minimum number of individuals PER SITE for the site to be included in the final dataset
# ploidy: Ploidy number
# missing: The maximum ALLOWED RATE of missing data per SNP position
# by.ind: If True, the missing data threshold is applied at the entire dataset, and potential per site missing frequencies are imputed.
#         If False, the threshold is applied by sampling site (more restrictive), which render no missing frequencies.

snps_freqs <- function(x, nmin, ploidy, missing, by.ind = T) {
  
  if (by.ind){
    
    x <- x[ ,colSums(is.na(x))/nrow(x) <= missing]
    snps.mean <- as.data.frame(matrix(ncol = (ncol(x)-2)))
    colnames(snps.mean) <- paste('SNP', seq(1, (ncol(x)-2)), sep = '')
    
    h = 0
    for (i in levels(x$Pop)){
      
      aux.df <- x[x$Pop == i, ]
      
      if (nrow(aux.df) >= nmin) {
        
        h = h + 1
        aux.df <- aux.df[,-c(1, 2)]
        mean.df <- colMeans(aux.df)
        snps.mean[h, ] <- mean.df
        rownames(snps.mean)[h] <- i
        
      }
      
    }
    
    snps.dataset.freq <- snps.mean/ploidy
    snps.dataset.freq <- apply(snps.dataset.freq, 2, function(y){
      y[is.na(y)] <- names(which.max(table(y)))
      return(y) })
    
    snps.dataset.freq <- as.data.frame(snps.dataset.freq)
    snps.dataset.freq <- snps.dataset.freq %>% mutate_at(vars(starts_with("SNP")), as.numeric)
    
  } else {
    
    snps.mean <- as.data.frame(matrix(ncol = (ncol(x)-2)))
    colnames(snps.mean) <- paste('SNP', seq(1, (ncol(x)-2)), sep = '')
    
    h = 0
    mis.vec <- vector()
    
    for (i in levels(x$Pop)){
      
      aux.df <- x[x$Pop == i, ]
      
      if (nrow(aux.df) >= nmin) {
        
        h = h + 1
        aux.df <- aux.df[,-c(1, 2)]
        
        mis.vec <- c(mis.vec, names(which(colSums(is.na(aux.df))/nrow(aux.df) > missing)))
        
        mean.df <- colMeans(aux.df, na.rm = T)
        snps.mean[h, ] <- mean.df
        rownames(snps.mean)[h] <- i
      }
      
      
    }  
    
    mis.vec <- unique(mis.vec)
    snps.dataset.freq <- snps.mean/ploidy
    snps.dataset.freq <- snps.dataset.freq[ , !(names(snps.dataset.freq) %in% mis.vec)]  
    
  }
  
  snps.dataset.freq <- snps.dataset.freq[, sapply(snps.dataset.freq, function(x) length(unique(snps.dataset.freq)) > 1)]
  return(snps.dataset.freq)
  
}
  
snps.dataset.freq <- snps_freqs(snps.dataset, nmin = 3, ploidy = 2, missing = 0.5, by.ind = F)
sum(is.na(snps.dataset.freq))
head(snps.dataset.freq[, 1:10])


# RDA
##---------------------------------------------------------------------------------------------------------------------------------------------
## Validating datasets for RDA
coffee.set <- coffee.set[ row.names(coffee.set) %in% row.names(snps.dataset.freq), ]
identical(row.names(coffee.set), row.names(snps.dataset.freq))

## Simple RDA
# Forward selection method
#rda0 <- rda(snps.dataset.freq ~1, coffee.set)
#rda.full <- rda(snps.dataset.freq ~., coffee.set)
#anova(rda.full)
#adjr2.rdafull <- RsquareAdj(rda.full)$adj.r.squared

#mod <- ordiR2step (rda0, scope = formula (rda.full), R2scope = adjr2.rdafull, direction = 'forward', permutations = 999)

#mod.anova <- mod$anova
#rm(rda0, rda.full, mod, adjr2.rdafull)

# VIF
rda.clim <- rda(snps.dataset.freq ~., coffee.set)
vif.cca(rda.clim) # Rerun until VIF < 5
rda.clim <- rda(snps.dataset.freq ~
                CHELSA_bio10_07 +
                CHELSA_bio10_10 +
                CHELSA_bio10_13 +
                Aridity +
                Roughness +
                Density +
                Landscape +
                GAEZ,
                coffee.set)

screeplot(rda.clim)
RsquareAdj(rda.clim)
clim.rda.summary <- summary(rda.clim)
capture.output(clim.rda.summary, file ='RDA_SummaryReport.txt')
rm(clim.rda.summary)

# RDA plot
transects.df <- unique(coffee.df[, 4:5])
transects.df <- aggregate(Site ~ Transect, data = transects.df, FUN = max)

transects.df <- cbind(transects.df, data.frame(
  c('Dark Red', 'Dark Blue', 'Dark Green', 'Purple', 'Dark Orange', 'Dark Cyan', 'Dark Magenta', 'Yellow'),
  c('Red', 'Blue', 'Green', 'Lavender', 'Orange', 'Cyan', 'Magenta', 'Light Yellow')))

colnames(transects.df)[3:4] <- c('Color1', 'Color2')
        
generate_gradient <- function(start_color, end_color, num_shades) {
  gradient_palette <- colorRampPalette(c(start_color, end_color))
  return(gradient_palette(num_shades))
}

concatenated_gradient <- character(0)

for (i in 1:nrow(transects.df)) {
  start_color <- transects.df[i, 3]
  end_color <- transects.df[i, 4]
  num_shades <- as.numeric(transects.df[i, 2])
  
  row_gradient <- generate_gradient(start_color, end_color, num_shades)
  concatenated_gradient <- c(concatenated_gradient, row_gradient)
}


plot(rda.clim, scaling=3, pch = 21, bg = concatenated_gradient, cex = 1.3)
points(rda.clim, display="sites", pch = 21,
       bg = concatenated_gradient, scaling = 3, cex = 1.3, choices = c(1, 2))
legend("topright", legend = paste(transects.df$Transect, transects.df$Site, sep = ' - '),
       bty="n", pch = 21, pt.bg = transects.df$Color1)


## variance Partition
spe.part <- varpart(snps.dataset.freq,
                    coffee.set[ , c('CHELSA_bio10_07', 'CHELSA_bio10_10', 'CHELSA_bio10_13', 'Aridity', 'Roughness', 'GAEZ')],
                    coffee.set[ , c('Density', 'Landscape')])
spe.part
plot(spe.part)

# Test of fractions [a+b+c]
permutation_sig.all <- anova.cca(rda.clim, step = 10000)

# Test of fractions [a+b]
anova.clim <- anova.cca(rda(snps.dataset.freq ~
                              CHELSA_bio10_07 +
                              CHELSA_bio10_10 +
                              CHELSA_bio10_13 +
                              Aridity +
                              Roughness +
                              GAEZ, coffee.set), step = 10000)

# Test of fractions [b+c]
anova.soc <- anova.cca(rda(snps.dataset.freq ~
                              Density + Landscape, coffee.set), step = 10000)

# Test of fraction [A]
anova.a <- anova.cca(rda(snps.dataset.freq ~
                           CHELSA_bio10_07 +
                           CHELSA_bio10_10 +
                           CHELSA_bio10_13 +
                           Aridity +
                           Roughness +
                           GAEZ + Condition(Density + Landscape), coffee.set),
                           step = 10000)

# Test of fraction [C]
anova.c <- anova.cca(rda(snps.dataset.freq ~
                           Density + Landscape +
                           Condition(CHELSA_bio10_07 +
                                     CHELSA_bio10_10 +
                                     CHELSA_bio10_13 +
                                     Aridity +
                                     Roughness +
                                     GAEZ), coffee.set),
                                     step = 10000)

# Partial RDA
##---------------------------------------------------------------------------------------------------------------------------------------------
## Neutral dataset
path_to_vcf <- "F:/C_Kew/B_Papers/D_GEAeth/A_Datasets/coffee_ethioProjOnly.postfilt.noIndels.vcf"
snps.vcf <- read.pcadapt(path_to_vcf, type='vcf')

snps.pc <- pcadapt(snps.vcf, K = 4)
summary(snps.pc)

pop.vec <- agris.df[agris.df$`Full_id_v1 match` %in% snps.ind$V1, ]
pop.vec <- pop.vec[match(snps.ind$V1, pop.vec$`Full_id_v1 match`), ]
pop.vec <- pop.vec$Transect
plot(snps.pc,option="screeplot")
plot(snps.pc, option = "scores", pop = pop.vec)

hist(snps.pc$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(snps.pc, option = "stat.distribution")

# Outliers
qval <- qvalue(snps.pc$pvalues)$qvalues
outliers <- which(qval < 0.05)
length(outliers)

# Neutral dataset and PCA
snps.neu.freq <- snps.dataset.freq[,-(outliers)]
freq.pca <- rda(snps.neu.freq, scale=T)
plot(freq.pca, scaling=3, pch = 21, bg = concatenated_gradient, cex = 1.3)
points(freq.pca, display="sites", pch = 21,
       bg = concatenated_gradient, scaling = 3, cex = 1.3, choices = c(1, 2))
legend("topleft", legend = paste(transects.df$Transect, transects.df$Site, sep = ' - '),
       bty="n", pch = 21, pt.bg = transects.df$Color1)

pc.df <- as.data.frame(scores(freq.pca, choices=c(1:2), display="sites", scaling=0))

rm(snps.pc, path_to_vcf, snps.neu.freq, snps.vcf, qval, freq.pca, outliers)

## Partial RDA
rda.pclim <- rda(snps.dataset.freq ~
                   CHELSA_bio10_07 +
                   CHELSA_bio10_10 +
                   CHELSA_bio10_13 +
                   Aridity +
                   Roughness +
                   Density +
                   Landscape +
                   GAEZ+
                   Condition(pc.df$PC1 + pc.df$PC2),
                   coffee.set)

plot(rda.pclim, scaling=3, pch = 21, bg = concatenated_gradient, cex = 1.3)
points(rda.pclim, display="sites", pch = 21,
       bg = concatenated_gradient, scaling = 3, cex = 1.3, choices = c(1, 2))
legend("bottomright", legend = paste(transects.df$Transect, transects.df$Site, sep = ' - '),
       bty="n", pch = 21, pt.bg = transects.df$Color1)

## variance Partition
spe.part <- varpart(snps.dataset.freq,
                    coffee.set[ , c('CHELSA_bio10_07', 'CHELSA_bio10_10', 'CHELSA_bio10_13', 'Aridity', 'Roughness', 'GAEZ')],
                    coffee.set[ , c('Density', 'Landscape')],
                    pc.df)
spe.part
plot(spe.part)

# Test of fractions [ALL]
rda.clim_land_pc <- rda(snps.dataset.freq ~
                          CHELSA_bio10_07 +
                          CHELSA_bio10_10 +
                          CHELSA_bio10_13 +
                          Aridity +
                          Roughness +
                          Density +
                          Landscape +
                          GAEZ +
                          pc.df$PC1 +
                          pc.df$PC2,
                          coffee.set)

permutation_sig_all <- anova.cca(rda.clim_land_pc, step = 10000)


# Test of fraction [Clim]
prda.A <- rda(snps.dataset.freq ~
                  CHELSA_bio10_07 +
                  CHELSA_bio10_10 +
                  CHELSA_bio10_13 +
                  Aridity +
                  Roughness +
                  GAEZ +
                  Condition(Density + Landscape + pc.df$PC1 + pc.df$PC2),
                  coffee.set)

anova.A <- anova.cca(prda.A, step = 10000)

# Test of fraction [Socioeconomic]
prda.B <- rda(snps.dataset.freq ~
                   Density + Landscape +
                   Condition(CHELSA_bio10_07 +
                             CHELSA_bio10_10 +
                             CHELSA_bio10_13 +
                             Aridity +
                             Roughness +
                             GAEZ +
                             pc.df$PC1 +
                             pc.df$PC2),
                   coffee.set)

anova.B <- anova.cca(prda.B, step = 10000)

# Test of fraction [Neutral structure]
prda.C <- rda(snps.dataset.freq ~
                   pc.df$PC1 + pc.df$PC2 +
                   Condition(Density +
                             Landscape +
                             CHELSA_bio10_07 +
                             CHELSA_bio10_10 +
                             CHELSA_bio10_13 +
                             Aridity +
                             Roughness +
                             GAEZ),
                   coffee.set)

anova.C <- anova.cca(prda.C, step = 10000)
