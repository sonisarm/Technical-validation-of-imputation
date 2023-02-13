``` r

##########################################################
### Author: Sarmiento Cabello, Sonia                   ###
### Version: 1.0.                                      ###
### Objective: Compare genotypes between two VCFs.     ###
##########################################################
``` 

This code compares the genetic information of two VCF files. It was designed to
assess the accuracy of imputed low-coverage sequences in comparison to high-coverage
sequences from the same individual. First, the distribution of differences is 
plotted for all replicates, then a closer examination is done on the differences
in specific sections of the genome to identify areas with more discrepancies.
IMPORTANT: Input GDS have to be previously filtered to only include intersect SNPs.

``` r

### Load libraries ### 
library(SNPRelate) 
library(stats)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# TO INSTALL : 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("SNPRelate")


# this is how you go from vcf to gds btw: 
#vcf <- 'your.vcf.name.vcf'
#snpgdsVCF2GDS(vcf,'my_file.gds')

### OPEN GDS, GET GENO & SAMPLE NAMES ### 
gds1 <- snpgdsOpen('LowCov_imputed_Super-Scaffold_14.gds') 
gds2 <- snpgdsOpen('HighCov_imputed_Super-Scaffold_14.gds')

SNPs_gds1 = snpgdsSNPList(gds1)
SNPs_gds2 = snpgdsSNPList(gds2)

geno1 <- snpgdsGetGeno(gds1)
geno2 <- snpgdsGetGeno(gds2)
sm <- snpgdsSummary(gds1)

closefn.gds(gds1)
closefn.gds(gds2)

rownames(geno1) <- sm$sample.id
rownames(geno2) <- sm$sample.id

# 1) Computing differences between the genotypes of the two VCFs (overall view) #
result <- matrix(nrow = 32)
mis <- matrix(nrow = 32)

for(i in 1:32){
  tmp <- geno1[i,] == geno2[i,] 
  result[i,] <- length(tmp[tmp==FALSE])
  mis[i,] <- length(tmp[is.na(tmp)])
  
}

# Change names and calculate percentages
colnames(result)[1] <- 'mismatch'
rownames(result) <- rownames(geno1)

colnames(mis)[1] <- 'missing genotypes'
rownames(mis) <- rownames(geno1)

result_percentage <- result/675410*100  # Here I would have to remove the number of missing per individual from the division 
missing_percentage <- mis/675410*100

# Plot
jpeg(paste0('Mismatches and missing genotypes distribution.jpg'), width=12,height=6,unit='in',quality = 1000,res=800)
par(mfrow=c(1,2))
# Plot distribution of mismatches
hist(result_percentage, main = 'Mismatched genotypes distribution in Super-Scaffold_14', xlab='Mismatching genotypes (%)', breaks = 40)
# Plot distribution of missing genotypes
hist(missing_percentage, main = 'Missing genotypes distribution in Super-Scaffold_14', xlab='Missing genotypes (%)', breaks = 40)
dev.off()

# Check if individuals having high mismatches also have high missingness
 # Prepare for merging
result <- as.data.frame(result)
result$indv <- rownames(result)
mis <- as.data.frame(mis)
mis$indv <- rownames(mis)

# merge by indv
df <- merge(result, mis)

``` 
Now, lets plot the differences along the scaffold for each individual. 
However, the output coming from this is too cluttered due to the high number of 
SNPs in the scaffold. To mitigate this, we established SNP windows, counted 
the differences within each window, and divided by the window size to calculate 
the proportion.

``` r

# 2) Computing differences between the genotypes of the two VCFs (per position) #
# Substract genotypes and square. 
geno_substract <- (geno1 - geno2)^2   # Mismatch is indicated with a 1.

### Define functions
  # Count mismatches
count.mismatch <- function(row){ length(which(row == 1))}

  # Create dataframe with counts
mismatch <- function(geno.chunk){
  return(apply(geno.chunk,1,FUN=count.mismatch))
}

# Count mismatches
count.NA <- function(row){ length(which(is.na(row)))}
       
# Create dataframe with counts
NonApp <- function(geno.chunk){
  return(apply(geno.chunk,1,FUN=count.NA))
}

### 

# Set SNP windows
win <- 1000
round.down <- ncol(geno1) - (ncol(geno1) %% win) # rounds down to previous 1000th
change <- (ncol(geno1) - (ncol(geno1) %% win)) / win
# makes dataset of intervals of SNPs 
int <- data.frame('start'=seq(1,round.down,win),
                  'end'=seq(1,round.down,win)+win-1)

# Get total no. of snps and positions (should be the same in both VCFs/GDS)
maxsnps <- nrow(SNPs_gds1)
pos <- SNPs_gds1$position

#empty dataframes with 1 row / window SNP interval and 1 col / sample 
mdf <- as.data.frame(matrix(NA,nrow=nrow(int), ncol=nrow(geno_substract)))
NonAppl <- as.data.frame(matrix(NA,nrow=nrow(int), ncol=nrow(geno_substract)))

# Loop through windows obtaining no. of NAs (NonAppl) and mismatches (mdf)
for(i in 1:nrow(int)){
  # Sum number of NAs
  NonAppl[i,] <-  NonApp(geno_substract[,seq(int$start[i],int$end[i])])
  # Sum number of mismatches
  mdf[i,] <- mismatch(geno_substract[,seq(int$start[i],int$end[i])])
  # Extract mean position of SNPs in the window
  pos[i] <- mean(SNPs_gds1$position[int$start[i]:int$end[i]])
}

# To make a percentage, multiply by 100 and divide by window size (minus NAs)
# NAs have not been accounted for (it's not a match) and could bias our results
mdf_final <- mdf * 100/(win-NonAppl)

names(mdf_final) <- sm$sample.id
head(mdf_final)


# Plot mismatches in snps windows:
jpeg(paste0('Mismatches between high- and imputed low-coverage.jpg'), width=12,height=6,unit='in',quality = 1000,res=800)
par(mfrow=c(1,1))
plot(0,pch='',xlab=paste0(win,'-snps window index'),ylab='% mismatching sites', xlim=c(0,nrow(int)),
     ylim=c(0,max(mdf_final)),sub='ss14', main = "Mismatches between high- and imputed low-coverage")
for(i in 1:ncol(mdf_final)){
  lines(mdf_final[,i])
}
dev.off()

``` 

In my dataset, I found a clear distinction between indidividuals with high (15%) 
and low (5%) proportion of mismatches. It is meaningful to plot the individuals 
separately and check if we see any differences.

``` r

# Get list of individuals
result_percentage <- as.data.frame(result_percentage)
topindv <- result_percentage %>% filter(result_percentage$mismatch > 14)
topindv <- rownames(topindv)
bottomindv <- result_percentage %>% filter(mismatch < 14)
bottomindv <- rownames(bottomindv)

# Get mismatching values accordingly
mdf_top <- mdf_final[, topindv]
mdf_bottom <- mdf_final[, bottomindv]

# Set colors with RColorBrewer package to differentiate between individuals
n <- ncol(mdf_top)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

# Function
# transparency color function
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Plot separating by individuals
jpeg(paste0('Mismatches between high- and imputed low-coverage - colors - separated indvs.jpg'), width=12,height=6,unit='in',quality = 1000,res=800)
par(mfrow=c(2,1))

plot(0,pch='',xlab=paste0(win,'-snps window index (ss14)'),ylab='% mismatches',
     xlim=c(0,nrow(int)), col = col_vector,
     ylim=c(0,max(mdf_top)), 
     main = "Mismatches between high- and imputed low-coverage (indv. with high missing rates)",
     cex.lab=0.7, cex.axis=0.7, cex.main=1)

for(i in 1:ncol(mdf_top)){
  lines(mdf_top[,i], col=add.alpha(col_vector[i],.8))
}
legend('topleft',col=col_vector,lwd=c(2,2),legend=topindv, cex = 0.25)

plot(0,pch='',xlab=paste0(win,'-snps window index (ss14)'),ylab='% mismatches',
     xlim=c(0,nrow(int)), col = col_vector,
     ylim=c(0,max(mdf_bottom)), 
     main = "Mismatches between high- and imputed low-coverage (indv. with low missing rates)",
     cex.lab=0.7, cex.axis=0.7, cex.main=1)
for(i in 1:ncol(mdf_bottom)){
  lines(mdf_bottom[,i], col=add.alpha(col_vector[i],.8))
}
legend('topleft',col=col_vector,lwd=c(2,2),legend=bottomindv, cex = 0.25)
dev.off()

```
