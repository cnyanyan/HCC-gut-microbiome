library(readxl)
library(tidyr)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(dplyr)
library(vegan)
library(tibble)
library(ggpubr)
library(wesanderson)
library(tidyverse)
library(reshape2)
library(gplots)
library(pheatmap)
library(momr)

## Loading data
metadata = read.csv("~/Desktop/Project/data/metadata.csv",header=T,sep=",")
merged_final = read.csv("~/Desktop/Project/data/merged_final.csv",sep="\t",header=T)
gene_lengths = read.csv("~/Desktop/Project/data/gene_lengths.tab",sep="\t",header=T)
counting = read.csv("~/Desktop/Project/data/counting_report.csv",sep=",",header=T)
load('~/Desktop/Project/data/MmMGS_April2014.RData')

merged = merged_final
row.names(merged) = merged$gene_id
merged$gene_id = NULL
merged = data.matrix(merged)
merged = round(merged,0)

# The downsizing does not want to work with  #
# the gene names as row names - this is very #
# strange. We set the row names to numbers   #
# and then restore the row names after the   #
# downsizing to the gene names               #
geneRowNames = row.names(merged)
row.names(merged) = c(1:nrow(merged)) 

gene_len = gene_lengths
gene_len = subset(gene_len,gene_len$gene %in% row.names(merged))
gene_len = tibble::deframe(gene_len)

# DOWNSIZING
min_nb_reads = summary(colSums(merged))
min_nb_reads = min_nb_reads["Min."]
downSiz10M = downsizeMatrix(data = merged,
                            level = min_nb_reads,
                            repetitions = 1,
                            silent = F)
row.names(downSiz10M) = geneRowNames
# NORMALIZATION
Norm10m = normFreqRPKM(dat=downSiz10M, cat=gene_len)

mgs = row.names(MmMGS$attr)[grep('MmMGS',row.names(MmMGS$attr))]  # Get only MGS and not MmCAG
MGS = MmMGS$sets
MGS = MGS[mgs]
for(i in 1:length(MGS)){
  MGS[[i]] <- MGS[[i]][1:50]   # only the first 50 genes (enough information with them)
}

# get all the genes
mgsGeneList = unique(do.call(c, MGS)) 
# subset normalised data to only have the genes that are in the reduced catalog
data = subset(Norm10m,row.names(Norm10m) %in% mgsGeneList)
data[is.na(data)] <- 0
genebag = row.names(data)
mgs = projectOntoMGS(genebag=genebag, list.mgs=MGS)
mgs.dat = extractProfiles(mgs, data)
mgs.med.vect = computeFilteredVectors(profile=mgs.dat, type="median")
mgs.med.vect = mgs.med.vect[rowSums(mgs.med.vect)>0,]


