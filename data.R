
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
library(knitr)

## Loading data
merged = read.csv("~/Desktop/Project/data/merged.csv",sep="\t",header=T)
metadata <- as.data.frame(read.csv("~/Desktop/Project/data/metadata.csv",header=T,sep=","))
merged_final = as.data.frame(read.csv("~/Desktop/Project/data/merged_final.csv",sep="\t",header=T))
gene_lengths = as.data.frame(read.csv("~/Desktop/Project/data/gene_lengths.tab",sep="\t",header=T))
counting <- as.data.frame(read.csv("~/Desktop/Project/data/counting_report.csv",sep=",",header=T))
mgs_species = as.data.frame(read.csv("~/Desktop/Project/data/mgs_species.csv",sep=",",header=T))
mgs.med = read.csv('~/Desktop/Project/data/mgs.csv')
taxo = read.csv('~/Desktop/Project/data/taxo.csv')



mgs.med1 <- mgs.med
for (i in 1:26){
  mgs.med1[,i+1] <- mgs.med[,index[i]+1]
}

df = merge(mgs.med4, taxo, by.x=c("X"),by.y=c("mgs"))
df = melt(df)
head(df)
ggplot(data=df,aes(x=variable,y=value,fill=phylum)) + geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "top") + labs(x = "variable", y = "relative abundance")

mgs.med2 <- as.data.frame.matrix(mgs.med1)
rownames(mgs.med2) = mgs.med2[,1]
mgs.med2 = mgs.med2[,2:27]

mgs.med3 <- sweep(mgs.med2,2,colSums(mgs.med2),`/`)
mgs.med4 <- as.data.frame(mgs.med3)
mgs.med4$X <- rownames(mgs.med4)



## regulate table
colnames(merged) = colnames(merged)[2:27]
merged = merged[1:26]
head(merged)

index <- match(metadata$SampleID,colnames(merged))
int_exp <- as.matrix(unlist(merged[3,index]))
metadata$int_exp=int_exp

metadata$labels = paste(metadata$SampleID, metadata$Experiment, metadata$Feedin, metadata$Tissue, sep='_')
index1 <- match(colnames(merged),metadata$SampleID)
index_labels <- as.matrix(unlist(metadata[index1,8]))
merged1 <- merged
colnames(merged1) <- index_labels

rownames(merged_final) =  merged_final[,1]
merged_final = merged_final[2:27]



## Mapped read count
counting[,2] <- metadata[,2]
counting[,2] <- metadata[,3]
counting[,2] <- metadata[,4]

ggplot(data=counting, aes(x=sample, y=X.mapped_read_count, fill = condition)) + # define which columns you will display in plot
  geom_bar(stat="identity") + # plotting bar plot
  labs(x = "Sample", y = "mapped_read_count (%)",
       title = "Barplot of %mapped_read_count", fill = "Condition") + # define names for each axis and title
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) + # define the style of plot title
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



## Normalization
data(gene_lengths) 
str(gene_lengths)
norm.data <- normFreqRPKM(dat=merged_final, cat=gene_lengths)



## Downsizing
data.genenb <- downsizeGC(data = merged_final, level= c(1000000), repetitions=1, silent=TRUE)



## Gene richness
metadata$Gene_richness <- as.matrix(unlist(data.genenb[1,index]))

ggboxplot(metadata, x = "Experiment", y = "Gene_richness", 
          color = "Experiment",palette =c("#00AFBB", "#E7B800"), add = "jitter") +
  stat_compare_means(method = "t.test", label.x = 1.3, label= "p")

ggboxplot(metadata, x = "Feeding", y = "Gene_richness", 
          color = "Feeding",palette =c("#00AFBB", "#E7B800"), add = "jitter") +
  stat_compare_means(method = "t.test", label.x = 1.3, label= "p")

ggboxplot(metadata, x = "Tissue", y = "Gene_richness", 
          color = "Tissue",palette =c("#00AFBB", "#E7B800"), add = "jitter") +
  stat_compare_means(method = "t.test", label.x = 1.3, label= "p")



## Species richness
Species.number <- c()
for (i in 1:26){
  n.species <- 0
  for (j in 1:168){
    if(mgs.med.vect[j,i]>0){
      n.species=n.species+1
    } else{}
  }
  Species.number <- c(Species.number,n.species)
}
Species.richness <- as.matrix(Species.number)
rownames(Species.richness) = colnames(mgs.med.vect)
colnames(Species.richness) = "Species_richness"
metadata$Species_richness <- as.matrix(unlist(Species.richness[index,1]))

ggboxplot(metadata, x = "Experiment", y = "Species_richness", 
          color = "Experiment",palette =c("#00AFBB", "#E7B800"), add = "jitter") +
  stat_compare_means(method = "t.test", label.x = 1.3, label= "p")

ggboxplot(metadata, x = "Feeding", y = "Species_richness", 
          color = "Feeding",palette =c("#00AFBB", "#E7B800"), add = "jitter") +
  stat_compare_means(method = "t.test", label.x = 1.3, label= "p")

ggboxplot(metadata, x = "Tissue", y = "Species_richness", 
          color = "Tissue",palette =c("#00AFBB", "#E7B800"), add = "jitter") +
  stat_compare_means(method = "t.test", label.x = 1.3, label= "p")



## cluster analysis
hc = hclust(dist(t(merged1)))
plot(hc,hang=1,cex=1,labels = NULL)



## pheatmap
pheatmap(merged.hc.data$mat.rho)
pheatmap(merged.hc.data$mat.rho, cellwidth = 10, cellheight = 10)
pheatmap(merged.hc.data$mat.rho, color = colorRampPalette(c("white", "firebrick3"))(50))

annot = metadata
row.names(annot) = annot$SampleID
annot$SampleID = NULL
annot$labels = NULL

pheatmap(merged.hc.data$mat.rho, cellwidth = 10, cellheight = 10,
         color = colorRampPalette(c("white", "darkred"))(50),
         annotation_col = annot)



## MSPs through Wilcoxon-ranked singned test
mgs.list <- as.matrix(rownames(mgs.med.vect))
mgs_list <- as.data.frame(mgs.list)
index_mgs <- match(mgs_list$mgs,mgs_species$mgs)
int_species <- as.matrix(unlist(mgs_species[index_mgs,2]))
mgs_list$species = int_species

taxo_mgs <- as.matrix(unlist(taxo[index_mgs,6]))
mgs_list$labels = taxo_mgs

### DEN - nonDEN
p.Experiment <- c()
n.Experiment <- 0
for (i in 1:168){
  metadata$mgs <- as.matrix(unlist(mgs.med.vect[i,index]))
  res.Experiment <- wilcox.test(mgs ~ Experiment, data = metadata, paired = FALSE)
  p.Experiment <- c(p.Experiment,res.Experiment$p.value)
  if(res.Experiment$p.value<0.01){
    print(mgs_list[i,1]);n.Experiment=n.Experiment+1;
    filename.E=paste("E",i,".jpeg",sep="")
    jpeg(file=filename.E)
    print(ggboxplot(metadata, x = "Experiment", y = "mgs", 
                    color = "Experiment",palette =c("#00AFBB", "#E7B800"), add = "jitter") +
            stat_compare_means(method = "wilcox.test", label.x = 1.3, label= "p") + 
          labs(x = "Experiment", y = paste(mgs_list[i,2])))
    dev.off()
  } else{}
}
print(n.Experiment)
p.E <- as.matrix(p.Experiment)
rownames(p.E) <- rownames(mgs.med.vect)

for (i in 1:168){
  metadata$mgs <- as.matrix(unlist(mgs.med.vect[i,index]))
  res.Experiment <- wilcox.test(mgs ~ Experiment, data = metadata, paired = FALSE)
  if(res.Experiment$p.value<0.01){
    print(paste(rownames(mgs.med.vect)[i],mgs_list[i,2],mgs_list[i,3]));
  } else{}
}

### Fasted -Fed
p.Feeding <- c()
n.Feeding <- 0
for (i in 1:168){
  metadata$mgs <- as.matrix(unlist(mgs.med.vect[i,index]))
  res.Feeding <- wilcox.test(mgs ~ Feeding, data = metadata, paired = FALSE)
  p.Feeding <- c(p.Feeding,res.Feeding$p.value)
  if(res.Feeding$p.value<0.01){
    print(mgs_list[i,1]);n.Feeding=n.Feeding+1;
    filename.F=paste("F",i,".jpeg",sep="")
    jpeg(file=filename.F)
    print(ggboxplot(metadata, x = "Feeding", y = "mgs", 
                    color = "Feeding",palette =c("#00AFBB", "#E7B800"), add = "jitter") +
            stat_compare_means(method = "wilcox.test", label.x = 1.3, label= "p") + 
            labs(x = "Feeding", y = paste(mgs_list[i,2])))
    dev.off()
  } else{}
}
print(n.Feeding)
p.F <- as.matrix(p.Feeding)
rownames(p.F) <- rownames(mgs.med.vect)

for (i in 1:168){
  metadata$mgs <- as.matrix(unlist(mgs.med.vect[i,index]))
  res.Feeding <- wilcox.test(mgs ~ Feeding, data = metadata, paired = FALSE)
  if(res.Feeding$p.value<0.01){
    print(paste(rownames(mgs.med.vect)[i],mgs_list[i,2],mgs_list[i,3]));
  } else{}
}

### Faeces - Cecum
p.Tissue <- c()
n.Tissue <- 0
for (i in 1:168){
  metadata$mgs <- as.matrix(unlist(mgs.med.vect[i,index]))
  res.Tissue <- wilcox.test(mgs ~ Tissue, data = metadata, paired = TRUE)
  p.Tissue <- c(p.Tissue,res.Tissue$p.value)
  if(res.Tissue$p.value<0.01){
    print(mgs_list[i,1]);n.Tissue=n.Tissue+1;
    filename.T=paste("T",i,".jpeg",sep="")
    jpeg(file=filename.T)
    print(ggboxplot(metadata, x = "Tissue", y = "mgs", 
                    color = "Tissue",palette =c("#00AFBB", "#E7B800"), add = "jitter") +
            stat_compare_means(paired = TRUE, method = "wilcox.test", label.x = 1.3, label= "p") + 
            labs(x = "Tissue", y = paste(mgs_list[i,2])))
    dev.off()
  } else{}
}
print(n.Tissue)
p.T <- as.matrix(p.Tissue)
rownames(p.T) <- rownames(mgs.med.vect)

for (i in 1:168){
  metadata$mgs <- as.matrix(unlist(mgs.med.vect[i,index]))
  res.Tissue <- wilcox.test(mgs ~ Tissue, data = metadata, paired = TRUE)
  if(res.Tissue$p.value<0.01){
    print(paste(rownames(mgs.med.vect)[i],mgs_list[i,2],mgs_list[i,3]));
  } else{}
}




save.image(file = "~/Desktop/Project/data.Rdata")
