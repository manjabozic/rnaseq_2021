setwd("/lustre2/scratch/mb62095/r")



install.packages("FactoMineR")
install.packages("factoextra")
install.packages("cluster")
install.packages("stats")
install.packages("fpc")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("baySeq")

install.packages("dendextend")
install.packages("gplots")
install.packages("pheatmap")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

install.packages("hopkins")
install.packages("NbClust")
install.packages("clvalid")
install.packages("pvclust")
install.packages("dbscan")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterSeq")

#############################################################################################################

library("cluster") # load the neccessary packages
library("FactoMineR")
library("factoextra") #x
library("stats")
library("fpc")
library("baySeq") #x
library("dendextend")
library("gplots")
library("pheatmap")
library("ComplexHeatmap") #x
library("hopkins") #x
library("NbClust") #x
library("clValid")
library("pvclust")
library("dbscan")
library("clusterSeq")

#############################################################################################################

counts_mor <- as.matrix(read.delim("/lustre2/scratch/mb62095/r/normalized_counts.txt", header=TRUE, row.names=1))

counts <- read.delim("/lustre2/scratch/mb62095/r/cleandata.txt", header=TRUE, row.names=1)

#1. Clusterseq package normalization method

ls_all <- getLibsizes(data = counts)
counts[counts == 0] <- 1
norm_counts <- log2(t(t(counts / ls_all)) * mean(ls_all))

#2. Scaling the data 

counts <- scale(counts)

###############################################################################################################

set.seed(123)

pdf(file = '/lustre2/scratch/mb62095/r/scaled.pdf')
fviz_dist(dist(counts), show_labels = FALSE) + labs(title = "Scaled, readcount")
dev.off()

pdf(file = '/lustre2/scratch/mb62095/r/normalized.pdf')
fviz_dist(dist(norm_counts), show_labels = FALSE) + labs(title = "Normalized, log2, readcount")
dev.off()

pdf(file = '/lustre2/scratch/mb62095/r/mean_of_ratios.pdf')
fviz_dist(dist(counts_mor), show_labels = FALSE) + labs(title = "Mean of ratios normalization method, readcount")
dev.off()
