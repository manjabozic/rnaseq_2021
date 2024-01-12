#########################################################################################################################################
############################### Weighted Gene Co-expression Network Analysis (WGCNA)  #########################################################
##########################################################################################################################################

# Script to perform weighted Gene Co-expression Network Analysis (WGCNA) 

setwd("E:/Manja/wgcna")

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

allowWGCNAThreads()          # allow multi-threading (optional)

################################################################################################################
# 1. Fetch data and metadata
###############################################################################################################

counts <- read.delim('E:/Manja/wgcna/normalized_counts.txt', header = T, row.names = 1)
counts <- read.delim('E:/Manja/wgcna/htseq_name_counts.txt', header = T, row.names = 1)

exp.data <- read.delim('E:/Manja/wgcna/phenodata.txt', header = T, row.names = 1)
exp.data <- exp.data[,c("Description","Time_point", "Genotype")]
exp.data$Description <- factor(exp.data$Description) #factor() encodes a vector as a factor
exp.data$Time_point <- factor(exp.data$Time_point)
exp.data$Genotype <- factor(exp.data$Genotype)
  

################################################################################################################
# 2.  QC - outlier detection
###############################################################################################################
# There multiple ways to check: goodSamplesGenes, hierarchical clustering, PCA...

# goodSamplesGenes
qc <- goodSamplesGenes(t(counts)) #rows should be samples, and columns should be genes 
summary(qc)
qc$allOK
table(qc$goodGenes)
table(qc$goodSamples)

counts <- counts[qc$goodGenes == TRUE,] #if there are outliers this removes them

# Hierarchical clustering 
htree <- hclust(dist(t(counts)), method = "average")
plot(htree)

# PCA
pca <- prcomp(t(counts))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

samples.to.be.excluded <- c('C_6_Susc', 'C_6_Tol') # exclude outlier samples
counts.subset <- counts[,!(colnames(counts) %in% samples.to.be.excluded)]

### NOTE: If there are batch effects observed, correct for them before moving ahead

################################################################################################################
# 3. Normalization
###############################################################################################################

# You can do this in DESeq2 or here. If it is being done here the data loaded in the first step must be raw counts, that
# were not transformed in any way (for example, htseq_name_counts.txt)

# Create a dds file
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = exp.data,
                              design = ~ Description) # not specifying model


## Remove all genes with counts < 15 in more than 75% of samples (8*0.75=6) - suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 8,]
nrow(dds75) # 13284 genes

# Perform variance stabilization
dds_norm <- vst(dds75)

# Get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

counts <- t(counts)
#counts.subset <- t(counts.subset)
#counts_scaled <- scale(counts)
################################################################################################################
# 4. Network Construction
###############################################################################################################

# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2)) # creates a vecotr of possible soft tresholds

# Call the network topology analysis function, to choose the treshold
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft <- pickSoftThreshold(counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

demo_soft <- pickSoftThreshold(demo_norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

# Metrics used to decde which power gives a scale free topology are found in the fitIndices section - R squared value (SFT.R.sq) 
# and mean connectivity (max.k). ideally it is the power that gives the biggest R square value and smallest mean connectivity.

sft.data <- sft$fitIndices
demo_soft.data <- demo_soft$fitIndices

# Visualization to pick power

a1 <- ggplot(sft.data , aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data , aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)

# Next step is to identify modules - start by creating an adjacency matrix and applying the soft threshold to create a 
# weighted correlation matrix. Then we calculate the proximity measures - the topological overlap measures.

# Convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)
counts[] <- sapply(counts, as.numeric)

soft_power20 <- 20 # soft powers chosen in the previous step
soft_power36 <- 36
temp_cor <- cor # assigns the correlation function to a temp variable, so the WGCNA uses the correlation function from their package
cor <- WGCNA::cor #this assigns the WGCNA correlation function to the correlation function


# memory estimate w.r.t blocksize
bwnet_30 <- blockwiseModules(norm.counts,
                          maxBlockSize = 8000, #depends on the RAM
                          TOMType = "signed",
                          power = 30,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

bwnet_20 <- blockwiseModules(counts,
                             maxBlockSize = 8000, #depends on the RAM
                             TOMType = "signed",
                             power = 20,
                             mergeCutHeight = 0.25,
                             numericLabels = FALSE,
                             randomSeed = 1234,
                             verbose = 3)

bwnet_demo <- blockwiseModules(demo_norm.counts,
                          maxBlockSize = 8000,
                          TOMType = "signed",
                          power = soft_power18,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------

module_eigengenes_og <- bwnet_20$MEs
module_eigengenes_newnorm <- bwnet_30$MEs
module_eigengenes_demo <- bwnet_demo$MEs

head(module_eigengenes) # Print out a preview

table(bwnet$colors) # get number of genes for each module
table(bwnet$unmergedColors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_20$dendrograms[[1]], cbind(bwnet_20$unmergedColors, bwnet_20$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

plotDendroAndColors(bwnet_30$dendrograms[[1]], cbind(bwnet_30$unmergedColors, bwnet_30$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


# grey module = all genes that doesn't fall into other modules were assigned to the grey module




# 6A. Relate modules to traits --------------------------------------------------
# module trait associations



# create traits file - binarize categorical variables
traits <- exp.data %>% 
  mutate(Control_vs_Treatment = ifelse(grepl('Treatment', Description), 1, 0)) %>% 
  select(4)

demo_traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl('COVID', disease_state), 1, 0)) %>% 
  select(8)


# binarize categorical variables

exp.data$Genotype<- factor(exp.data$Genotype, levels = c("Susceptible", "Tolerant"))

genotype.out <- binarizeCategoricalColumns(exp.data$Genotype,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1,
                                           nameForAll = "Susceptible",)


traits <- cbind(traits, genotype.out)


# Define numbers of genes and samples
nSamples_norm <- nrow(norm.counts)
nGenes_norm <- ncol(norm.counts) 

nSamples <- nrow(counts)
nGenes <- ncol(counts)

module.trait.corr_og <- cor(module_eigengenes_og, traits, use = 'p')
module.trait.corr.pvals_og <- corPvalueStudent(module.trait.corr_og, nSamples)

module.trait.corr_newnorm <- cor(module_eigengenes_newnorm, traits, use = 'p')
module.trait.corr.pvals_newnorm <- corPvalueStudent(module.trait.corr_newnorm, nSamples_norm)


# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes_og, traits, by = 'row.names') #combine these two data sets by row names
heatmap.data_norm<- merge(module_eigengenes_newnorm, traits, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

heatmap.data_norm <- heatmap.data_norm %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[60:61], #choose the column names that contain the trait data (condition and genotype)
             y = names(heatmap.data)[1:59], #choose the column names that contain all the eigen gene names
             col = c("blue1", "skyblue", "white", "pink", "red"))

CorLevelPlot(heatmap.data_norm,
             x = names(heatmap.data_norm)[54:55], #choose the column names that contain the trait data (condition and genotype)
             y = names(heatmap.data_norm)[1:53], #choose the column names that contain all the eigen gene names
             col = c("blue1", "skyblue", "white", "pink", "red"))

# first column lists gene clusters associated with the the treated state (Control_vs_Treatment) 
# second column lists gene clusters associated with the the tolerant genotype (data.Tolerant.vs.Susceptible)
# focus only on the statistically significant ones (*, **, ***)

# From the heatmap, you can see what modules are significantly associated with the traits of interest - they are labeled ***,
# **, and * according to the p-values. You can see which modules are significantly associated with the condition - Treatment 
# and which are significantly associated with the tolerant genotype

module.gene.mapping <- as.data.frame(bwnet_20$colors)
module.gene.mapping_norm <- as.data.frame(bwnet_30$colors)

write.table(module.gene.mapping, "E:/Manja/wgcna/module.gene.mapping_og.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(module.gene.mapping_norm, "E:/Manja/wgcna/module.gene.mapping_newnorm.txt", sep = "\t", row.names = TRUE, col.names = NA)

module01 <- module.gene.mapping %>% 
  filter(`bwnet_20$colors` == 'paleturquoise') %>% 
  rownames()

module02 <- module.gene.mapping %>% 
  filter(`bwnet_20$colors` == 'salmon') %>% 
  rownames()

module03 <- module.gene.mapping %>% 
  filter(`bwnet_20$colors` == 'lightyellow') %>% 
  rownames()

module04 <- module.gene.mapping %>% 
  filter(`bwnet_20$colors` == 'tan') %>% 
  rownames()

module05 <- module.gene.mapping %>% 
  filter(`bwnet_20$colors` == 'pink') %>% 
  rownames()

module06 <- module.gene.mapping %>% 
  filter(`bwnet_20$colors` == 'green') %>% 
  rownames()



# 6B. Intramodular analysis: Identifying driver genes ---------------



# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes_og, counts, use = 'p') # correlation method is Pearson ('p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure_newnorm <- cor(module_eigengenes_newnorm, norm.counts, use = 'p')
module.membership.measure.pvals_newnorm <- corPvalueStudent(module.membership.measure_newnorm, nSamples_norm)

write.table(module.membership.measure.pvals, "E:/Manja/wgcna/module.membership.measure.pvals_og.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(module.membership.measure.pvals_newnorm, "E:/Manja/wgcna/module.membership.measure.pvals_newnorm.txt", sep = "\t", row.names = TRUE, col.names = NA)

module.membership.measure.pvals[1:10,1:10]


# Calculate the gene significance and associated p-values

gene.signf.corr_cvst <- cor(counts, traits$Control_vs_Treatment, use = 'p') # correlation method is Pearson ('p')
gene.signf.corr.pvals_cvst <- corPvalueStudent(gene.signf.corr_cvst, nSamples)

gene.signf.corr_tolerant <- cor(counts, traits$data.Tolerant.vs.Susceptible, use = 'p')
gene.signf.corr.pvals_tolerant <- corPvalueStudent(gene.signf.corr_tolerant, nSamples)

gene.signf.corr_cvst_newnorm <- cor(norm.counts, traits$Control_vs_Treatment, use = 'p') # correlation method is Pearson ('p')
gene.signf.corr.pvals_cvst_newnorm <- corPvalueStudent(gene.signf.corr_cvst_newnorm, nSamples_norm)

gene.signf.corr_tolerant_newnorm <- cor(norm.counts, traits$data.Tolerant.vs.Susceptible, use = 'p')
gene.signf.corr.pvals_tolerant_newnorm <- corPvalueStudent(gene.signf.corr_tolerant_newnorm, nSamples_norm)

gene.signf.corr.pvals_cvst  %>% 
  as.data.frame() %>% #transform into a data frame
  arrange(V1) %>% #arrange in ascending order
  write.table("E:/Manja/wgcna/gene.signf.corr.pvals_cvst_og.txt", sep = "\t", row.names = TRUE, col.names = NA) %>% 
  head(25) # show the top 25

gene.signf.corr.pvals_tolerant  %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  write.table("E:/Manja/wgcna/gene.signf.corr.pvals_tolerant_og.txt", sep = "\t", row.names = TRUE, col.names = NA) %>% 
  head(25)

gene.signf.corr.pvals_cvst_newnorm  %>% 
  as.data.frame() %>% #transform into a data frame
  arrange(V1) %>% #arrange in ascending order
  write.table("E:/Manja/wgcna/gene.signf.corr.pvals_cvst_newnorm.txt", sep = "\t", row.names = TRUE, col.names = NA) %>% 
  head(25) # show the top 25

gene.signf.corr.pvals_tolerant_newnorm  %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  write.table("E:/Manja/wgcna/gene.signf.corr.pvals_tolerant_newnorm.txt", sep = "\t", row.names = TRUE, col.names = NA) %>% 
  head(25)


# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.