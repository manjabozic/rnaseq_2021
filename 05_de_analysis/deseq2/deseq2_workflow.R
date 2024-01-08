######################################################################################################
##################################### DESEQ2 WORKFLOW ################################################
######################################################################################################





################################# DATA LOADING AND PREPARATION #######################################

#Set the working directory and load the necessary packages

setwd("E:/Manja/deseq2/fake_reps")
library(tidyverse)
library(DESeq2)

#DESeq2 offers different options for loading data depending on its origin. Since we are using htseq, we need 
#to merge all the counts files into one first... To do that, we will use an R script through the LINUX 
#command line - htseq-combine_all.R

Rscript htseq-combine_all.R <arg1> <arg2>
  
#arg1 is the path to where the files are and arg2 is the name of the final output
  
#Rscript htseq-combine_all.R ~/rnaseq_2021/04_read_counts/htseq/ counts

######################################################################################################

#Then we can start preparing the files for the DESeq2 analysis. First, load the merged count data:

counts <- read.delim("E:/Manja/deseq2/counts.csv", header=TRUE, row.names="ENSEMBL_GeneID")
cts=as.matrix(counts)

total <- read.csv("E:/Manja/deseq2/fake_reps/counts_total.csv", header=TRUE, row.names="ENSEMBL_GeneID")
counts_total=as.matrix(total)


counts_T_6 <- as.matrix(read.csv("E:/Manja/deseq2/fake_reps/counts_II_6.csv", header=TRUE, row.names="ENSEMBL_GeneID"))
counts_S_6 <- as.matrix(read.csv("E:/Manja/deseq2/fake_reps/counts_III_6.csv", header=TRUE, row.names="ENSEMBL_GeneID"))

counts_T_24 <- as.matrix(read.csv("E:/Manja/deseq2/fake_reps/counts_II_24.csv", header=TRUE, row.names="ENSEMBL_GeneID"))
counts_S_24 <- as.matrix(read.csv("E:/Manja/deseq2/fake_reps/counts_III_24.csv", header=TRUE, row.names="ENSEMBL_GeneID"))

#Watch out for the file you are loading - add the absolute path to it, and see how to rows are named in it
#(it doesn't have to be "ENSEMBL_GeneID") - check this and change the parameters if necessary.


#After that, we can create the coldata data set. It is necessary to create a CSV file that holds all 
#the sample information before that.This files needs to have the sample names in different rows, and
#information about the treatment conditions, type of sequencing, etc - in separate columns.

coldata <- read.csv("E:/Manja/deseq2/fake_reps/coldata.csv", row.names=1)
coldata <- coldata[,c("condition","time_point", "genotype")]
coldata$condition <- factor(coldata$condition) #factor() encodes a vector as a factor
coldata$time_point <- factor(coldata$time_point)
coldata$genotype <- factor(coldata$genotype)

cd_T_6 <- read.delim("E:/Manja/deseq2/fake_reps/coldata_II_6.txt", sep = "\t", row.names=1)
cd_T_6$condition <- factor(cd_T_6$condition)

cd_T_24 <- read.delim("E:/Manja/deseq2/fake_reps/coldata_II_24.txt", sep = "\t", row.names=1)
cd_T_24$condition <- factor(cd_T_24$condition)

cd_S_6 <- read.csv("E:/Manja/deseq2/fake_reps/coldata_III_6.txt", sep = "\t", row.names=1)
cd_S_6$condition <- factor(cd_S_6$condition)

cd_S_24 <- read.csv("E:/Manja/deseq2/fake_reps/coldata_III_24.txt", sep = "\t", row.names=1)
cd_S_24$condition <- factor(cd_S_24$condition)

#Again, be careful of the path to the CSV file and the namings of its columns.


#It is absolutely critical that the columns of the count matrix and the rows of the column data are in 
#the same order. Additionally, naming of the columns in the count matrix and rows in coldata must be
# consistent. You can check if this is true using these commands:

all(rownames(coldata) %in% colnames(counts))
all(rownames(coldata) == colnames(counts))

all(rownames(cd_S_24) %in% colnames(counts_S_24))
all(rownames(cd_S_24) == colnames(counts_S_24))

all(rownames(cd_T_24) %in% colnames(counts_T_24))
all(rownames(cd_T_24) == colnames(counts_T_24))

all(rownames(cd_S_6) %in% colnames(counts_S_6))
all(rownames(cd_S_6) == colnames(counts_S_6))

all(rownames(cd_T_6) %in% colnames(counts_T_6))
all(rownames(cd_T_6) == colnames(counts_T_6))

#Both have to give the output "TRUE"

#Then just, load the data into DESeq2

dds <- DESeqDataSetFromMatrix(countData = counts_total,
                              colData = coldata,
                              design = ~ condition)

dds_T6 <- DESeqDataSetFromMatrix(countData = counts_T_6,
                              colData = cd_T_6,
                              design = ~ condition)

dds_T24 <- DESeqDataSetFromMatrix(countData = counts_T_24,
                                   colData = cd_T_24,
                                   design = ~ condition)

dds_S6 <- DESeqDataSetFromMatrix(countData = counts_S_6,
                                 colData = cd_S_6,
                                 design = ~ condition)

dds_S24 <- DESeqDataSetFromMatrix(countData = counts_S_24,
                                  colData = cd_S_24,
                                  design = ~ condition)

#########################################################################################################

#There is also an option to do it without previously merging the files before:

directory <- "/path/to/your/files/" #Set a directory where the individual htseq files are.
sampleFiles <- grep("cut_",list.files(directory),value=TRUE) #Makes a list of sample names

#But be mindful of how you grep it - use some pattern found in all of the files, that you can also use 
#in the next step.

sampleCondition <- sub("(.*treated).*","\\1",sampleFiles) #Makes a list of conditions for each of your files.
#If you cannot do it with sub, write them out on your own - be careful of the order you are writing them in.
#sampleCondition = c("untreated", "treated")

sampleNames <- gsub('.tsv', '', sampleFiles)

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

######################################### PREFILTERING #################################################

#Remove rows with a minimum amount of reads (in this case 10). Smallest group size iz the number of different
# samples (not including technical replicates). Since we have no replicates, it will always be one.

smallestGroupSize <- 1
keep <- rowSums(counts(dds_all) >= 10) >= smallestGroupSize
dds_all <- dds_all[keep,]

keep <- rowSums(counts(dds_T6) >= 10) >= smallestGroupSize
dds_T6 <- dds_T6[keep,]
keep <- rowSums(counts(dds_T24) >= 10) >= smallestGroupSize
dds_T24 <- dds_T24[keep,]

keep <- rowSums(counts(dds_S6) >= 10) >= smallestGroupSize
dds_S6 <- dds_S6[keep,]
keep <- rowSums(counts(dds_S24) >= 10) >= smallestGroupSize
dds_S24 <- dds_S24[keep,]

################################# SETTING THE FACTOR LEVELS ############################################

# By default, R will choose a reference level for factors based on alphabetical order, so the
# comparisons will be based on the alphabetical order of the levels.
#This means, you need to tell DESeq2 which level you want to compare against (e.g. which level represents 
#the control group). 
#You can explicitly set the factors levels.

dds_all$condition <- factor(dds_all$condition, levels = c("untreated","treated"))
dds_all$genotype <- factor(dds_all$genotype, levels = c("II","III"))
dds_all$time_point <- factor(dds_all$time_point, levels = c("6","24"))
dds_all$condition <- relevel(dds_all$condition, ref = "untreated")


dds_T6$condition <- factor(dds_T6$condition, levels = c("untreated","treated"))
dds_T6$condition <- relevel(dds_T6$condition, ref = "untreated")
dds_T24$condition <- factor(dds_T24$condition, levels = c("untreated","treated"))
dds_T24$condition <- relevel(dds_T24$condition, ref = "untreated")

dds_S6$condition <- factor(dds_S6$condition, levels = c("untreated","treated"))
dds_S6$condition <- relevel(dds_S6$condition, ref = "untreated")
dds_S24$condition <- factor(dds_S24$condition, levels = c("untreated","treated"))
dds_S24$condition <- relevel(dds_S24$condition, ref = "untreated")

#In order to see the change of reference levels reflected in the results names, you need to either run DESeq 
#or nbinomWaldTest/nbinomLRT after the re-leveling operation.

#Second option is to explicitly tell results (a DESeq2 function) which comparison to make using 
#the contrast argument in a later step.

############################################ DE ANALYSIS ################################################

dds_DESEQ2 <- DESeq(dds)
res <- results(dds)

Tol_6h_DE <- DESeq(dds_T6)
Tol_6h_res <- results(Tol_6h_DE, alpha = 0.05)

Tol_24h_DE <- DESeq(dds_T24)
Tol_24h_res <- results(Tol_24h_DE, alpha = 0.05)

Susc_6h_DE <- DESeq(dds_S6)
Susc_6h_res <- results(Susc_6h_DE, alpha = 0.05)

Susc_24h_DE <- DESeq(dds_S24)
Susc_24h_res <- results(Susc_24h_DE, alpha = 0.05)

###################################### ADDITIONAL ANALYSIS PARAMETERS ###################################
summary(Tol_6h_res)#The results summary
summary(Tol_24h_res)

summary(Susc_6h_res)
summary(Susc_24h_res)

#You can sort the results by p-value:

resOrdered <- res[order(res$pvalue),]

Tol_6h_res_sort <- Tol_6h_res[order(Tol_6h_res$padj),]
Tol_24h_res_sort <- Tol_24h_res[order(Tol_24h_res$padj),]

Susc_6h_res_sort <- Susc_6h_res[order(Susc_6h_res$padj),]
Susc_24h_res_sort <- Susc_24h_res[order(Susc_24h_res$padj),]

sum(Tol_6h_res$padj < 0.05, na.rm=TRUE) #Number of DE genes with adjusted p-value less than 0.1
sum(Tol_24h_res$padj < 0.05, na.rm=TRUE)

sum(Susc_6h_res$padj < 0.05, na.rm=TRUE)
sum(Susc_24h_res$padj < 0.05, na.rm=TRUE)
#Note that the results function automatically performs independent filtering based on the mean of 
#normalized counts for each gene, optimizing the number of genes which will have an adjusted p value 
#below a given FDR cutoff, alpha. 

#By default the argument alpha is set to 0.1, but you can change it:

res05 <- results(dds, alpha=0.05)

#If the adjusted p value cutoff will be a value other than 0.1, alpha should be set to that value


########################################## RESULTS ######################################################

#Log fold change shrinkage

resLFC_Tol_6h <- lfcShrink(Tol_6h_DE, coef="condition_treated_vs_untreated", type="apeglm")
resLFC_Tol_24h <- lfcShrink(Tol_24h_DE, coef="condition_treated_vs_untreated", type="apeglm")

resLFC_Susc_6h <- lfcShrink(Susc_6h_DE, coef="condition_treated_vs_untreated", type="apeglm")
resLFC_Susc_24h <- lfcShrink(Susc_24h_DE, coef="condition_treated_vs_untreated", type="apeglm")

#MA-plot (The log2 fold changes attributable to a given variable over the mean of normalized counts for 
#all the samples in the DESeqDataSet.)

plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2)) #shrunken log2 fold changes
idx <- identify(res$baseMean, res$log2FoldChange) #detect the row number of individual genes - gene identifiers
rownames(res)[idx]

#You can also change the shrinkage estimators

resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")


#Plot counts (examine the counts of reads for a single gene across the groups; The counts are grouped 
#by the variables in intgroup)

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

#For customized plotting, an argument returnData specifies that the function should only return a 
#data.frame for plotting with ggplot.

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#ReportingTools - HTML report of the results with plots and sortable/filterable columns can be generated 
#using the ReportingTools package on a DESeqDataSet that has been processed by the DESeq function

#regionReport - HTML and PDF summary of the results with plots can also be generated using the regionReport package, 
#after finishing the DESeq2Report.

#Glimma - Interactive visualization of DESeq2 output, including MA-plots (also called MD-plots) 
#can be generated using the Glimma package.

#pcaExplorer - Interactive visualization of DESeq2 output, including PCA plots, boxplots of counts and 
#other useful summaries can be generated using the pcaExplorer package.

#iSEE - Provides functions for creating an interactive Shiny-based graphical user interface for 
#exploring data stored in SummarizedExperiment objects, including row- and column-level metadata. 

#DEvis DEvis is a powerful, integrated solution for the analysis of differential expression data. 

#Exporting results to CSV files

write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")

############################# DATA TRANSFORMATION AND VISUALIZATION ####################################
# variance stabilizing transformations (VST)
vsd <- vst(dds, blind=FALSE)

#regularized logarithm TRANSFORMATIONS (rlog)
rld <- rlog(dds, blind=FALSE)

#Heatmap of the count matrix
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#Principal component plot of the samples
plotPCA(vsd, intgroup=c("condition", "type"))
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()









