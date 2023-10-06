######################################################################################################
##################################### DESEQ2 WORKFLOW ################################################
######################################################################################################





################################# DATA LOADING AND PREPARATION #######################################

#Set the working directory and load the necessary packages

setwd("E:/Manja/deseq2")
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

#Watch out for the file you are loading - add the absolute path to it, and see how to rows are named in it
#(it doesn't have to be "ENSEMBL_GeneID") - check this and change the parameters if necessary.


#After that, we can create the coldata data set. It is necessary to create a CSV file that holds all 
#the sample information before that.This files needs to have the sample names in different rows, and
#information about the treatment conditions, type of sequencing, etc - in separate columns.

coldata <- read.csv("E:/Manja/deseq2/coldata.csv", row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

#Again, be careful of the path to the CSV file and the namings of its columns.


#It is absolutely critical that the columns of the count matrix and the rows of the column data are in 
#the same order. Additionally, naming of the columns in the count matrix and rows in coldata must be
# consistent. You can check if this is true using these commands:

all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

#Both have to give the output "TRUE"

#Then just, load the data into DESeq2

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
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
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


################################# SETTING THE FACTOR LEVELS ############################################

# By default, R will choose a reference level for factors based on alphabetical order, so the
# comparisons will be based on the alphabetical order of the levels.
#This means, you need to tell DESeq2 which level you want to compare against (e.g. which level represents 
#the control group). 
#You can explicitly set the factors levels.

dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
dds$condition <- relevel(dds$condition, ref = "untreated")

#In order to see the change of reference levels reflected in the results names, you need to either run DESeq 
#or nbinomWaldTest/nbinomLRT after the re-leveling operation.

#Second option is to explicitly tell results (a DESeq2 function) which comparison to make using 
#the contrast argument in a later step.

############################################ DE ANALYSIS ################################################

dds_DESEQ2 <- DESeq(dds)
res <- results(dds)

###################################### ADDITIONAL ANALYSIS PARAMETERS ###################################

#You can sort the results by p-value:

resOrdered <- res[order(res$pvalue),]

summary(res) #The results summary

sum(res$padj < 0.1, na.rm=TRUE) #Number of DE genes with adjusted p-value less than 0.1

#Note that the results function automatically performs independent filtering based on the mean of 
#normalized counts for each gene, optimizing the number of genes which will have an adjusted p value 
#below a given FDR cutoff, alpha. 

#By default the argument alpha is set to 0.1, but you can change it:

res05 <- results(dds, alpha=0.05)

#If the adjusted p value cutoff will be a value other than 0.1, alpha should be set to that value


########################################## RESULTS ######################################################

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









