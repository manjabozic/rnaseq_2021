######################################################################################################
##################################### edgeR WORKFLOW ################################################
######################################################################################################





################################# DATA LOADING AND PREPARATION #######################################

#Set the working directory and load the necessary packages

setwd("E:/Manja/edgeR")
library(tidyverse)
library(edgeR)

#Since we are using htseq, we need to merge all the counts files into one first... To do that, we will use an R script through the LINUX 
#command line - htseq-combine_all.R

#Rscript htseq-combine_all.R <arg1> <arg2>
  
#arg1 is the path to where the files are and arg2 is the name of the final output
  
#Rscript htseq-combine_all.R ~/rnaseq_2021/04_read_counts/htseq/ counts
  
######################################################################################################

#Then we can start preparing the files for the DESeq2 analysis. First, load the merged count data:

counts <- read.delim("E:/Manja/deseq2/counts.csv", header=TRUE, row.names="ENSEMBL_GeneID")

#This is for a a data set with fake replicates; also an example of how to properly format the input file

counts_fake<- read.delim("E:/Manja/deseq2/counts_fake.csv", header=TRUE)
rownames(counts_fake) <- counts_fake[,1]
counts_fake[,1] <- NULL

#######################################################################################################

#edgeR stores data in a simple list-based data object called a DGEList.The function readDGE makes a
#DGEList object directly. If the table of counts is already available as a matrix or a data.frame, then 
#a DGEList object can be made by:

group <- c(1, 2) #Adds a grouping factor
group_fake = c (1, 1, 2, 2)

dgelist <- DGEList(counts=counts, group=group) # Creates a DGEList object
dgelist_fake <- DGEList(counts=counts_fake, group=group_fake)

#######################################################################################################
#Filtering genes with very low counts

keep <- filterByExpr(dgelist) #keeps rows that have worthwhile counts in a minumum number of samples
dgelist <- dgelist[keep, , keep.lib.sizes=FALSE] #Recalculates the library sizes

keep_fake <- filterByExpr(dgelist_fake)
dgelist_fake <- dgelist_fake[keep, , keep.lib.sizes=FALSE]

########################################################################################################
#Normalizing library sizes 

dgelist <- normLibSizes(dgelist)
dgelist_fake <- normLibSizes(dgelist_fake)


##################################### DE ESTIMATION - qCML ##############################################
###################################(SINGLE FACTOR EXPERIMENTS)###########################################


# edgeR uses the quantile-adjusted conditional maximum likelihood (qCML) method for experiments with 
# single factor. qCML is the most reliable in terms of bias on a wide range of conditions and specifically 
# performs best in the situation of many small samples with a common dispersion, the model which is 
# applicable to Next-Gen sequencing data.

dgelist <- estimateDisp(dgelist)
dgelist_fake <- estimateDisp(dgelist_fake)

et <- exactTest(dgelist)
et_fake <- exactTest(dgelist_fake)

##################################### DE ESTIMATION - GLM ##############################################
#################################### (COMPLEX EXPERIMENTS) #############################################


# Generalized linear models (GLMs) are an extension of classical linear models to nonnormally distributed 
# response data.For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
# likelihood (CR) method in estimating dispersions.It takes care of multiple factors by fitting generalized 
# linear models (GLM) with a design matrix.

group <- factor(c(1, 2)) #Adds a grouping factor
group_fake = factor(c (1, 1, 2, 2))

design <- model.matrix(~group)
design_fake <- model.matrix(~group_fake)

#it is recommended to use estimateDisp(y, design)

dgelist <- estimateDisp(dgelist, design) #common dispersion, trended dispersions, tagwise dispersions
dgelist <- estimateGLMCommonDisp(dgelist, design) #common dispersion
dgelist <- estimateGLMTrendedDisp(dgelist, design) #trended dispersions
dgelist <- estimateGLMTagwiseDisp(dgelist, design) #tagwise dispersions

dgelist_fake <- estimateDisp(dgelist_fake, design_fake) 

#We can proceed with fitting negative binomial generalized linear models and testing procedures for 
#determining differential expression using either quasi-likelihood (QL) F-test.

fit <- glmQLFit(dgelist, design) #fitting the glm model
fit_fake <- glmQLFit(dgelist_fake, design_fake)

qlf.2vs1 <- glmQLFTest(fit, coef=2) 
qlf.2vs1_fake <- glmQLFTest(fit_fake, coef=2) #comparing baseline of group 01 (C) and group (T)
topTags(qlf.2vs1_fake) # to view the best DEGs

# If you have an experiment with more factors (groups), this has to be taken into consideration 
# when writing the function. The default is to compare everything to the baseline of group 1, for example:

qlf.2vs1 <- glmQLFTest(fit, coef=2) # this compares the second to the first group
qlf.3vs1 <- glmQLFTest(fit, coef=3) # this compares the third to the first group

#If you want to compare different groups with each other, not with the first one:

qlf.3vs2 <- glmQLFTest(fit, contrast=c(0,-1,1)) # this compares the third to second group

# To find genes different between any of the three groups:

qlf <- glmQLFTest(fit, coef=2:3)
topTags(qlf)


#########################################################################################################
########################################## NO REPLICATES ################################################
#########################################################################################################

#You need to estimate an appropriate dispersion value
bcv <- 0.2
counts <- matrix( rnbinom(40,size=1/bcv^2,mu=10), 20,2) #I dont this step is necessary
de <- DGEList(counts=counts, group=1:2)
et <- exactTest(de, dispersion=bcv^2)

##########################################################################################################
##########################################################################################################

#Differential expression above a fold-change threshold
qlf <- glmQLFTest(fit, coef=2)
tr <- glmTreat(fit, coef=2, lfc=1)
topTags(tr)

#########################################################################################################
#Gene onthology and KEGG enrichment analysis

#you need to have certain Bioconductor packages installed: GO.db, AnnotationForge

#You also need to create a Zea mays org.package since it doesnt't exist in GO.db:
library(AnnotationForge)
makeOrgPackageFromNCBI(version = "1",
                       author = "Manja Bozic",
                       maintainer = "mb62095@uga.edu",
                       outputDir = ".",
                       tax_id = "4577",
                       genus = "Zea",
                       species = "mays")

qlf <- glmQLFTest(fit, coef=2)
go <- goana(fit_fake, species="Zm")
topGO(go, sort="up")
keg <- kegga(qlf, species="Zm")
topKEGG(keg, sort="up")
