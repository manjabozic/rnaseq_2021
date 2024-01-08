######################################################################################################
##################################### edgeR WORKFLOW ################################################
######################################################################################################





################################# DATA LOADING AND PREPARATION #######################################

#Set the working directory and load the necessary packages

setwd("E:/Manja/edgeR")
library(tidyverse)
library(edgeR)
library(ggplot2)
library(EnhancedVolcano)

#Since we are using htseq, we need to merge all the counts files into one first... To do that, we will use an R script through the LINUX 
#command line - htseq-combine_all.R

#Rscript htseq-combine_all.R <arg1> <arg2>
  
#arg1 is the path to where the files are and arg2 is the name of the final output
  
#Rscript htseq-combine_all.R ~/rnaseq_2021/04_read_counts/htseq/ counts
  
######################################################################################################

#Then we can start preparing the files for the DESeq2 analysis. First, load the merged count data:

#for all counts
counts<- read.delim("E:/Manja/edgeR/data/htseq_name_counts.txt", header=TRUE, row.names="ENSEMBL_GeneID")

# control vs treatment, susceptible genotype
counts_Susc_6h <- read.delim("E:/Manja/edgeR/data/Susc_6h.txt", header=TRUE, row.names="ENSEMBL_GeneID")
counts_Susc_24h <- read.delim("E:/Manja/edgeR/data/Susc_24h.txt", header=TRUE, row.names="ENSEMBL_GeneID")

# control vs treatment, tolerant genotype
counts_Tol_6h <- read.delim("E:/Manja/edgeR/data/Tol_6h.txt", header=TRUE, row.names="ENSEMBL_GeneID")
counts_Tol_24h <- read.delim("E:/Manja/edgeR/data/Tol_24h.txt", header=TRUE, row.names="ENSEMBL_GeneID")

# tolerant vs susceptible genotype
counts_c_6h <- read.delim("E:/Manja/edgeR/data/control_6h.txt", header=TRUE, row.names="ENSEMBL_GeneID")
counts_c_24h <- read.delim("E:/Manja/edgeR/data/control_24h.txt", header=TRUE, row.names="ENSEMBL_GeneID")

counts_t_6h <- read.delim("E:/Manja/edgeR/data/treatment_6h.txt", header=TRUE, row.names="ENSEMBL_GeneID")
counts_t_24h <- read.delim("E:/Manja/edgeR/data/treatment_24h.txt", header=TRUE, row.names="ENSEMBL_GeneID")


#######################################################################################################

#edgeR stores data in a simple list-based data object called a DGEList.The function readDGE makes a
#DGEList object directly. If the table of counts is already available as a matrix or a data.frame, then 
#a DGEList object can be made by:

group <- c(1, 2) #Adds a grouping factor

# Creates a DGEList object
dgelist_Susc_6h <- DGEList(counts=counts_Susc_6h, group=group)
dgelist_Susc_24h <- DGEList(counts=counts_Susc_24h, group=group)
dgelist_Tol_6h <- DGEList(counts=counts_Tol_6h, group=group)
dgelist_Tol_24h <- DGEList(counts=counts_Tol_24h, group=group)

dgelist_c_6h <- DGEList(counts=counts_c_6h, group=group)
dgelist_c_24h <- DGEList(counts=counts_c_24h, group=group)
dgelist_t_6h <- DGEList(counts=counts_t_6h, group=group)
dgelist_t_24h <- DGEList(counts=counts_t_24h, group=group)
#######################################################################################################
#Filtering genes with very low counts

keep <- filterByExpr(dgelist_Susc_6h, min.count = 15, min.total.count = 30) #keeps rows that have worthwhile counts in a minumum number of samples
dgelist_Susc_6h <- dgelist_Susc_6h[keep, , keep.lib.sizes=FALSE] #Recalculates the library sizes

keep <- filterByExpr(dgelist_Susc_24h, min.count = 15, min.total.count = 30) 
dgelist_Susc_24h <- dgelist_Susc_24h[keep, , keep.lib.sizes=FALSE]

keep <- filterByExpr(dgelist_Tol_6h, min.count = 15, min.total.count = 30) 
dgelist_Tol_6h <- dgelist_Tol_6h[keep, , keep.lib.sizes=FALSE]

keep <- filterByExpr(dgelist_Tol_24h, min.count = 15, min.total.count = 30) 
dgelist_Tol_24h <- dgelist_Tol_24h[keep, , keep.lib.sizes=FALSE]

keep <- filterByExpr(dgelist_c_6h, min.count = 15, min.total.count = 30) 
dgelist_c_6h <- dgelist_c_6h[keep, , keep.lib.sizes=FALSE] 

keep <- filterByExpr(dgelist_c_24h, min.count = 15, min.total.count = 30) 
dgelist_c_24h <- dgelist_c_24h[keep, , keep.lib.sizes=FALSE]

keep <- filterByExpr(dgelist_t_6h, min.count = 15, min.total.count = 30) 
dgelist_t_6h <- dgelist_t_6h[keep, , keep.lib.sizes=FALSE]

keep <- filterByExpr(dgelist_t_24h, min.count = 15, min.total.count = 30) 
dgelist_t_24h <- dgelist_t_24h[keep, , keep.lib.sizes=FALSE]

########################################################################################################
#Normalizing library sizes 

dgelist_Susc_6h <- normLibSizes(dgelist_Susc_6h)
dgelist_Susc_24h <- normLibSizes(dgelist_Susc_24h)
dgelist_Tol_6h <- normLibSizes(dgelist_Tol_6h)
dgelist_Tol_6h <- normLibSizes(dgelist_Tol_6h)

dgelist_c_6h <- normLibSizes(dgelist_c_6h)
dgelist_c_24h <- normLibSizes(dgelist_c_24h)
dgelist_t_6h <- normLibSizes(dgelist_t_6h)
dgelist_t_6h <- normLibSizes(dgelist_t_6h)


##################################### DE ESTIMATION - qCML ##############################################
###################################(SINGLE FACTOR EXPERIMENTS)###########################################


# edgeR uses the quantile-adjusted conditional maximum likelihood (qCML) method for experiments with 
# single factor. qCML is the most reliable in terms of bias on a wide range of conditions and specifically 
# performs best in the situation of many small samples with a common dispersion, the model which is 
# applicable to Next-Gen sequencing data.

#Since there are no replicates dispersion values can be set manually individually.

qcml_Susc_6h <- estimateDisp(dgelist_Susc_6h)
qcml_Susc_6h$common.dispersion = 0.09
qcml_Susc_6h$tagwise.dispersion = 0.09
qcml_Susc_6h$trended.dispersion = 0.09

qcml_Susc_24h <- estimateDisp(dgelist_Susc_24h)
qcml_Susc_24h$common.dispersion = 0.09
qcml_Susc_24h$tagwise.dispersion = 0.09
qcml_Susc_24h$trended.dispersion = 0.09

qcml_Tol_6h <- estimateDisp(dgelist_Tol_6h)
qcml_Tol_6h$common.dispersion = 0.09
qcml_Tol_6h$tagwise.dispersion = 0.09
qcml_Tol_6h$trended.dispersion = 0.09

qcml_Tol_24h <- estimateDisp(dgelist_Tol_24h)
qcml_Tol_24h$common.dispersion = 0.09
qcml_Tol_24h$tagwise.dispersion = 0.09
qcml_Tol_24h$trended.dispersion = 0.09

qcml_c_6h <- estimateDisp(dgelist_c_6h)
qcml_c_6h$common.dispersion = 0.09
qcml_c_6h$tagwise.dispersion = 0.09
qcml_c_6h$trended.dispersion = 0.09

qcml_c_24h <- estimateDisp(dgelist_c_24h)
qcml_c_24h$common.dispersion = 0.09
qcml_c_24h$tagwise.dispersion = 0.09
qcml_c_24h$trended.dispersion = 0.09

qcml_t_6h <- estimateDisp(dgelist_t_6h)
qcml_t_6h$common.dispersion = 0.09
qcml_t_6h$tagwise.dispersion = 0.09
qcml_t_6h$trended.dispersion = 0.09

qcml_t_24h <- estimateDisp(dgelist_t_24h)
qcml_t_24h$common.dispersion = 0.09
qcml_t_24h$tagwise.dispersion = 0.09
qcml_t_24h$trended.dispersion = 0.09

#Or just by setting the BCV value

bcv = 0.4

#Depending on what you choose the dispersion value can be deleted from the command

qcml_Susc_6h <- exactTest(qcml_Susc_6h, dispersion=bcv^2)
qcml_Susc_24h <- exactTest(qcml_Susc_24h, dispersion=bcv^2)
qcml_Tol_6h <- exactTest(qcml_Tol_6h, dispersion=bcv^2)
qcml_Tol_24h <- exactTest(qcml_Tol_24h, dispersion=bcv^2)

qcml_Susc_6h <- exactTest(qcml_Susc_6h)
qcml_Susc_24h <- exactTest(qcml_Susc_24h)
qcml_Tol_6h <- exactTest(qcml_Tol_6h)
qcml_Tol_24h <- exactTest(qcml_Tol_24h)

qcml_c_6h <- exactTest(qcml_c_6h)
qcml_c_24h <- exactTest(qcml_c_24h)
qcml_t_6h <- exactTest(qcml_t_6h)
qcml_t_24h <- exactTest(qcml_t_24h)

#Extracting the results

DE_Susc_6h <- as.data.frame(topTags(qcml_Susc_6h, n = 25640, adjust.method = "BH", sort.by = "PValue", p.value = 0.1))
DE_Susc_24h <- as.data.frame(topTags(qcml_Susc_24h, n = 25198, adjust.method = "BH", sort.by = "PValue", p.value = 0.1))

DE_Tol_6h <- as.data.frame(topTags(qcml_Tol_6h, n = 25396, adjust.method = "BH", sort.by = "PValue", p.value = 0.1))
DE_Tol_24h <- as.data.frame(topTags(qcml_Tol_24h, n = 25436, adjust.method = "BH", sort.by = "PValue", p.value = 0.1))

write.table(DE_Susc_6h, "E:/Manja/edgeR/results/DE_Susc_6h.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(DE_Susc_24h, "E:/Manja/edgeR/results/DE_Susc_24h.txt", sep = "\t", row.names = TRUE, col.names = NA)

write.table(DE_Tol_6h, "E:/Manja/edgeR/results/DE_Tol_6h.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(DE_Tol_24h, "E:/Manja/edgeR/results/DE_Tol_24h.txt", sep = "\t", row.names = TRUE, col.names = NA)


DE_c_6h <- as.data.frame(topTags(qcml_c_6h, n = 26126, adjust.method = "BH", sort.by = "PValue", p.value = 0.05))
DE_c_24h <- as.data.frame(topTags(qcml_c_24h, n = 26345, adjust.method = "BH", sort.by = "PValue", p.value = 0.05))

DE_t_6h <- as.data.frame(topTags(qcml_t_6h, n = 25516, adjust.method = "BH", sort.by = "PValue", p.value = 0.05))
DE_t_24h <- as.data.frame(topTags(qcml_t_24h, n = 25091, adjust.method = "BH", sort.by = "PValue", p.value = 0.05))

write.table(DE_c_6h, "E:/Manja/edgeR/results/DE_control_6h.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(DE_c_24h, "E:/Manja/edgeR/results/DE_control_24h.txt", sep = "\t", row.names = TRUE, col.names = NA)

write.table(DE_t_6h, "E:/Manja/edgeR/results/DE_treatment_6h.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(DE_t_24h, "E:/Manja/edgeR/results/DE_treatment_24h.txt", sep = "\t", row.names = TRUE, col.names = NA)

##################################### DE ESTIMATION - GLM ##############################################
#################################### (COMPLEX EXPERIMENTS) #############################################


# Generalized linear models (GLMs) are an extension of classical linear models to nonnormally distributed 
# response data.For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
# likelihood (CR) method in estimating dispersions.It takes care of multiple factors by fitting generalized 
# linear models (GLM) with a design matrix.

counts_htseq <- read.delim("E:/Manja/edgeR/data/htseq_name_counts.txt", header=TRUE, row.names="ENSEMBL_GeneID")
group_individual <- c(1, 2, 3, 4, 5, 6, 7, 8)
dge <- DGEList(counts = counts_htseq, group = group_individual)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- normLibSizes(dge)

counts_fr <- read.delim("E:/Manja/edgeR/data/counts_total.txt", header=TRUE, row.names="ENSEMBL_GeneID")
group_fr <- c("s_c_6", "s_c_6", "s_t_6", "s_t_6","s_c_24", "s_c_24", "s_t_24", "s_t_24", "t_c_6", "t_c_6", "t_t_6", "t_t_6","t_c_24", "t_c_24", "t_t_24", "t_t_24")
dge_fr <- DGEList(counts = counts_fr, group = group_fr)
keep <- filterByExpr(dge_fr)
dge_fr <- dge_fr[keep, , keep.lib.sizes=FALSE]
dge_fr <- normLibSizes(dge_fr)


design <- model.matrix(~0+group, data=dge$samples)
colnames(design) <- levels(dge$samples$group)
attr(design,"assign")
attr(design,"contrasts")
attr(design,"contrasts")$group

design_fr <- model.matrix(~0+group, data=dge_fr$samples)
colnames(design_fr) <- levels(dge_fr$samples$group)
attr(design,"assign")
attr(design,"contrasts")
attr(design,"contrasts")$group

#it is recommended to use estimateDisp(y, design)

dge_disp <- estimateDisp(dge, design) #common dispersion, trended dispersions, tagwise dispersions
dge_dispfr <- estimateDisp(dge_fr, design_fr)

dgelist <- estimateGLMCommonDisp(dgelist, design) #common dispersion
dgelist <- estimateGLMTrendedDisp(dgelist, design) #trended dispersions
dgelist <- estimateGLMTagwiseDisp(dgelist, design) #tagwise dispersions


#We can proceed with fitting negative binomial generalized linear models and testing procedures for 
#determining differential expression using either quasi-likelihood (QL) F-test.

fit <- glmQLFit(dge, design, dispersion = bcv^2) #fitting the glm model if you have no replicates. doesn't work very well
fit_fr <- glmQLFit(dge_dispfr, design_fr) #fitting the model with replicates

#Then, if you have multiple groups to compare you can create a contrasts-variable that contains, all the
# information on what comparissons you want to make:

contrasts <- makeContrasts(Susc_6h=s_t_6-s_c_6, Susc_24h=s_t_24-s_c_24, Tol_6h=t_t_6-t_c_6, Tol_24h=t_t_24-t_c_24, levels=design_fr)

#After that, you can get the DE genes using the glmQLFtest() and topTags()
qlf.Susc_6h <- glmQLFTest(fit_fr, contrast=contrasts[,"Susc_6h"])
qlf.Susc_24h <- glmQLFTest(fit_fr, contrast=contrasts[,"Susc_24h"])

qlf.Tol_6h <- glmQLFTest(fit_fr, contrast=contrasts[,"Tol_6h"])
qlf.Tol_24h <- glmQLFTest(fit_fr, contrast=contrasts[,"Tol_24h"])

DE_Susc_6h <- as.data.frame(topTags(qlf.Susc_6h, n = 25879, adjust.method = "BH", sort.by = "PValue", p.value = 0.05))
DE_Susc_24h <- as.data.frame(topTags(qlf.Susc_24h, n = 25879, adjust.method = "BH", sort.by = "PValue", p.value = 0.05))

DE_Tol_6h <- as.data.frame(topTags(qlf.Tol_6h, n = 25879, adjust.method = "BH", sort.by = "PValue", p.value = 0.05))
DE_Tol_24h <- as.data.frame(topTags(qlf.Tol_24h, n = 25879, adjust.method = "BH", sort.by = "PValue", p.value = 0.05))

write.table(DE_Susc_6h, "E:/Manja/edgeR/results/fr_DE_Susc_6h.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(DE_Susc_24h, "E:/Manja/edgeR/results/fr_DE_Susc_24h.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(DE_Tol_6h, "E:/Manja/edgeR/results/fr_DE_Tol_6h.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(DE_Tol_24h, "E:/Manja/edgeR/results/fr_DE_Tol_24h.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

#You can also do this step by step

# If you have an experiment with more factors (groups), this has to be taken into consideration 
# when writing the function. The default is to compare everything to the baseline of group 1, for example:

qlf.2vs1 <- glmQLFTest(fit, coef=2) # this compares the second to the first group
qlf.3vs1 <- glmQLFTest(fit, coef=3) # this compares the third to the first group

#If you want to compare different groups with each other, not with the first one:

qlf.3vs2 <- glmQLFTest(fit, contrast=c(0,-1,1)) # this compares the third to second group

# To find genes different between any of the three groups:

qlf <- glmQLFTest(fit, coef=2:3)
topTags(qlf)

#Additionaly, you can do this with a more rigurous TREAT statistical test. 

glmtreat_Susc_6h <- glmTreat(fit_fr, contrast=contrasts[,"Susc_6h"], coef="s_t_6", lfc=1)
glmtreat_DE_Susc_6h <- as.data.frame(topTags(glmtreat_Susc_6h, n = 25879, adjust.method = "BH", sort.by = "PValue", p.value = 0.05))
write.table(glmtreat_DE_Susc_6h, "E:/Manja/edgeR/results/glmtreat_DE_Susc_6h.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
#########################################################################################################
######################################### PAIRED SAMPLES ################################################
#########################################################################################################

# Paired samples occur whenever we compare two treatments and each independent subject
# in the experiment receives both treatments.

#create a design matrix
info <- read.delim("E:/Manja/edgeR/data/info.txt", header=TRUE, row.names="Sample")
Subject <- factor(info$Subject)
Treatment <- factor(info$Treatment, levels=c("C-6", "T-6", "C-24", "T-24"))
design_paired <- model.matrix(~Subject+Treatment)

dgelist <- DGEList(counts_htseq, design_paired) #???

##########################################################################################################
##########################################################################################################

#Differential expression above a fold-change threshold
qlf <- glmQLFTest(fit, coef=2)
tr <- glmTreat(fit, 
topTags(tr)

#########################################################################################################
#Volcano plot
#########################################################################################################

#The first step is not necessary - it is to select genes whose names will show on the plot

upreg <- c('Zm00001eb027780','Zm00001eb160010','Zm00001eb397030', 'Zm00001eb389190', 'Zm00001eb344730')
downreg <- c('Zm00001eb087980', 'Zm00001eb283460', 'Zm00001eb145080', 'Zm00001eb012570', 'Zm00001eb357740', 'Zm00001eb229190','Zm00001eb065210', 'Zm00001eb383440', 'Zm00001eb208380', 'Zm00001eb272040')

#Next few steps are responsible for setting the rules for coloring the points on the plot

keyvals <- ifelse(DE_Susc_6h$logFC < -1 & DE_Susc_6h$FDR < 10e-4, 'darkred',
  ifelse(DE_Susc_6h$logFC > 1 & DE_Susc_6h$FDR < 10e-4, 'darkgreen','grey')) 

# this tells R that genes with logFC higher than 1 and low enough FDR values will be colored in darkgreen; 
# genes with logFC lower than 1 and low enough FDR values will be colored in darkred, and the rest will be grey

names(keyvals)[keyvals == 'darkgreen'] <- 'Up-regulated'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'darkred'] <- 'Down-regulated'


EnhancedVolcano(DE_Susc_6h, lab = rownames(DE_Susc_6h),  
                x = 'logFC', y ='FDR', selectLab = c(upreg, downreg),
                pCutoff = 10e-4, ylim = c(0, -log10(10e-15)),
                title = 'Susceptible genotype, 6h',
                xlim = c(-10, 10),legendLabSize = 15,legendIconSize = 7,
                colCustom = keyvals, colAlpha = 1, pointSize = 5, 
                labSize = 7,legendPosition = 'left', shape = c (18, 15, 17, 16), titleLabSize = 25)


##########################################################################################################

upreg <- c('Zm00001eb409750', 'Zm00001eb272230', 'Zm00001eb187130', 'Zm00001eb188750', 'Zm00001eb274870',
           'Zm00001eb168900', 'Zm00001eb006310', 'Zm00001eb031210', 'Zm00001eb284350')

downreg <- c('Zm00001eb339420', 'Zm00001eb371530', 'Zm00001eb056970', 'Zm00001eb394330', 'Zm00001eb340550', 
             'Zm00001eb382310', 'Zm00001eb053150', 'Zm00001eb056110', 'Zm00001eb047740', 'Zm00001eb054810')
            

keyvals <- ifelse(DE_Susc_24h$logFC < -1 & DE_Susc_24h$FDR < 10e-4, 'darkred',
                 ifelse(DE_Susc_24h$logFC > 1 & DE_Susc_24h$FDR < 10e-4, 'darkgreen','grey')) 

names(keyvals)[keyvals == 'darkgreen'] <- 'Up-regulated'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'darkred'] <- 'Down-regulated'


EnhancedVolcano(DE_Susc_24h, lab = rownames(DE_Susc_24h),  
                x = 'logFC', y ='FDR', selectLab = c(upreg, downreg),
                pCutoff = 10e-4, ylim = c(0, -log10(10e-15)),
                title = 'Susceptible genotype, 24h',
                xlim = c(-10, 10),legendLabSize = 15,legendIconSize = 7,
                colCustom = keyvals, colAlpha = 1, pointSize = 5, 
                labSize = 7,legendPosition = 'left', shape = c (18, 15, 17, 16), titleLabSize = 25)

##########################################################################################################

upreg <- c('Zm00001eb397030', 'Zm00001eb099460', 'Zm00001eb250940', 'Zm00001eb011990', 'Zm00001eb087860', 
           'Zm00001eb284770', 'Zm00001eb072870')

downreg <- c('Zm00001eb401330', 'Zm00001eb110240', 'Zm00001eb409030', 'Zm00001eb229190', 
             'Zm00001eb357740', 'Zm00001eb067380', 'Zm00001eb159020')


keyvals <- ifelse(DE_Tol_6h$logFC < -1 & DE_Tol_6h$FDR < 10e-4, 'darkred',
                  ifelse(DE_Tol_6h$logFC > 1 & DE_Tol_6h$FDR < 10e-4, 'darkgreen','grey')) 

names(keyvals)[keyvals == 'darkgreen'] <- 'Up-regulated'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'darkred'] <- 'Down-regulated'


EnhancedVolcano(DE_Tol_6h, lab = rownames(DE_Tol_6h),  
                x = 'logFC', y ='FDR', selectLab = c(upreg, downreg),
                pCutoff = 10e-4, ylim = c(0, -log10(10e-15)),
                title = 'Tolerant genotype, 6h',
                xlim = c(-10, 10),legendLabSize = 15,legendIconSize = 7,
                colCustom = keyvals, colAlpha = 1, pointSize = 5, 
                labSize = 7,legendPosition = 'left', shape = c (18, 15, 17, 16), titleLabSize = 25)

##########################################################################################################

upreg <- c('Zm00001eb187130', 'Zm00001eb221790', 'Zm00001eb356910', 'Zm00001eb006310',
           'Zm00001eb201580', 'Zm00001eb027960', 'Zm00001eb113780')

downreg <- c('Zm00001eb027780', 'Zm00001eb200910', 'Zm00001eb293270', 'Zm00001eb009940', 'Zm00001eb382310',
             'Zm00001eb186620', 'Zm00001eb238440', 'Zm00001eb004670', 'Zm00001eb291720', 'Zm00001eb056970')
             


keyvals <- ifelse(DE_Tol_24h$logFC < -1 & DE_Tol_24h$FDR < 10e-4, 'darkred',
                  ifelse(DE_Tol_24h$logFC > 1 & DE_Tol_24h$FDR < 10e-4, 'darkgreen','grey')) 

names(keyvals)[keyvals == 'darkgreen'] <- 'Up-regulated'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'darkred'] <- 'Down-regulated'


EnhancedVolcano(DE_Tol_24h, lab = rownames(DE_Tol_24h),  
                x = 'logFC', y ='FDR', selectLab = c(upreg, downreg),
                pCutoff = 10e-4, ylim = c(0, -log10(10e-15)),
                title = 'Tolerant genotype, 24h',
                xlim = c(-10, 10),legendLabSize = 15,legendIconSize = 7,
                colCustom = keyvals, colAlpha = 1, pointSize = 5, 
                labSize = 7,legendPosition = 'left', shape = c (18, 15, 17, 16), titleLabSize = 25)

##########################################################################################################

upreg <- c('Zm00001eb292880', 'Zm00001eb317310', 'Zm00001eb391290', 'Zm00001eb160010', 'Zm00001eb304830', 
           'Zm00001eb352960', 'Zm00001eb008040', 'Zm00001eb072050', 'Zm00001eb253360', 'Zm00001eb160010',
           'Zm00001eb188570')

downreg <- c('Zm00001eb210320', 'Zm00001eb415430', 'Zm00001eb141650', 'Zm00001eb267220', 'Zm00001eb029940', 
             'Zm00001eb259830', 'Zm00001eb364840', 'Zm00001eb238830', 'Zm00001eb115210', 'Zm00001eb422590')

keyvals <- ifelse(DE_c_6h$logFC < -1 & DE_c_6h$FDR < 10e-4, 'darkred',
                  ifelse(DE_c_6h$logFC > 1 & DE_c_6h$FDR < 10e-4, 'darkgreen','grey')) 

names(keyvals)[keyvals == 'darkgreen'] <- 'Up-regulated'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'darkred'] <- 'Down-regulated'


EnhancedVolcano(DE_c_6h, lab = rownames(DE_c_6h),  
                x = 'logFC', y ='FDR', selectLab = c(upreg, downreg),
                pCutoff = 10e-4, ylim = c(0, -log10(10e-27)),
                title = 'Control, 6h',
                xlim = c(-15, 15),legendLabSize = 15,legendIconSize = 7,
                colCustom = keyvals, colAlpha = 1, pointSize = 5, 
                labSize = 7,legendPosition = 'left', shape = c (18, 15, 17, 16), titleLabSize = 25)
                
##########################################################################################################

upreg <- c('Zm00001eb265350', 'Zm00001eb330320', 'Zm00001eb367250', 'Zm00001eb332200', 'Zm00001eb291040', 
           'Zm00001eb374980', 'Zm00001eb317310', 'Zm00001eb168900', 'Zm00001eb238030', 'Zm00001eb188570',
           'Zm00001eb160010')

downreg <- c('Zm00001eb293030', 'Zm00001eb210320', 'Zm00001eb415430', 'Zm00001eb410980', 'Zm00001eb238830', 
             'Zm00001eb364840', 'Zm00001eb221790', 'Zm00001eb011990', 'Zm00001eb266530', 'Zm00001eb115210', 
             'Zm00001eb422590', 'Zm00001eb077780')

keyvals <- ifelse(DE_c_24h$logFC < -1 & DE_c_24h$FDR < 10e-4, 'darkred',
                  ifelse(DE_c_24h$logFC > 1 & DE_c_24h$FDR < 10e-4, 'darkgreen','grey')) 

names(keyvals)[keyvals == 'darkgreen'] <- 'Up-regulated'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'darkred'] <- 'Down-regulated'


EnhancedVolcano(DE_c_24h, lab = rownames(DE_c_24h),  
                x = 'logFC', y ='FDR', selectLab = c(upreg, downreg),
                pCutoff = 10e-4, ylim = c(0, -log10(10e-30)),
                title = 'Control, 24h',
                xlim = c(-16, 14),legendLabSize = 15,legendIconSize = 7,
                colCustom = keyvals, colAlpha = 1, pointSize = 5, 
                labSize = 7,legendPosition = 'left', shape = c (18, 15, 17, 16), titleLabSize = 25)

##########################################################################################################

upreg <- c('Zm00001eb168900', 'Zm00001eb396560', 'Zm00001eb178570', 'Zm00001eb299830', 'Zm00001eb317310', 
           'Zm00001eb304830', 'Zm00001eb423520', 'Zm00001eb301080', 'Zm00001eb241410', 'Zm00001eb292880')

downreg <- c('Zm00001eb006490', 'Zm00001eb293030', 'Zm00001eb410970', 'Zm00001eb410980', 'Zm00001eb2611377220', 
             'Zm00001eb392900', 'Zm00001eb415430', 'Zm00001eb238830', 'Zm00001eb115210', 'Zm00001eb221790')


keyvals <- ifelse(DE_t_6h$logFC < -1 & DE_t_6h$FDR < 10e-4, 'darkred',
                  ifelse(DE_t_6h$logFC > 1 & DE_t_6h$FDR < 10e-4, 'darkgreen','grey')) 

names(keyvals)[keyvals == 'darkgreen'] <- 'Up-regulated'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'darkred'] <- 'Down-regulated'


EnhancedVolcano(DE_t_6h, lab = rownames(DE_t_6h),  
                x = 'logFC', y ='FDR', selectLab = c(upreg, downreg),
                pCutoff = 10e-4, ylim = c(0, -log10(10e-30)),
                title = 'Treatment, 6h',
                xlim = c(-16, 16),legendLabSize = 15,legendIconSize = 7,
                colCustom = keyvals, colAlpha = 1, pointSize = 5, 
                labSize = 7,legendPosition = 'left', shape = c (18, 15, 17, 16), titleLabSize = 25)

##########################################################################################################

upreg <- c('Zm00001eb265350', 'Zm00001eb292880', 'Zm00001eb291040', 'Zm00001eb391290', 'Zm00001eb047740', 
           'Zm00001eb072050', 'Zm00001eb044850', 'Zm00001eb152140', 'Zm00001eb395670', 'Zm00001eb374980',
           'Zm00001eb160010')

downreg <- c('Zm00001eb259830','Zm00001eb235080', 'Zm00001eb410970', 'Zm00001eb267220', 'Zm00001eb410980', 'Zm00001eb266530', 'Zm00001eb029940', 'Zm00001eb238830', 'Zm00001eb115210', 'Zm00001eb422590')

keyvals <- ifelse(DE_t_24h$logFC < -1 & DE_t_24h$FDR < 10e-4, 'darkred',
                  ifelse(DE_t_24h$logFC > 1 & DE_t_24h$FDR < 10e-4, 'darkgreen','grey')) 

names(keyvals)[keyvals == 'darkgreen'] <- 'Up-regulated'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'darkred'] <- 'Down-regulated'


EnhancedVolcano(DE_t_24h, lab = rownames(DE_t_24h),  
                x = 'logFC', y ='FDR', selectLab = c(upreg, downreg),
                pCutoff = 10e-4, ylim = c(0, -log10(10e-30)),
                title = 'Treatment, 24h',
                xlim = c(-16, 16),legendLabSize = 15,legendIconSize = 7,
                colCustom = keyvals, colAlpha = 1, pointSize = 5, 
                labSize = 7,legendPosition = 'left', shape = c (18, 15, 17, 16), titleLabSize = 25)
############################################################################################################
#PCA plot
############################################################################################################

#There is a separate script :)

############################################################################################################
#Cluster analysis and heatmaps
############################################################################################################

#There is a separate script :)

############################################################################################################
#Venn's diagram
############################################################################################################

library(VennDiagram)
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
library(ggvenn)
library(ggVenDiagram)
library(gplots)
library(VennDetail)

setwd("E:/Manja/edgeR/results")

susc.6h = readLines("genes_susc.6h.txt")
susc.24h = readLines("genes_susc.24h.txt")
tol.6h = readLines("genes_tol.6h.txt")
tol.24h = readLines("genes_tol.24h.txt")

control.6h = readLines("genes_control.6h.txt")
control.24h = readLines("genes_control.24h.txt")
treatment.6h = readLines("genes_treatment.6h.txt")
treatment.24h = readLines("genes_treatment.24h.txt")

t_vs_c <-  venndetail(list('Susc.genotype, 6h' = susc.6h, 'Susc.genotype, 24h' = susc.24h, 'Tol.genotype, 6h' = tol.6h, 'Tol.genotype, 24h' = tol.24h))

tol_vs_susc <- venndetail(list('Control, 6h' = control.6h, 'Control, 24h' = control.24h, 'Treatment, 6h' = treatment.6h, 'Treatment, 24h' = treatment.24h))

plot(t_vs_c, type = "venn", col = c("#89c2fa", "#ffd844", "#ffac2a", "#b8e17f"),
     mycol = c("#89c2fa", "#ffd844", "#ffac2a", "#b8e17f"),
     sep = "_", cat.cex = 2, alpha = 0.1, cex = 2,
     cat.fontface = "bold", margin = 0.05, text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5), 
     filename = NULL, piecolor = NULL,
     revcolor = "lightgrey", any = NULL, show.number = TRUE,
     show.x = TRUE, log = FALSE, base = NULL, percentage = FALSE,
     sets.x.label = "Set Size", mainbar.y.label = "Intersection Size",
     nintersects = 40, abbr = FALSE, abbr.method = "both.sides",
     minlength = 3)

dev.off()

plot(tol_vs_susc, type = "venn", col = c("#89c2fa", "#ffd844", "#ffac2a", "#b8e17f"),
     mycol = c("#89c2fa", "#ffd844", "#ffac2a", "#b8e17f"),
     sep = "_", cat.cex = 2, alpha = 0.1, cex = 2,
     cat.fontface = "bold", margin = 0.05, text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5), 
     filename = NULL, piecolor = NULL,
     revcolor = "lightgrey", any = NULL, show.number = TRUE,
     show.x = TRUE, log = FALSE, base = NULL, percentage = FALSE,
     sets.x.label = "Set Size", mainbar.y.label = "Intersection Size",
     nintersects = 40, abbr = FALSE, abbr.method = "both.sides",
     minlength = 3)

dev.off()


############################################################################################################
#Gene onthology
############################################################################################################


#Gene onthology and KEGG enrichment analysis


