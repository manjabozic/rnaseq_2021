################################################################################################################
################################################ PCA ###########################################################
################################################################################################################

#Principal component analysis (PCA) allows us to summarize and to visualize the information in a data
#set containing individuals/observations described by multiple inter-correlated quantitative variables. 

#Principal component analysis is used to extract the important information from a multivariate data table
#and to express this information as a set of few new variables called principal components.

#Install neccessary packages:

install.packages("FactoMineR")
install.packages("factoextra")

# Load those packages:

library("FactoMineR")
library("factoextra")
library("baySeq")

#Set the working directory:

setwd("E:/Manja/pca")

#Load the data

counts <- read.delim("E:/Manja/pca/cleandata.txt", header=TRUE, row.names=1)
counts <- t(counts)

#you can normalize the counts in DESeq2 and use that instead

library(DESeq2)
counts <- as.matrix(read.delim("E:/Manja/edgeR/data/htseq_name_counts.txt", header=TRUE, row.names="ENSEMBL_GeneID"))

coldata <- read.delim("E:/Manja/pca/coldata.txt", row.names=1)
coldata <- coldata[,c("condition","time_point", "genotype")]
coldata$condition <- factor(coldata$condition) #factor() encodes a vector as a factor
coldata$time_point <- factor(coldata$time_point)
coldata$genotype <- factor(coldata$genotype)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- estimateSizeFactors(dds)
smallestGroupSize <- 1
keep <- rowSums(counts(dds) >= 20) >= smallestGroupSize
dds <- dds[keep,]

normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="E:/Manja/pca/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

norm_counts <- read.delim("E:/Manja/pca/normalized_counts.txt", header=TRUE, row.names=1)
norm_counts <- t(norm_counts)


############################################################################################################
# ALL GENES COUNTED - PCA
#############################################################################################################

res.pca <- PCA(counts, scale.unit = TRUE, graph = FALSE)
res.pca_norm <- PCA(norm_counts, scale.unit = TRUE, graph = FALSE)

print(res.pca)

#############################################################################################################
#Eigenvalues / Variances
#############################################################################################################

#The eigenvalues measure the amount of variation retained by each principal component. 
#Eigenvalues are large for the first PCs and small for the subsequent PCs. That is, the first PCs corresponds 
#to the directions with the maximum amount of variation in the data set. 

#The sum of all the eigenvalues give a total variance.
#The proportion of variation explained by each eigenvalue is given in the second column.

eig.val <- get_eigenvalue(res.pca)
eig.val_norm <- get_eigenvalue(res.pca_norm)

#You can use the eigenvalues to determine the number of components to retain, based on the percentage 
#of variation explained. A scree plot can also be used for that.

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
fviz_eig(res.pca_norm, addlabels = TRUE, ylim = c(0, 50), title = 'Percentage of explained variances', 
         choice = c("variance", "eigenvalue"),
         geom = c("bar", "line"),
         barfill = "#E7B800",
         barcolor = "#FC4E07",
         linecolor = "#00AFBB")

#############################################################################################################
# Graph of variables
#############################################################################################################

var <- get_pca_var(res.pca) 
var_norm <- get_pca_var(res.pca_norm) 

## 1 "$coord" "Coordinates for the variables"
## 2 "$cor" "Correlations between variables and dimensions"
## 3 "$cos2" "Cos2 for the variables"
## 4 "$contrib" "contributions of the variables"

fviz_pca_var(res.pca, col.var = "black") # too many variables
library("corrplot")
corrplot(var_norm$cos2, is.corr=FALSE) # too many variables

#############################################################################################################
# Graph of individuals
#############################################################################################################

ind <- get_pca_ind(res.pca)
ind_norm <- get_pca_ind(res.pca_norm)

fviz_pca_ind(res.pca)

fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) #Avoid text overlapping (slow if many points
             
fviz_pca_ind(res.pca, pointsize = "cos2", pointshape = 21, fill = "#E7B800", repel = TRUE)

fviz_pca_ind(res.pca, col.ind = "cos2", pointsize = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

fviz_pca_ind(res.pca_norm, col.ind = "contrib", pointsize = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, title = 'PCA plot - normalized counts')

fviz_pca_ind(res.pca, col.ind = "cos2", pointsize = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
#############################################################################################################
# Quality of representation
#############################################################################################################

fviz_cos2(res.pca, choice = "ind")
fviz_cos2(res.pca_norm, choice = "ind") 

#############################################################################################################
# Contributions of variables to PCs
#############################################################################################################

fviz_contrib(res.pca_norm, choice = "var", axes = 1, top = 10)
fviz_contrib(res.pca_norm, choice = "var", axes = 2, top = 10)

fviz_contrib(res.pca_norm, choice = "var", axes = 1:2, top = 10)
fviz_pca_var(res.pca_norm, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

set.seed(123)
res.km <- kmeans(var_norm$coord, centers = 8, nstart = 25)
grp <- as.factor(res.km$cluster)

fviz_pca_var(res.pca_norm, col.var = grp,
             palette = c("#0073C2FF", "#EFC000FF", "#fff139", "#41aa5c", "#00AFBB", "#FC4E07","#a46ee4", "#8fce00"),
             legend.title = "Cluster")

#############################################################################################################
# Filtering results
#############################################################################################################

# If you have many individuals/variable, it's possible to visualize only some of them using the arguments
# select.ind and select.var.

fviz_pca_var(res.pca_norm, select.var = list(cos2 = 10))
fviz_pca_var(res.pca_norm, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, 
             select.var = list(cos2 = 20))

fviz_pca_var(res.pca_norm, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, 
             select.var = list(cos2 = 20))


