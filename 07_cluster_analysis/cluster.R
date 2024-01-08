################################################################################################################
######################################## CLUSTERING & HEATMAPS ################################################
################################################################################################################

# The classification of observations into groups requires some methods for computing the distance or the 
# (dis)similarity between each pair of observations. The result of this computation is known as a 
# dissimilarity or distance matrix.

# The classical methods for distance measures are Euclidean and Manhattan distances.Other dissimilarity 
# measures exist such as correlation-based distances, which is widely used for gene expression data analyses.
# This includes the Pearson correlation distance, Eisen cosine correlation distance, Spearman correlation 
# distance, and Kendall correlation distance.

#Most widely used methods are Euclidean, pearson and Spearman. 

#WORKING DIRECTORY, PACKAGES, DATA

setwd("E:/Manja/cluster") #set the working directory


library("cluster") # load the neccessary packages
library("FactoMineR")
library("factoextra")
library("stats")
library("fpc")
library("baySeq")
library("dendextend")
library("gplots")
library("pheatmap")
library("ComplexHeatmap")
library("hopkins")
library("NbClust")
library("clValid")
library("pvclust")
library("dbscan")
library("clusterSeq")

# You can load already normalized data (this was noralized in deseq2 median of ratios method)

counts_mor <- as.matrix(read.delim("E:/Manja/cluster/normalized_counts.txt", header=TRUE, row.names=1))
#counts_mor <- t(counts_mor)

# this is the data for all DE genes, selected from the normalized counts

de <- as.matrix(read.delim("E:/Manja/cluster/norm.data_pca.txt", header=TRUE, row.names="Gene_ID"))
de <- as.matrix(read.delim("E:/Manja/cluster/de_genes.txt", header=TRUE, row.names="Gene_ID"))
de <-scale(de)
#de <- t(de)

# or you can load raw read counts and normalize them here...

counts <- read.delim("E:/Manja/cluster/cleandata.txt", header=TRUE, row.names=1)

# The value of distance measures is intimately related to the scale on which measurements are made. 
# Therefore, variables are often scaled (i.e. standardized) before measuring the inter-observation 
# dissimilarities.The goal is to make the variables comparable. Generally variables are scaled to have
# i) standard deviation one and ii) mean zero.

#1. Clusterseq package normalization method

ls_all <- getLibsizes(data = counts)
counts[counts == 0] <- 1
norm_counts <- log2(t(t(counts / ls_all)) * mean(ls_all))
#norm_counts <- t(norm_counts)

#2. Scaling the data 

# Since we are comparing gene expression patterns we need to scale the data otherwise all of the highly 
# expressed genes will cluster together even if they have different patterns among the samples. 
# The function scale scales collumns. So to scale per gene (rows) data needs to be transposed.

counts <- t(scale(t(counts)))
counts_mor <- t(scale(t(counts_mor)))
genes <- t(scale(t(de)))
norm_counts <- t(scale(t(norm_counts)))

#rev_de <- t(de)

################################################################################################################
#Clustering tendency Validation
################################################################################################################

# Clustering tendency can be evaluated in multiple ways: visual (pca plot, VAT algorithm) and statistical. PCA plots are 
# explained in a separate script. Most common statistical method for evaluating the clustering tendency 
# is the Hopkins statistic.

set.seed(123)

hopkins(counts, m = nrow(counts)-1, d = ncol(counts), k = 1, method = "simple") # H=0.9962365
hopkins(counts_mor, m = nrow(counts_mor)-1, d = ncol(counts_mor), k = 1, method = "simple") # H=0.9958561
hopkins(norm_counts, m = nrow(norm_counts)-1, d = ncol(norm_counts), k = 1, method = "simple") # H=0.9999993
hopkins(de, m = nrow(de)-1, d = ncol(de), k = 1, method = "simple") # H=0.9998826

res <- get_clust_tendency(counts, n = nrow(counts)-1, graph = FALSE)
res_mor <- get_clust_tendency(counts_mor, n = nrow(counts_mor)-1, graph = FALSE)
norm_res <- get_clust_tendency(norm_counts, n = nrow(norm_counts)-1, graph = FALSE)
res_de <- get_clust_tendency(de, n = nrow(de)-1, graph = FALSE)

#If the samples were uniformly distributed, and H would be about 0.5. However, if clusters are present, then 
# the distances for artificial points would be substantially larger than for the real ones, and thus the value 
# of H will increase. A value for H higher than 0.75 indicates a clustering tendency at the 90% confidence 
# level.

# Visual assessment of cluster tendency (VAT) approach

fviz_dist(dist(counts), show_labels = FALSE) + labs(title = "Scaled, readcount")
fviz_dist(dist(norm_counts), show_labels = FALSE) + labs(title = "Normalized, log2, readcount")
fviz_dist(dist(counts_mor), show_labels = FALSE) + labs(title = "Mean of ratios normalization method, readcount")

fviz_dist(dist(de), show_labels = FALSE) + labs(title = "Mean of ratios normalization method, DE readcount")

################################################################################################################
#Determining the number of clusters
################################################################################################################
#Elbow, Silhouhette and Gap statistic methods

#Elbow method
fviz_nbclust(counts, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)+ labs(subtitle = "Elbow method")
fviz_nbclust(norm_counts, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)+ labs(subtitle = "Elbow method")
fviz_nbclust(counts_mor, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)+ labs(subtitle = "Elbow method")
fviz_nbclust(de, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)+ labs(subtitle = "Elbow method")

# Silhouette method 
fviz_nbclust(counts, kmeans, method = "silhouette")+labs(subtitle = "Silhouette method")
fviz_nbclust(norm_counts, kmeans, method = "silhouette")+labs(subtitle = "Silhouette method")
fviz_nbclust(counts_mor, kmeans, method = "silhouette")+labs(subtitle = "Silhouette method")
fviz_nbclust(de, kmeans, method = "silhouette")+labs(subtitle = "Silhouette method")

# Gap statistic
set.seed(123)
fviz_nbclust(counts, kmeans, nstart = 25, method = "gap_stat", nboot = 500) + labs(subtitle = "Gap statistic method")
fviz_nbclust(norm_counts, kmeans, nstart = 25, method = "gap_stat", nboot = 500) + labs(subtitle = "Gap statistic method")
fviz_nbclust(counts_mor, kmeans, nstart = 25, method = "gap_stat", nboot = 500) + labs(subtitle = "Gap statistic method")
fviz_nbclust(de, kmeans, nstart = 25, method = "gap_stat", nboot = 500) + labs(subtitle = "Gap statistic method")

# The disadvantage of elbow and average silhouette methods is that, they measure a global clustering 
# characteristic only. A more sophisticated method is to use the gap statistic which provides a statistical 
# procedure to formalize the elbow/silhouette heuristic in order to estimate the optimal number of clusters.

nbclust <- NbClust(de, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans")

#distance: the distance measure to be used to compute the dissimilarity matrix. Possible values include 
# “euclidean”, “manhattan” or “NULL

################################################################################################################
#Distance matrix computation
################################################################################################################

# The classification of observations into groups requires some methods for computing the distance or the 
# (dis)similarity between each pair of observations. The result of this computation is known as a 
# dissimilarity or distance matrix.

# Methods for creating the distance matrix include classical (Euclidean, Manhattan) and correlation-based 
# methods (Pearson, Kendall or Spearman correlation distance).Correlation-based distances are widely used for gene 
# expression data analyses

#Euclidean
e_dist <- dist(counts, method = "euclidean")
e_dist_norm <- dist(norm_counts, method = "euclidean")
e_dist_mor <- dist(counts_mor, method = "euclidean")
e_dist_de <- dist(de, method = "euclidean")

#Pearson
p_dist <- get_dist(counts, method = "pearson")
p_dist_norm <- get_dist(norm_counts, method = "pearson")
p_dist_mor <- get_dist(counts_mor, method = "pearson")
p_dist_de <- get_dist(de, method = "pearson")

#Speraman
s_dist <- get_dist(counts, method = "spearman")
s_dist_norm <- get_dist(norm_counts, method = "spearman")
s_dist_mor <- get_dist(counts_mor, method = "spearman")
s_dist_de <- get_dist(de, method = "spearman")

#rev_e_dist_de <- dist(rev_de, method = "euclidean")
#rev_p_dist_de <- get_dist(rev_de, method = "pearson")
#rev_s_dist_de <- get_dist(rev_de, method = "spearman")

# Visualisation

fviz_dist(e_dist_de, show_labels = FALSE)
fviz_dist(p_dist_de, show_labels = FALSE)
fviz_dist(s_dist_de, show_labels = FALSE)

fviz_dist(e_dist_mor, show_labels = FALSE)
fviz_dist(p_dist_mor, show_labels = FALSE)
fviz_dist(s_dist_mor, show_labels = FALSE)

################################################################################################################
# Partitioning Clustering
################################################################################################################

# Partitioning clustering are clustering methods used to classify observations, within a data set, into 
# multiple groups based on their similarity. The algorithms require the analyst to specify the number of 
# clusters to be generated.

# K-Means Clustering

# K-means clustering (MacQueen, 1967) is the most commonly used unsupervised machine learning algorithm for 
# partitioning a given data set into a set of k groups (i.e.k clusters), where k represents the number of 
# groups pre-specified by the analyst.

km4 <- kmeans(counts, 4, nstart = 100) # numerix matrix first, number of clusters
km4_norm <- kmeans(norm_counts, 4, nstart = 100) #nstart - number of random starting partitions when centers is a number
km4_mor <- kmeans(counts_mor, 4, nstart = 1000)

km2 <- kmeans(de, 2, nstart = 500)
km4 <- kmeans(de, 4, nstart = 500)
km8 <- kmeans(de, 8, nstart = 500)

# Available components of the kmeans file are: 
# "cluster" - vector of integers (from 1:k) indicating the cluster to which each point is allocated
# "centers" - matrix of cluster centers (cluster means)
# "totss" - total sum of squares (TSS), which measures the total variance in the data
# "withinss"  - vector of within-cluster sum of squares, one component per cluster
# "tot.withinss" - total within-cluster sum of squares, i.e. sum(withinss)
# "betweenss" - the between-cluster sum of squares, i.e. totss ≠ tot.withinss
# "size" - the number of observations in each cluster
# "iter" 
# "ifault"

print(km4) #prints all the results of the kmeans function
aggregate(de, by=list(cluster=km4$cluster), mean) #computes the mean of each variables by clusters using the OG data
km4$cluster

# Visualization of the kmeans results

fviz_cluster(km2, data = de, geom = c("point"), 
  palette = c("#00AFBB", "#E7B800"), 
  main = "K-means method, k=2", pointsize = 1.0,labelsize = 12, 
  ellipse.type = "euclid", star.plot = TRUE, repel = TRUE,ggtheme = theme_minimal())

fviz_cluster(km4, data = de, geom = c("point"), 
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#6aa84f"), 
             main = "K-means method, k=4", pointsize = 1.0,labelsize = 12, 
             ellipse.type = "euclid", star.plot = TRUE, repel = TRUE,ggtheme = theme_minimal())

fviz_cluster(km8, data = de, geom = c("point"), 
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", "#6aa84f", "#dbf259", "#cd3f85", "#674ea7"), 
             main = "K-means method, k=8", pointsize = 1.0,labelsize = 12, 
             ellipse.type = "euclid", star.plot = TRUE, repel = TRUE,ggtheme = theme_minimal())

# K-Medoids

# The k-medoids algorithm is a clustering approach related to k-means clustering  for partitioning a data set
# into k groups or clusters. In k-medoids clustering, each cluster is represented by one of the data point in 
# the cluster. These points are named cluster medoids. The most common k-medoids clustering methods is the 
# PAM algorithm (Partitioning Around Medoids, Kaufman & Rousseeuw, 1990).

# K-medoid is a robust alternative to k-means clustering. This means that, the algorithm is less sensitive 
# to noise and outliers, compared to k-means, because it uses medoids as cluster centers instead of means 
# (used in k-means).

# Input files for pam can be the count files (data matrix or numeric data frame) or a dissimilarity matrix

pam4 <- pam(counts, 4)
pam4_norm <- pam(norm_counts, 4)
pam4_mor <- pam(counts_mor, 4)
pam4_de <- pam(de, 4)

pam2 <- pam(de, 2, metric = "euclidean", stand = FALSE) #counts file, cluster no, no standardization
pam4 <- pam(de, 4, metric = "euclidean", stand = FALSE)
pam8 <- pam(de, 8, metric = "euclidean", stand = FALSE)


fviz_cluster(pam2, data = de, geom = c("point"), 
             palette = c("#00AFBB", "#E7B800"), 
             main = "K-medioids, PAM method, k=2", pointsize = 1.0,labelsize = 12, 
             ellipse.type = "euclid", star.plot = TRUE, repel = TRUE,ggtheme = theme_minimal())

fviz_cluster(pam4, data = de, geom = c("point"), 
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#6aa84f"), 
             main = "K-medioids, PAM method, k=4", pointsize = 1.0,labelsize = 12, 
             ellipse.type = "euclid", star.plot = TRUE, repel = TRUE,ggtheme = theme_minimal())

fviz_cluster(pam8, data = de, geom = c("point"), 
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", "#6aa84f", "#dbf259", "#cd3f85", "#674ea7"), 
             main = "K-medioids, PAM method, k=8", pointsize = 1.0,labelsize = 12, 
             ellipse.type = "euclid", star.plot = TRUE, repel = TRUE,ggtheme = theme_minimal())


################################################################################################################
# Hierarchical clustering
################################################################################################################

# Hierarchical clustering (HCA) is an alternative approach to partitioning clustering for grouping objects 
# based on their similarity. In contrast to partitioning clustering, hierarchical clustering does not
# require you to pre-specify the number of clusters to be produced. Hierarchical clustering can be subdivided 
# into two types: agglomerative and divise clustering. In agglomerative clustering, each observation is 
# initially considered as a cluster of its own (leaf). Then, the most similar clusters are successively 
# merged until there is just one single big cluster (root).
# Divise clustering begins with the root, in witch all objects are included in one cluster. Then the most 
# heterogeneous clusters are successively divided until all observation are in their own cluster. The result 
# of hierarchical clustering is a tree-based representation of the objects, also known as dendrogram.


# Agglomerative Clustering - agglomerative nesting (AGNES):

# Preparing the data
# Computing (dis)similarity information between every pair of objects in the data set
# Using linkage function to group objects into hierarchical cluster tree, based on the distance information 
# generated in the previous step. 
# Determining where to cut the hierarchical tree into clusters. This creates a partition of the data.

# While using the hclust() function, you can choose the linkage method:

# Maximum or complete linkage: The distance between two clusters is defined as the maximum value of all 
# pairwise distances between the elements in cluster 1 and the elements in cluster 2. It tends to produce 
# more compact clusters.

# Minimum or single linkage: The distance between two clusters is defined as the minimum value of all 
# pairwise distances between the elements in cluster 1 and the elements in cluster 2. It tends to produce 
# long, "loose" clusters.

# Mean or average linkage: The distance between two clusters is defined as the average distance between the 
# elements in cluster 1 and the elements in cluster 2.

# Centroid linkage: The distance between two clusters is defined as the distance between the centroid for 
# cluster 1 (a mean vector of length p variables) and the centroid for cluster 2.

# Ward’s minimum variance method: It minimizes the total within-cluster variance. At each step the pair of 
# clusters with minimum between-cluster distance are merged.

# The data is already loaded and standardized, and the distance matrices are created.

de_e.hc <- hclust(d = e_dist_de, method = "ward.D2")
de_p.hc <- hclust(d = p_dist_de, method = "ward.D2")
de_s.hc <- hclust(d = s_dist_de, method = "ward.D2")

genes_p.hc <- hclust(d = p_dist_genes, method = "ward.D2")
genes_s.hc <- hclust(d = s_dist_genes, method = "ward.D2")

rev_de_e.hc <- hclust(d = rev_e_dist_de, method = "ward.D2")
rev_de_p.hc <- hclust(d = rev_p_dist_de, method = "ward.D2")
rev_de_s.hc <- hclust(d = rev_s_dist_de, method = "ward.D2")

com_de_e.hc <- hclust(d = e_dist_de, method = "complete")
com_de_p.hc <- hclust(d = p_dist_de, method = "complete")
com_de_s.hc <- hclust(d = s_dist_de, method = "complete")

# Create a dendogram
# Dendrograms correspond to the graphical representation of the hierarchical tree. 

# The height of the fusion, provided on the vertical axis, indicates the (dis)similarity/distance between 
# two objects/clusters. The higher the height of the fusion, the less similar the objects are. This height 
# is known as the cophenetic distance between the two objects.

fviz_dend(de_e.hc, cex = 0.5, main = "Dendrogram - ward.D2",
          xlab = "Objects", ylab = "Distance", sub = "", horiz = TRUE)
fviz_dend(de_p.hc, cex = 0.5, main = "Dendrogram - ward.D2",
          xlab = "Objects", ylab = "Distance", sub = "", horiz = TRUE)
fviz_dend(de_s.hc, cex = 0.5)

colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", 
           "#6aa84f", "#cd3f85", "#674ea7", "#ea9999",
           "#dbf259", "#20124d", "#660000","#9795f0" )

col4 = c("#00AFBB", "#E7B800", "#FC4E07", "#6aa84f")

fviz_dend(rev_de_p.hc, k=4, cex = 1, 
          k_colors = col4,
          show_labels = TRUE,
          color_labels_by_k = TRUE, lwd = 1,
          rect = TRUE, rect_fill = TRUE, main = "Dendrogram - ward.D2",
          xlab = "Objects", ylab = "Distance", horiz = TRUE)

fviz_dend(rev_de_p.hc, k=4, cex = 1, 
          k_colors = col4, type = "phylogenic", phylo_layout = "layout_with_lgl",
          show_labels = TRUE,
          color_labels_by_k = TRUE, lwd = 1,
          rect = TRUE, rect_fill = TRUE, main = "Dendrogram - ward.D2",
          xlab = "Objects", ylab = "Distance", horiz = TRUE)

require("igraph")
fviz_dend(de_e.hc, k=4, cex = 0.3, 
          k_colors = col4,
          show_labels = TRUE,
          color_labels_by_k = TRUE, lwd = 0.5,
          rect = TRUE, rect_fill = TRUE, main = "Dendrogram - ward.D2",
          xlab = "Objects", ylab = "Distance", horiz = TRUE, repel = TRUE,
          )

# You can zoom in on certain parts of the dendogram using xlim and ylim; the graph cannot be horizontal for this to work 

fviz_dend(de_e.hc, k=4, cex = 0.3, 
          k_colors = col4, xlim = c(1, 3), ylim = c(1, 20),
          show_labels = TRUE,
          color_labels_by_k = TRUE, lwd = 0.5,
          rect = TRUE, rect_fill = TRUE, main = "Dendrogram - ward.D2",
          xlab = "Objects", ylab = "Distance", repel = TRUE)

# Plotting a sub-tree of dendrograms  - first create the whole dendrogram using fviz_dend() and save the result into an
# object, and then use cut.dendrogram() to cut the dendrogram, at a given height (h), into multiple sub-trees

de_p_dend.plot <- fviz_dend(rev_de_p.hc, k=4, cex = 1, # Create the whole dendrogram
                            k_colors = col4, type = "phylogenic", phylo_layout = "layout_with_lgl",
                            show_labels = TRUE, color_labels_by_k = TRUE, lwd = 1,
                            rect = TRUE, rect_fill = TRUE, main = "Dendrogram - ward.D2",
                            xlab = "Objects", ylab = "Distance", horiz = TRUE)

dend_data <- attr(de_p_dend.plot, "dendrogram") #Extract dendrogram data
dend_cuts <- cut(dend_data, h = 10)
fviz_dend(dend_cuts$upper)


# Verify the cluster tree by computing the correlation between the cophenetic distances and the original 
# distance data generated by the dist() function.  The closer the value of the correlation coefficient is 
# to 1, the more accurately the clustering solution reflects your data. Values above 0.75 are felt to be 
# good. The “average” linkage method appears to produce high values of this statistic. T

res.coph <- cophenetic(de_e.hc)
cor(e_dist_de, res.coph) # 0.9824879 and 0.8717125

res.coph <- cophenetic(de_p.hc)
cor(p_dist_de, res.coph) # 0.9972443 and 0.880627

res.coph <- cophenetic(de_s.hc)
cor(s_dist_de, res.coph) # 0.9975393 and 0.8707527

# You can cut the hierarchical tree at a given height in order to partition your data into clusters. 
# The R base function cutree() can be used to cut a tree, generated by the hclust() function, into several 
# groups either by specifying the desired number of groups or the cut height.

grp_e <- cutree(de_e.hc, k = 4)
head(grp_e, n = 4)

grp_p <- cutree(de_p.hc, k = 13) #maybe 13?
head(grp_p, n = 4)

grp_s <- cutree(de_s.hc, k = 4)
head(grp_s, n = 4)

grp_genes <- cutree(genes.hc, k = 12)
head(grp_s, n = 4)

table(grp_e) # get the number of members in each cluster
table(grp_p)
table(grp_s)
table(grp_genes)

rownames(de)[grp_e == 1]# get the names for the members of cluster 1
rownames(de)[grp_p == 1]
rownames(de)[grp_s == 1]

clusters <- data.frame (first_column  = cluster1,
                  second_column = cluster2
)

fviz_dend(rev_de_e.hc, k = 4, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE) # Add rectangle around groups

fviz_cluster(list(data = de, cluster = grp_e),
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
             ellipse.type = "convex", # Concentration ellipse
             repel = TRUE, # Avoid label overplotting (slow)
             show.clust.cent = FALSE, ggtheme = theme_minimal())



# The R package cluster makes it easy to perform cluster analysis in R. It provides the function agnes() and 
# diana() for computing agglomerative and divisive clustering, respectively. These functions perform all the 
# necessary steps for you. You don’t need to execute the scale(), dist() and hclust() function separately.


#Compare two dendograms:

dend1 <- as.dendrogram (de_s.hc) #Create two dendrograms
dend2 <- as.dendrogram (com_de_s.hc)

dend_list <- dendlist(dend1, dend2) #Create a list to hold dendrograms

tanglegram(dend1, dend2) #Draw a tanglegram to compare the two dendograms

cor.dendlist(dend_list, method = "cophenetic") #Correlation matrix between a list of dendrograms

cor_cophenetic(dend1, dend2) # Cophenetic correlation coefficient

# You can also compare multiple dendograms, for example:

dend3 <- as.dendrogram(hclust(d = s_dist_de, method = "average"))
dend4 <- as.dendrogram(hclust(d = s_dist_de, method = "single"))
dend5 <- as.dendrogram(hclust(d = s_dist_de, method = "centroid"))

dend_list <- dendlist("Ward" = dend1, "Complete" = dend2, "Average" = dend3, "Single" = dend4,
                       "Centroid" = dend5)
cors <- cor.dendlist(dend_list)
round(cors, 2)
library(corrplot)
corrplot(cors, "pie", "lower")

################################################################################################################
# Heatmaps
################################################################################################################

# A heatmap (or heat map) is another way to visualize hierarchical clustering. First hierarchical clustering 
# is done of both the rows and the columns of the data matrix. The columns/rows of the data matrix are 
# re-ordered according to the hierarchical clustering result, putting similar observations close to each other.

# There are a multiple numbers of R packages and functions for drawing interactive and static heatmaps, 
# including: heatmap(), heatmap2(), pheatmap(), d3heatmap(), Heatmap()...

col<- colorRampPalette(c("red", "white", "blue"))(256)

# colorRamp2 is recommended because it makes comparing multiple heatmaps possible - the same color always corresponds to 
# the same value

library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun02 = colorRamp2(c(-2, 0, 2), c("darkred", "white", "darkblue"))
col_fun(seq(-3, 3))

heatmap(de, scale = "none", col = col_fun)

heatmap.2(de, scale = "none", col = bluered(100), distfun = get_dist,
          trace = "none", density.info = "none")


Heatmap(de, name = "All DE expressed genes", column_title = "Genotypes and treatments", 
        row_title = "Genes", row_names_gp = gpar(fontsize = 7), 
        clustering_distance_rows = "spearman", clustering_method_columns = "ward.D", row_dend_width = unit(70, "mm"), 
        row_dend_gp = gpar(),show_row_names = FALSE, cluster_columns = TRUE)

Heatmap(genes, name = "C vs. T, DE expressed genes", 
        column_title = "Genotypes and treatments", column_title_gp = gpar(fontsize = 50),clustering_method_columns = "ward.D",
        row_title = "Genes", clustering_distance_rows = "pearson", row_title_gp = gpar(fontsize = 50),
        row_dend_width = unit(30, "mm"),  column_dend_height = unit(5, "mm"),
        column_names_gp = gpar(fontsize = 35),
        row_dend_gp = gpar(), show_row_names = FALSE, cluster_columns = TRUE, cluster_rows = TRUE, col = col_fun02)
