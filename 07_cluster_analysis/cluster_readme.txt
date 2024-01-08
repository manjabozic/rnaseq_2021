Data:
- data.xlsx - all the data combined int one excel file
- cleandata.txt - Total data, not normalized 
- data_pca.txt - All the DE genes, selected from the clean data
- de_genes.txt - Selected 241 DE genes, from the normalized data
- norm.data_pca.txt - All the DE genes, selected from the normalized data
- normalized_counts.txt - Total data, normalized in deseq2 median of ratios method (>= 20)

R scripts:
- cluster.R - clustering workflow
- vat.R - attempt of creating a script to do this on the cluster

Images:
- K-Means, k=2/4/8
- K-Medioids, k=2/4/8
- Choosing k - elbow, silhouette, gap statistic method
- Spearman DE
- NBClus choose K
- heatmaps
- dendrograms