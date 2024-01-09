Data:
- data.xlsx - all the data combined int one excel file
- cleandata.txt - Total data, not normalized 
- data_pca.txt - All the DE genes, selected from the clean data
- norm.data_pca.txt - All the DE genes, selected from the normalized data
- normalized_counts.txt - Total data, normalized in deseq2 median of ratios method (>= 20)
- coldata.txt - metadata on experimental setup

R scripts:
- pca.R - pca workflow

Results:
- cos2.txt
- pca_contrib.txt

Images:
- PCA plot, normalized counts (two sizes)
- Percentage of explained variances
- variables, pca - nepotrebno