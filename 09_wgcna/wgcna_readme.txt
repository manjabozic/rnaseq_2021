Data:
- htseq_name_counts.txt - raw data, not normalized
- normalized_counts.txt - Total data, normalized in deseq2 median of ratios method (>= 20)
- phenodata.txt - metadata on the experimental setup
- clusters_olddata.txt - lists of significant clusters for both conditions using normalized_counts.txt
- clusters_newnorm.txt - lists of significant clusters for both conditions using newly normalized
			 htseq_name_counts.txt
- module.gene.mapping_og.txt - all the genes and the modules they belong to, OG
- module.gene.mapping_newnorm.txt - all the genes and the modules they belong to, newnorm
- significant_clusters.xlsx/significant_clusters01.xlsx  - significant clusters/modules and the genes 
			 that belong to them

Scripts:
- wgcna.R

Plots:
- correlationplots.pdf
- corrplot_ogdata.pdf - correlation plot, used normalized_counts.txt
- corrplot_newnorm.pdf - correlation plot, used htseq_name_counts.txt and then normalized them in R according to
			 instructions from a YT tutorial (Bioinformagician)