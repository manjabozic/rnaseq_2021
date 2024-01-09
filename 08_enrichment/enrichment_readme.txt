Data:
- cleandata.txt - Total data, not normalized 
- norm.data_pca.txt - All the DE genes, selected from the normalized data
- Zm-B73-REFERENCE-NAM.UniProt.proteins - annotation containing all info on the coding proteins
- Zm-B73-REFERENCE-NAM.UniProt.GO - annotatation containing the ensemblID and and corresponding UNIPROT and GO id
- zea_mays_uniprot_annotation.txt 
- de_susc_6h_sign - significant DEG, susceptible genotype, C vs. T, 6h
- de_susc_6h_24h_sign - significant DEG, susceptible genotype, C vs. T, 24h
- de_tol_6h_sign - significant DEG, tolerant genotype, C vs. T, 6h
- de_tol_6h_24h_sign - significant DEG, tolerant genotype, C vs. T, 24h
- all_genes - list of all the DE genes (for the universe in GO)
- species - list of all species and their IDs available for GO
- entrezid_de_susc_6h - DEGs that have available entrezid instead of ensembl,susceptible genotype, C vs. T, 6h 
- entrezid_de_susc_24h - DEGs that have available entrezid instead of ensembl,susceptible genotype, C vs. T, 24h 
- entrezid_de_tol_6h - DEGs that have available entrezid instead of ensembl,tolerant genotype, C vs. T, 6h 
- entrezid_de_tol_24h - DEGs that have available entrezid instead of ensembl,tolerant genotype, C vs. T, 24h
- de_genes.txt - Selected 241 DE genes, from the normalized data
- de_sign.xlsx - excel file containing most of the previous information
- info.txt - metadata on the experimental setup
- diff_genes_swissprot.txt - ?
- diff_genes_swissprot_all.txt - ?


R scripts:
- funct_enrichment.R - clustering workflow

Subdirectories:
results/