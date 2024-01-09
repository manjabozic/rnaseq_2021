#########################################################################################################################
#################################### FUNCTIONAL ENRICHMENT ANALYSIS #####################################################
#########################################################################################################################

# A functional enrichement analysis will determine whether some functions are enriched in your set of differentially 
# expressed genes.One important goal is to gain a higher view and not only deal with individual genes but understand
# which pathways are involved in the response.

# Once we obtain a list of genes, we have multiple analysis to perform to go beyond a simple list of genes:
# Annotating our list of genes with cross-databases identifiers and descriptions (Entrezid, Uniprot, KEGG, etc.).
# Performing Over-Representation Analysis (ORA) or Gene Set Enrichment Analysis (GSEA) using R or webtools.
# Interpreting the results 

setwd("E:/Manja/go") #set the working directory

genes <- read.delim("E:/Manja/go/de_genes.txt", header=TRUE) # do not use gene id as row.names
genes <- read.delim("E:/Manja/go/norm.data_pca.txt", header=TRUE)

library("biomartr")
library("clusterProfiler")
library("tidyverse")
library("enrichplot")
library("biomaRt")
library("AnnotationHub")

ah <- AnnotationHub() # all these steps are to load an org.db file needed for go annotation analysis
query(ah, "OrgDb")
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Zea mays")[[1]]

################################################################################################################
#Over Representation Analysis (ORA)
################################################################################################################

# Over Representation Analysis is searching for biological functions or pathways that are enriched in a list obtained through 
# experimental studies compared to the complete list of functions/pathways. Usually, ORA makes use of gene ontologies (GO) 
# where each gene receives one or multiple layers of information on their function, cellular localization, etc.

# ORA analysis rely on a mathematical equation to compute a p-value for a given gene set classified under a certain GO.
# Elements included n this equation are:
# N - the total number of genes in the background distribution (the “universe” of our transcriptome); all genes with annotations
# M - the number of genes within that distribution that are annotated (directly or indirectly) to the gene set of interest.
# n - the size of the list of genes of interest (the size of your “drawing”).
# k - the number of genes “drawn” within that list which are annotated to the gene set.


# The Gene Ontology (GO) produces a bird’s-eye view of biological systems by building a tree of terms related to biological 
# functions. This is particularly helpful when dealing with results from genome-wide experiments (e.g. transcriptomics) 
# since classifying genes into groups of related functions can assist in the interpretation of results. 
# Rather than focusing on each gene, one by one, the researcher gets access to metabolic pathways, functions related to 
# development, etc.

# AmiGO2 (https://amigo.geneontology.org/amigo/landing) - GO data base

# AgriGO (http://systemsbiology.cau.edu.cn/agriGOv2/index.php) - a webtool accessible here to perform gene ontology analyses.

# Panther (https://www.pantherdb.org/) - a webtool accessible here to perform gene ontology analyses.


#Annotation and ORA in R
##############################################################################################################

# Annotating your DE genes with Ensembl and biomartr

#Drost HG, Paszkowski J. Biomartr: genomic data retrieval with R. Bioinformatics (2017) 33(8): 1216-1217. doi:10.1093/bioinformatics/btw821.

# If biomart refuses to query Ensembl again, run this command

biomaRt::biomartCacheClear() # to solve a known bug https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/335

#load zea mays info from biomart

biomartr::organismBM(organism = "Zea mays") # retrieval of genomes and their annotation from ENSEMBL

zm_attributes <- biomartr::organismAttributes("Zea mays") %>% #retrieve all the available information contained in the dataset
                 filter(dataset == "zmays_eg_gene")

attributes_to_retrieve = c("ensembl_transcript_id", "ensembl_peptide_id","gene_biotype",
                            "entrezgene_id", "go_id", 	
                           "uniprotswissprot")

attributes_to_retrieve_go = c("go_id",
                           "name_1006", "definition_1006", "go_linkage_type",
                           "namespace_1003", "goslim_goa_accession", "goslim_goa_description",
                           "plant_reactome_pathway", "plant_reactome_reaction")

#create a file needed for go analysis, consisting of genes of interest

annotation_de <- biomartr::biomart( genes      = genes$Gene_ID,                  # genes were retrieved using biomartr::getGenome()
                                mart       = "plants_mart",                     # marts were selected with biomartr::getMarts()
                                dataset    = "zmays_eg_gene",               # datasets were selected with biomartr::getDatasets()
                                attributes = attributes_to_retrieve,            # attributes were selected with biomartr::getAttributes()
                                filters =   "ensembl_gene_id" ) # query key

annotation_de_go <- biomartr::biomart( genes      = genes$Gene_ID,                  
                                    mart       = "plants_mart",                     
                                    dataset    = "zmays_eg_gene",              
                                    attributes = attributes_to_retrieve_go,        
                                    filters =   "ensembl_gene_id" )

# Building the background/universe!

all_zm_genes <- read.delim("cleandata.txt", header = T, stringsAsFactors = F)[,1] # Directly selects the gene column

# Select the attributes to retrieve we want the correspondence of TAIR/Ensembl symbols with NCBI Entrez gene ids

all_zm_genes_annotated <- biomartr::biomart(genes = all_zm_genes, # Query the Ensembl API
                                            mart       = "plants_mart",                 
                                            dataset    = "zmays_eg_gene",           
                                            attributes = attributes_to_retrieve,        
                                            filters =  "ensembl_gene_id" )  

# for compatibility with enrichGO universe genes in the universe need to be characters and not integers (Entrez gene id)
all_zm_genes_annotated$entrezgene_id = as.character(
  all_zm_genes_annotated$entrezgene_id)

annotation_de$entrezgene_id = as.character(
  annotation_de$entrezgene_id)

write.table(all_zm_genes, "E:/Manja/go/all_genes.txt", sep = "\t", col.names = NA)
################################################################################################################################
# clusterProfiler
################################################################################################################################

# In clusterProfiler, the groupGO() function is designed for gene classification based on GO distribution at a specific level. 

go_annotation_bp <- groupGO(gene     = annotation_de$entrezgene_id,
               OrgDb    = orgdb,
               ont      = "BP",
               keyType = "ENTREZID",
               level    = 3)

go_annotation_cc <- groupGO(gene     = annotation_de$entrezgene_id,
               OrgDb    = orgdb,
               ont      = "CC",
               keyType = "ENTREZID",
               level    = 3)

go_annotation_mf <- groupGO(gene     = annotation_de$entrezgene_id,
               OrgDb    = orgdb,
               ont      = "MF",
               keyType = "ENTREZID",
               level    = 3)

write.table(go_annotation_bp, "E:/Manja/go/results/go_annotation_bp.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(go_annotation_cc, "E:/Manja/go/results/go_annotation_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(go_annotation_mf, "E:/Manja/go/results/go_annotation_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)


# Performing the ORA for Gene Ontology Biological Process class

ora_zm_bp <- enrichGO(gene = annotation_de$entrezgene_id, 
                            universe = all_zm_genes_annotated$entrezgene_id, 
                            OrgDb = orgdb,  
                            keyType = "ENTREZID",
                            ont = "BP",              # either "BP", "CC" or "MF",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = TRUE, 
                            pool = FALSE)

ora_zm_cc <- enrichGO(gene = annotation_de$entrezgene_id, 
                      universe = all_zm_genes_annotated$entrezgene_id, 
                      OrgDb = orgdb,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                      keyType = "ENTREZID",
                      ont = "CC",              # either "BP", "CC" or "MF",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE, 
                      pool = FALSE)

ora_zm_mf <- enrichGO(gene = annotation_de$entrezgene_id, 
                      universe = all_zm_genes_annotated$entrezgene_id, 
                      OrgDb = orgdb,  
                      keyType = "ENTREZID",
                      ont = "MF",              # either "BP", "CC" or "MF",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE, 
                      pool = FALSE)

write.table(ora_zm_bp@result, "E:/Manja/go/results/goresults_bp_de.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_zm_cc@result, "E:/Manja/go/results/goresults_cc_de.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_zm_mf@result, "E:/Manja/go/results/goresults_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)

# The Gene Ontology classification is very redundant meaning that parental terms overlap a lot with their related child terms. 
# The clusterProfiler package comes with a dedicated function called simplify to solve this issue.

ora_zm_bp_simplified <- clusterProfiler::simplify(ora_zm_bp) 
ora_zm_cc_simplified <- clusterProfiler::simplify(ora_zm_cc) 
ora_zm_mf_simplified <- clusterProfiler::simplify(ora_zm_mf) 

write.table(ora_zm_bp_simplified@result, "E:/Manja/go/results/go_results_bp.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_zm_cc_simplified@result, "E:/Manja/go/results/go_results_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_zm_mf_simplified@result, "E:/Manja/go/results/go_results_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)

# GO Gene Set Enrichment Analysis
# The clusterProfiler package provides the gseGO() function for gene set enrichment analysis using gene ontology. GSEA analysis 
# requires a ranked gene list, which contains three features:
#       - numeric vector: fold change or other type of numerical variable
#       - named vector: every number has a name, the corresponding gene ID
#       - sorted vector: number should be sorted in decreasing order

# The file should contain two columns, one for gene ID (no duplicated ID allowed) and another one for fold change. 
# Now this is a bit complicated
# You can prepare your own geneList via the following command:

#First, load the file where the fisrt column is Gene Ensembl ID, second is fold change

de_susc_6h = read.delim("E:/Manja/go/de_susc_6_sign.txt", header=TRUE)  # assume 1st column is ID, 2nd column is FC
de_susc_24h = read.delim("E:/Manja/go/de_susc_24_sign.txt", header=TRUE)
de_tol_6h = read.delim("E:/Manja/go/de_tol_6_sign.txt", header=TRUE) 
de_tol_24h = read.delim("E:/Manja/go/de_tol_24_sign.txt", header=TRUE)

# then you need to switch from the Gene Ensembl ID to Entrez ID

attr = c("entrezgene_id")

de_susc_6h = biomartr::biomart( genes      = de_susc_6h$Gene,                 
                                mart       = "plants_mart",                    
                                dataset    = "zmays_eg_gene",              
                                attributes = attr,            
                                filters =   "ensembl_gene_id" )

de_susc_24h = biomartr::biomart( genes      = de_susc_24h$Gene,                 
                                mart       = "plants_mart",                    
                                dataset    = "zmays_eg_gene",              
                                attributes = attr,            
                                filters =   "ensembl_gene_id" )

de_tol_6h = biomartr::biomart( genes      = de_tol_6h$Gene,                 
                                mart       = "plants_mart",                    
                                dataset    = "zmays_eg_gene",              
                                attributes = attr,            
                                filters =   "ensembl_gene_id" )

de_tol_24h = biomartr::biomart( genes      = de_tol_24h$Gene,                 
                                 mart       = "plants_mart",                    
                                 dataset    = "zmays_eg_gene",              
                                 attributes = attr,            
                                 filters =   "ensembl_gene_id" )

#this deletes the fold change, so you need to write.table with this and do it by hand in excell, so now the first column is 
# Entrez ID and the second is fold change.

de_susc_6h = read.delim("E:/Manja/go/entrezid_de_susc_6h.txt", header=TRUE)  # assume 1st column is ID, 2nd column is FC
de_susc_24h = read.delim("E:/Manja/go/entrezid_de_susc_24h.txt", header=TRUE)
de_tol_6h = read.delim("E:/Manja/go/entrezid_de_tol_6h.txt", header=TRUE) 
de_tol_24h = read.delim("E:/Manja/go/entrezid_de_tol_24h.txt", header=TRUE)

geneList_susc6h = de_susc_6h[,2] # feature 1: numeric vector
geneList_susc24h = de_susc_24h[,2] 
geneList_tol6h = de_tol_6h[,2] 
geneList_tol24h = de_tol_24h[,2] 

names(geneList_susc6h) = as.character(de_susc_6h[,1]) # feature 2: named vector
names(geneList_susc24h) = as.character(de_susc_24h[,1])
names(geneList_tol6h) = as.character(de_tol_6h[,1]) 
names(geneList_tol24h) = as.character(de_tol_24h[,1])

geneList_susc6h = sort(geneList_susc6h, decreasing = TRUE) # feature 3: decreasing orde
geneList_susc24h = sort(geneList_susc24h, decreasing = TRUE)
geneList_tol6h = sort(geneList_tol6h, decreasing = TRUE) 
geneList_tol24h = sort(geneList_tol24h, decreasing = TRUE)

gsea_s6_bp <- gseGO(geneList  = geneList_susc6h,
              OrgDb        = orgdb,
              ont          = "BP",
              minGSSize    = 25,
              maxGSSize    = 50,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

gsea_s6_cc <- gseGO(geneList  = geneList_susc6h,
                    OrgDb        = orgdb,
                    ont          = "CC",
                    minGSSize    = 25,
                    maxGSSize    = 50,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)

gsea_s6_mf <- gseGO(geneList  = geneList_susc6h,
                    OrgDb        = orgdb,
                    ont          = "MF",
                    minGSSize    = 25,
                    maxGSSize    = 50,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)

gsea_s24_bp <- gseGO(geneList  = geneList_susc24h,
                    OrgDb        = orgdb,
                    ont          = "BP",
                    minGSSize    = 25,
                    maxGSSize    = 50,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)

gsea_s24_cc <- gseGO(geneList  = geneList_susc24h,
                    OrgDb        = orgdb,
                    ont          = "CC",
                    minGSSize    = 25,
                    maxGSSize    = 50,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)

gsea_s24_mf <- gseGO(geneList  = geneList_susc24h,
                    OrgDb        = orgdb,
                    ont          = "MF",
                    minGSSize    = 25,
                    maxGSSize    = 50,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)

# GO analysis for non-model organisms

# Both the enrichGO() and gseGO() functions require an OrgDb object as the background annotation. For organisms that don’t have 
# OrgDb provided by Bioconductor, users can query one (if available) online via AnnotationHub. If there is no OrgDb available, 
# users can obtain GO annotation from other sources, e.g. from biomaRt or Blast2GO. Then the enricher() or GSEA() functions can 
# be used to perform GO analysis for these organisms, similar to the examples using wikiPathways and MSigDB. Another solution is 
# to create an OrgDb on your own using AnnotationForge package.

# The clusterProfiler package provides enricher() function for hypergeometric test and GSEA() function for gene set enrichment 
# analysis that are designed to accept user defined annotation. 

# They accept two additional parameters TERM2GENE and TERM2NAME. As indicated in the parameter names, TERM2GENE is a data.frame 
# with first column of term ID and second column of corresponding mapped gene and TERM2NAME is a data.frame with first column of 
# term ID and second column of corresponding term name. TERM2NAME is optional.

# input data for the analysis is the same as before

genes <- read.delim("E:/Manja/go/de_genes.txt", header=TRUE) #you just need the first column with the gene ids

#or you can do create a gene list object like in the previous step (entrezid transformation not necessary)

#annotation data - one column gene_ide, and the other the term we want, Go or UNIPROT

zm_annotation_data <- read.delim("E:/Manja/go/zea_mays_uniprot_annotation.txt", header=TRUE) #this needs to have two columns

x <- enricher(genes$Gene_ID, TERM2GENE = zm_annotation_data)

head(x)

write.table(x@result, "E:/Manja/go/results/x.txt", sep = "\t", row.names = TRUE, col.names = NA)
dotplot(x)
################################################################################################################################
# Visualization
################################################################################################################################

#BAR PLOT

barplot(ora_zm_cc, showCategory=20, col = 'p.adjust', title = "Cellular component", label_format = 50, font.size = 20)
barplot(ora_zm_bp, showCategory=20, col = 'p.adjust', title = "Biological process", label_format = 50, font.size = 15)
barplot(ora_zm_mf, showCategory=20, col = 'p.adjust', title = "Molecular function", label_format = 50, font.size = 20)

goplot(ora_zm_cc)
goplot(ora_zm_bp)
goplot(ora_zm_mf)

dotplot(ora_zm_bp, showCategory=20, font.size = 20) + ggtitle("Biological process, ORA")
dotplot(ora_zm_cc, showCategory=20, font.size = 20) + ggtitle("Cellular component, ORA")
dotplot(ora_zm_mf, showCategory=20, font.size = 20) + ggtitle("Molecular function, ORA")


ora_analysis_bp <- pairwise_termsim(ora_zm_bp)
ora_analysis_cc <- pairwise_termsim(ora_zm_cc)
ora_analysis_mf <- pairwise_termsim(ora_zm_mf)

emapplot(ora_analysis_bp, color = "qvalue", cex.params = list(category_node = 1.5), layout.params = list(layout = 'dh'), 
         repel = TRUE, cex_label_group = 20) + ggtitle("Biological process, ORA")
emapplot(ora_analysis_cc, color = "qvalue", cex.params = list(category_node = 1.5), layout.params = list(layout = 'dh'), cex_label_group = 50, max.overlaps = 5)+ ggtitle("Cellular component, ORA")
emapplot(ora_analysis_mf, color = "qvalue", cex.params = list(category_node = 1.5), layout.params = list(layout = 'dh'), repel = TRUE, cex_label_group = 20)+ ggtitle("Molecular function, ORA")

# layout parameters: 'circle', 'gem', 'dh', 'graphopt', 'mds', 'fr', 'kk', 'lgl'


################################################################################################################
#KEGG
################################################################################################################

# To see if your organism is referenced in the KEGG database, you can search this page: https://www.genome.jp/kegg/catalog/org_list.html
# You can also do this programmatically using R and the clusterProfiler package.

zma <- search_kegg_organism('zma', by='kegg_code')
search_kegg_organism('Arabidopsis thaliana', by='scientific_name')
genelist <- sort(annotation_de$entrezgene_id, decreasing = TRUE)
genelist <- unique(genelist)
# do the analysis

ora_analysis_kegg <- enrichKEGG(gene = annotation_de$entrezgene_id,
                                universe = all_zm_genes_annotated$entrezgene_id,
                                organism = 'zma',
                                keyType = "ncbi-geneid",
                                minGSSize = 10,
                                maxGSSize = 750,
                                pAdjustMethod = "BH",
                                qvalueCutoff = 1,
                                use_internal_data = FALSE)
dotplot(ora_analysis_kegg, 
        color = "qvalue", 
        showCategory = 10, 
        size = "Count")


ora_analysis_kegg_modules <- enrichMKEGG(gene = annotation_de$entrezgene_id,
                                         universe = all_zm_genes_annotated$entrezgene_id,
                                         organism = 'zma',
                                         keyType = "ncbi-geneid",
                                         minGSSize = 10,           # minimal size of genes annotated by Ontology term for testing.
                                         maxGSSize = 750,          # maximal size of genes annotated for testing
                                         pAdjustMethod = "BH",
                                         qvalueCutoff = 1)
dotplot(ora_analysis_kegg_modules, 
        color = "qvalue", 
        showCategory = 10, 
        size = "Count")

################################################################################################################
# Transcriptomic and metabolomic data integration
################################################################################################################

# Create a list of all DE genes containing Swissprot IDs (UniProt IDs) and use that in the webtool

annotation_de %>% 
  filter(uniprotswissprot != "") %>%                                       # to remove genes with no matching Uniprot entries
  unique() %>% 
  mutate(id_for_ipath = paste("UNIPROT",uniprotswissprot,sep = ":")) %>%   # to create an ID that iPath can use
  dplyr::select(id_for_ipath) %>%  # we keep only the relevant ID for further copy-pasting 
  unique() %>%
  write.table(., 
              file = "diff_genes_swissprot.tsv", 
              row.names = FALSE, 
              quote = FALSE)

# https://pathways.embl.de/ - Interactive Pathways Explorer v3


