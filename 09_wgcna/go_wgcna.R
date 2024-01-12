ora_treat.module01_cc <- enrichGO(gene = treat.module01$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "CC",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)
ora_treat.module02_cc <- enrichGO(gene = treat.module02$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "CC",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)

ora_treat.module05_cc <- enrichGO(gene = treat.module05$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "CC",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)

ora_treat.module03_cc <- enrichGO(gene = treat.module03$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "CC",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)
ora_treat.module04_cc <- enrichGO(gene = treat.module04$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "CC",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)

ora_treat.module06_cc <- enrichGO(gene = treat.module06$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "CC",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)


ora_treat.module01_mf <- enrichGO(gene = treat.module01$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "MF",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)
ora_treat.module02_mf <- enrichGO(gene = treat.module02$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "MF",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)

ora_treat.module05_mf <- enrichGO(gene = treat.module05$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "MF",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)

ora_treat.module03_mf <- enrichGO(gene = treat.module03$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "MF",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)
ora_treat.module04_mf <- enrichGO(gene = treat.module04$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "MF",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)

ora_treat.module06_mf <- enrichGO(gene = treat.module06$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "MF",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)

ora_tol.module01_cc <- enrichGO(gene = tol.module01$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "CC",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)
ora_tol.module02_cc <- enrichGO(gene = tol.module02$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "CC",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)

ora_tol.module05_cc <- enrichGO(gene = tol.module05$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "CC",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)

ora_tol.module03_cc <- enrichGO(gene = tol.module03$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "CC",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)
ora_tol.module04_cc <- enrichGO(gene = tol.module04$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "CC",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)



ora_tol.module01_mf <- enrichGO(gene = tol.module01$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "MF",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)
ora_tol.module02_mf <- enrichGO(gene = tol.module02$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "MF",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)

ora_tol.module05_mf <- enrichGO(gene = tol.module05$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "MF",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)

ora_tol.module03_mf <- enrichGO(gene = tol.module03$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "MF",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)
ora_tol.module04_mf <- enrichGO(gene = tol.module04$entrezgene_id, 
                                  universe = all_zm_genes_annotated$entrezgene_id, 
                                  OrgDb = orgdb,  
                                  keyType = "ENTREZID",
                                  ont = "MF",              # either "BP", "CC" or "MF",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05,
                                  readable = TRUE, 
                                  pool = FALSE)



write.table(ora_treat.module01_cc@result, "E:/Manja/wgcna/modules_go/treat.module01_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_treat.module02_cc@result, "E:/Manja/wgcna/modules_go/treat.module02_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_treat.module03_cc@result, "E:/Manja/wgcna/modules_go/treat.module03_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_treat.module04_cc@result, "E:/Manja/wgcna/modules_go/treat.module04_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_treat.module05_cc@result, "E:/Manja/wgcna/modules_go/treat.module05_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_treat.module06_cc@result, "E:/Manja/wgcna/modules_go/treat.module06_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_tol.module01_cc@result, "E:/Manja/wgcna/modules_go/tol.module01_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_tol.module02_cc@result, "E:/Manja/wgcna/modules_go/tol.module02_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_tol.module03_cc@result, "E:/Manja/wgcna/modules_go/tol.module03_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_tol.module04_cc@result, "E:/Manja/wgcna/modules_go/tol.module04_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_tol.module05_cc@result, "E:/Manja/wgcna/modules_go/tol.module05_cc.txt", sep = "\t", row.names = TRUE, col.names = NA)

write.table(ora_treat.module01_mf@result, "E:/Manja/wgcna/modules_go/treat.module01_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_treat.module02_mf@result, "E:/Manja/wgcna/modules_go/treat.module02_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_treat.module03_mf@result, "E:/Manja/wgcna/modules_go/treat.module03_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_treat.module04_mf@result, "E:/Manja/wgcna/modules_go/treat.module04_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_treat.module05_mf@result, "E:/Manja/wgcna/modules_go/treat.module05_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_treat.module06_mf@result, "E:/Manja/wgcna/modules_go/treat.module06_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_tol.module01_mf@result, "E:/Manja/wgcna/modules_go/tol.module01_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_tol.module02_mf@result, "E:/Manja/wgcna/modules_go/tol.module02_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_tol.module03_mf@result, "E:/Manja/wgcna/modules_go/tol.module03_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_tol.module04_mf@result, "E:/Manja/wgcna/modules_go/tol.module04_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(ora_tol.module05_mf@result, "E:/Manja/wgcna/modules_go/tol.module05_mf.txt", sep = "\t", row.names = TRUE, col.names = NA)


#BAR PLOT

barplot(ora_tol.module01_bp, showCategory=20, col = 'p.adjust', title = "Biological process - Module 01, associated with tolerance", label_format = 50, font.size = 20)
barplot(ora_tol.module02_bp, showCategory=20, col = 'p.adjust', title = "Biological process - Module 02, associated with tolerance", label_format = 50, font.size = 20)
barplot(ora_tol.module03_bp, showCategory=20, col = 'p.adjust', title = "Biological process - Module 03, associated with tolerance", label_format = 50, font.size = 20)
barplot(ora_tol.module04_bp, showCategory=20, col = 'p.adjust', title = "Biological process - Module 04, associated with tolerance", label_format = 50, font.size = 20)
barplot(ora_tol.module05_bp, showCategory=20, col = 'p.adjust', title = "Biological process - Module 05, associated with tolerance", label_format = 50, font.size = 20)


barplot(ora_treat.module01_bp, showCategory=20, col = 'p.adjust', title = "Biological process - Module 01, associated with the treatment", label_format = 50, font.size = 20)
barplot(ora_treat.module02_bp, showCategory=20, col = 'p.adjust', title = "Biological process - Module 02, associated with the treatment", label_format = 50, font.size = 20)
barplot(ora_treat.module03_bp, showCategory=20, col = 'p.adjust', title = "Biological process - Module 03, associated with the treatment", label_format = 50, font.size = 20)
barplot(ora_treat.module04_bp, showCategory=20, col = 'p.adjust', title = "Biological process - Module 04, associated with the treatment", label_format = 50, font.size = 20)
barplot(ora_treat.module05_bp, showCategory=20, col = 'p.adjust', title = "Biological process - Module 05, associated with the treatment", label_format = 50, font.size = 20)
barplot(ora_treat.module06_bp, showCategory=20, col = 'p.adjust', title = "Biological process - Module 06, associated with the treatment", label_format = 50, font.size = 20)


barplot(ora_tol.module01_cc, showCategory=20, col = 'p.adjust', title = "Cellular component - Module 01, associated with tolerance", label_format = 50, font.size = 20)
barplot(ora_tol.module04_cc, showCategory=20, col = 'p.adjust', title = "Cellular component - Module 04, associated with tolerance", label_format = 50, font.size = 20)

barplot(ora_treat.module01_cc, showCategory=20, col = 'p.adjust', title = "Cellular component - Module 01, associated with the treatment", label_format = 30, font.size = 15)
barplot(ora_treat.module03_cc, showCategory=20, col = 'p.adjust', title = "Cellular component - Module 03, associated with the treatment", label_format = 30, font.size = 15)
barplot(ora_treat.module04_cc, showCategory=20, col = 'p.adjust', title = "Cellular component - Module 04, associated with the treatment", label_format = 30, font.size = 15)

barplot(ora_tol.module01_mf, showCategory=20, col = 'p.adjust', title = "Molecular function - Module 01, associated with tolerance", label_format = 30, font.size = 15)
barplot(ora_tol.module04_mf, showCategory=20, col = 'p.adjust', title = "Molecular function - Module 04, associated with tolerance", label_format = 30, font.size = 15)

barplot(ora_treat.module01_mf, showCategory=20, col = 'p.adjust', title = "Molecular function - Module 01, associated with the treatment", label_format = 30, font.size = 15)
barplot(ora_treat.module03_mf, showCategory=20, col = 'p.adjust', title = "Molecular function - Module 03, associated with the treatment", label_format = 30, font.size = 15)
barplot(ora_treat.module04_mf, showCategory=20, col = 'p.adjust', title = "Molecular function - Module 04, associated with the treatment", label_format = 30, font.size = 15)



goplot(ora_tol.module01_bp)
goplot(ora_tol.module04_bp)
goplot(ora_treat.module03_bp)
goplot(ora_treat.module01_bp)
goplot(ora_treat.module04_bp)

dotplot(ora_tol.module01_bp, showCategory=20, font.size = 15) + ggtitle("Biological process - Module 01, associated with tolerance")
dotplot(ora_tol.module04_bp, showCategory=20, font.size = 15) + ggtitle("Biological process - Module 04, associated with tolerance")
dotplot(ora_treat.module01_bp, showCategory=20, font.size = 15) + ggtitle("Biological process - Module 01, associated with treatment conditions")
dotplot(ora_treat.module04_bp, showCategory=20, font.size = 15) + ggtitle("Biological process - Module 04, associated with treatment conditions")
dotplot(ora_treat.module03_bp, showCategory=20, font.size = 15) + ggtitle("Biological process - Module 03, associated with treatment conditions")

dotplot(ora_zm_cc, showCategory=20, font.size = 20) + ggtitle("Cellular component, ORA")
dotplot(ora_zm_mf, showCategory=20, font.size = 20) + ggtitle("Molecular function, ORA")


ora_analysis_bp <- pairwise_termsim(ora_zm_bp)
ora_analysis_cc <- pairwise_termsim(ora_zm_cc)
ora_analysis_mf <- pairwise_termsim(ora_zm_mf)


emapplot(ora_treat.module01_bp, color = "qvalue", cex.params = list(category_node = 1.5), layout.params = list(layout = 'circle'), 
         repel = TRUE, cex_label_group = 20) + ggtitle("Biological process, module 01, ORA")


emapplot(ora_analysis_cc, color = "qvalue", cex.params = list(category_node = 1.5), layout.params = list(layout = 'dh'), cex_label_group = 50, max.overlaps = 5)+ ggtitle("Cellular component, ORA")
emapplot(ora_analysis_mf, color = "qvalue", cex.params = list(category_node = 1.5), layout.params = list(layout = 'dh'), repel = TRUE, cex_label_group = 20)+ ggtitle("Molecular function, ORA")

# layout parameters: 'circle', 'gem', 'dh', 'graphopt', 'mds', 'fr', 'kk', 'lgl'