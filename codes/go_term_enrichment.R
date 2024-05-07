# -----------------------------------------------------------
# Script Name: go_term_enrichment.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

mh=read.table("results/top.table_final_sva_dream_ipsc_mh2_new_FALSE_TRUE.txt")
mh_aut=mh[which(mh$chrs != "X" & mh$chrs!="Y" & mh$chrs!="MT"),]

mhe=read.table("results/top.table_final_sva_dream_ipsc_mhe2_new_FALSE_TRUE.txt")
mhe_aut=mhe[which(mhe$chrs != "X" & mhe$chrs!="Y" & mhe$chrs!="MT"),]

ml=read.table("results/top.table_final_sva_dream_ipsc_ml2_new_FALSE_TRUE.txt")
ml_aut=ml[which(ml$chrs != "X" & ml$chrs!="Y" & ml$chrs!="MT"),]

mh_up=which(mh_aut$adj.P.Val<0.05 & mh_aut$logFC<0)
mh_down=which(mh_aut$adj.P.Val<0.05 & mh_aut$logFC>0)
mhe_up=which(mhe_aut$adj.P.Val<0.05 & mhe_aut$logFC<0)
mhe_down=which(mhe_aut$adj.P.Val<0.05 & mhe_aut$logFC>0)
ml_up=which(ml_aut$adj.P.Val<0.05 & ml_aut$logFC<0)
ml_down=which(ml_aut$adj.P.Val<0.05 & ml_aut$logFC>0)

G3_up=setdiff(ml_up,c(ml_down,mhe_up,mhe_down,mh_down,mh_up))
G3_down=setdiff(ml_down,c(ml_up,mhe_up,mhe_down,mh_down,mh_up))
G2_up=setdiff(mhe_up,c(mhe_down,ml_up,ml_down,mh_down,mh_up))
G2_down=setdiff(mhe_down,c(mhe_up,ml_up,ml_down,mh_down,mh_up))
G3_G2_up=setdiff(intersect(mhe_up,ml_up),intersect(mhe_up,intersect(ml_up,mh_down)))
G3_G2_down=setdiff(intersect(mhe_down,ml_down),intersect(mhe_down,intersect(ml_down,mh_up)))

source("codes/extract_text_before_first_dot.R")

ind=G3_down

res=ml_aut
#BiocManager::install("clusterProfiler")
#BiocManager::install(organism, character.only = TRUE)
library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

original_gene_list <- res$logFC
names(original_gene_list)=extract_text_before_first_dot(rownames(res))
gene_list = sort(original_gene_list, decreasing = TRUE)
#sig_genes_df = subset(res, adj.P.Val < 0.05)
sig_genes_df = res[ind,]
genes <- sig_genes_df$logFC
names(genes) <- extract_text_before_first_dot(rownames(sig_genes_df))
genes <- names(genes)
#genes <- names(genes)[genes > 0]
go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = TRUE,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
onlyG3_down=go_enrich@result
onlyG3_down$Category="uniqueG3:down"
onlyG3_down=onlyG3_down[c(11,1:10)]



ind=G2_down
res=mhe_aut
original_gene_list <- res$logFC
names(original_gene_list)=extract_text_before_first_dot(rownames(res))
gene_list = sort(original_gene_list, decreasing = TRUE)
#sig_genes_df = subset(res, adj.P.Val < 0.05)
sig_genes_df = res[ind,]
genes <- sig_genes_df$logFC
names(genes) <- extract_text_before_first_dot(rownames(sig_genes_df))
genes <- names(genes)
#genes <- names(genes)[genes > 0]
go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = TRUE,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
onlyG2_down=go_enrich@result
onlyG2_down$Category="uniqueG2:down"
onlyG2_down=onlyG2_down[c(11,1:10)]


ind=G3_G2_up
res=ml_aut
original_gene_list <- res$logFC
names(original_gene_list)=extract_text_before_first_dot(rownames(res))
gene_list = sort(original_gene_list, decreasing = TRUE)
#sig_genes_df = subset(res, adj.P.Val < 0.05)
sig_genes_df = res[ind,]
genes <- sig_genes_df$logFC
names(genes) <- extract_text_before_first_dot(rownames(sig_genes_df))
genes <- names(genes)
#genes <- names(genes)[genes > 0]
go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = TRUE,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
G3_G2_up=go_enrich@result
G3_G2_up$Category="G2_and_G3:up"
G3_G2_up=G3_G2_up[c(11,1:10)]


ind=G3_G2_down
res=ml_aut
original_gene_list <- res$logFC
names(original_gene_list)=extract_text_before_first_dot(rownames(res))
gene_list = sort(original_gene_list, decreasing = TRUE)
#sig_genes_df = subset(res, adj.P.Val < 0.05)
sig_genes_df = res[ind,]
genes <- sig_genes_df$logFC
names(genes) <- extract_text_before_first_dot(rownames(sig_genes_df))
genes <- names(genes)
#genes <- names(genes)[genes > 0]
go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = TRUE,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
G3_G2_down=go_enrich@result
G3_G2_down$Category="G2_and_G3:down"
G3_G2_down=G3_G2_down[c(11,1:10)]

df_combined=rbind(onlyG3_down,onlyG2_down,G3_G2_down,G3_G2_up)
write_xlsx(df_combined,"results/go_term_enrichment_autosomal_G2_G3.xlsx")



#G3_down=go_enrich@result

# res=hl_aut
# original_gene_list <- res$logFC
# names(original_gene_list)=extract_text_before_first_dot(rownames(res))
# gene_list = sort(original_gene_list, decreasing = TRUE)
# sig_genes_df = subset(res, adj.P.Val < 0.05)
# genes <- sig_genes_df$logFC
# names(genes) <- extract_text_before_first_dot(rownames(sig_genes_df))
# genes <- names(genes)[genes < 0]
# go_enrich <- enrichGO(gene = genes,
#                       universe = names(gene_list),
#                       OrgDb = organism, 
#                       keyType = 'ENSEMBL',
#                       readable = T,
#                       ont = "ALL",
#                       pvalueCutoff = 0.05, 
#                       qvalueCutoff = 0.10)

hl=read.table("results/top.table_final_sva_dream_ipsc_mhe2_new_FALSE_TRUE.txt")
hl_aut=hl[which(hl$chrs != "X" & hl$chrs!="Y" & hl$chrs!="MT"),]
res=hl_aut
original_gene_list <- res$logFC
names(original_gene_list)=extract_text_before_first_dot(rownames(res))
gene_list = sort(original_gene_list, decreasing = TRUE)
sig_genes_df = subset(res, adj.P.Val < 0.05)
genes <- sig_genes_df$logFC
names(genes) <- extract_text_before_first_dot(rownames(sig_genes_df))
genes <- names(genes)[genes > 0]
go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism,
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)
G2_down=go_enrich@result

G3_down_common_in_G2_down=G3_down[match(intersect(G2_down$ID,G3_down$ID),G3_down$ID),]
G2_down_common_in_G3_down=G2_down[match(intersect(G2_down$ID,G3_down$ID),G2_down$ID),]

library(writexl)
write_xlsx(G3_down,"results/go_term_enrichment_autosomal_G3_down.xlsx")
write_xlsx(G2_down,"results/go_term_enrichment_autosomal_G2_down.xlsx")
write_xlsx(G3_down_common_in_G2_down,"results/go_term_enrichment_autosomal_G3_down_G2G3down.xlsx")
write_xlsx(G2_down_common_in_G3_down,"results/go_term_enrichment_autosomal_G2_down_G2G3down.xlsx")

G3_down$Category="G3:down"
G2_down$Category="G2:down"

G3_down=G3_down[c(11,1:10)]
G2_down=G2_down[c(11,1:10)]

G2G3=rbind(G3_down,G2_down)
G2G3$common_in_G2down_G3down=(G2G3$ID %in% G3_down_common_in_G2_down$ID)
write_xlsx(G2G3,"results/go_term_enrichment_autosomal_G2_down_G3_down.xlsx")






hl=read.table("results/top.table_final_sva_dream_ipsc_mh2_new_FALSE_TRUE.txt")
hl_aut=hl[which(hl$chrs != "X" & hl$chrs!="Y" & hl$chrs!="MT"),]
res=hl_aut
original_gene_list <- res$logFC
names(original_gene_list)=extract_text_before_first_dot(rownames(res))
gene_list = sort(original_gene_list, decreasing = TRUE)
sig_genes_df = subset(res, adj.P.Val < 0.05)
genes <- sig_genes_df$logFC
names(genes) <- extract_text_before_first_dot(rownames(sig_genes_df))
genes <- names(genes)[genes > 0]
go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism,
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)
G1_down=go_enrich@result
