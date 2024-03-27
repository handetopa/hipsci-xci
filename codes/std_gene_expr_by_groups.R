# -----------------------------------------------------------
# Script Name: std_gene_expr_by_groups.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

std_gene_expr_by_groups <- function() {
  
  path = "/Users/topah/Desktop/hipsci_codes"
  source(file.path(path,"codes/get_ase_matrix.R"))
  res=get_ase_matrix(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=0)
  source(file.path(path,"codes/plot_heatmap.R"))
  res_10=plot_heatmap(res=res, min_nonna_num=10)
  D=res_10$D
  source("/Users/topah/Desktop/hipsci_codes/codes/extract_text_before_first_dot.R")
  load("/Users/topah/Desktop/hipsci_codes/data/data_for_DE_new.RData")
  
  clusters=read.table("/Users/topah/Desktop/hipsci_codes/data/female_clusters.txt",header=TRUE,sep="\t")
  D$xist_group=clusters$clusterName[match(D$lines,clusters$lines)]
  D$xist_group[which(D$sex=="male")]="male"
  
  cpm_log_X=cpm_log[match(extract_text_before_first_dot(res_10$gene_info$gene_id),extract_text_before_first_dot(rownames(cpm_log))),match(D$lines,sub('_[^_]*$', '', colnames(cpm_log)))]
  
  std_cpm_log_X=t(scale(t(cpm_log_X)))
  xist_ind=which(extract_text_before_first_dot(rownames(cpm_log_X))=="ENSG00000229807")
  D_simple=D
  
  my_comparisons <- list(c("G1", "G2"),  c("G1","G3"), c("G2","G3"))
  df_stdexpr=data.frame(expr=c(rowMeans(std_cpm_log_X[-xist_ind,which(D_simple$xist_group=="high")]),
                             rowMeans(std_cpm_log_X[-xist_ind,which(D_simple$xist_group=="highescape")]),
                             rowMeans(std_cpm_log_X[-xist_ind,which(D_simple$xist_group=="low")])),
                      xist_group=c(rep("G1",138),rep("G2",138),rep("G3",138)))
  df_stdexpr$xist_group=factor(df_stdexpr$xist_group,levels=c("G1","G2","G3"))
  
  p.expr=ggplot(df_stdexpr,aes(x=xist_group,y=expr,color=xist_group)) +
    geom_boxplot(outlier.shape = NA,size=1) +
    geom_point(position=position_jitterdodge()) +
    stat_compare_means(comparisons = my_comparisons,method="t.test",paired=TRUE) +
    theme_minimal() +
    scale_color_manual(name = "Female group", breaks=c("G1","G2","G3"), values=c("#CC0033", "mediumorchid","#FF99CC")) +
    theme(legend.position = "none") +
    xlab("Female groups") +
    ylab("Mean standardized gene expression") +
    scale_x_discrete(breaks=c("G1","G2","G3"),labels=c("Group 1","Group 2","Group 3")) +
    theme(text = element_text(size=20)) 
  
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/std_gene_expr_by_groups.pdf", p.expr, width=15,height=15,units="cm",limitsize = FALSE)

  # Paired t.test
  t.test(df_stdexpr$expr[which(df_stdexpr$xist_group=="G1")],df_stdexpr$expr[which(df_stdexpr$xist_group=="G2")],paired=TRUE)$p.value
  #[1] 1.672785e-21  # p-value for G1-G2
  t.test(df_stdexpr$expr[which(df_stdexpr$xist_group=="G1")],df_stdexpr$expr[which(df_stdexpr$xist_group=="G3")],paired=TRUE)$p.value
  #[1] 6.760063e-26  # p-value for G1-G3
  t.test(df_stdexpr$expr[which(df_stdexpr$xist_group=="G2")],df_stdexpr$expr[which(df_stdexpr$xist_group=="G3")],paired=TRUE)$p.value
  #[1] 3.403488e-16  # p-value for G2-G3
  
  # Paired wilcox.test
  wilcox.test(df_stdexpr$expr[which(df_stdexpr$xist_group=="G1")],df_stdexpr$expr[which(df_stdexpr$xist_group=="G2")],paired=TRUE)$p.value
  # [1] 3.496911e-17  # p-value for G1-G2
  wilcox.test(df_stdexpr$expr[which(df_stdexpr$xist_group=="G1")],df_stdexpr$expr[which(df_stdexpr$xist_group=="G3")],paired=TRUE)$p.value
  # [1] 5.336343e-20  # p-value for G1-G3
  wilcox.test(df_stdexpr$expr[which(df_stdexpr$xist_group=="G2")],df_stdexpr$expr[which(df_stdexpr$xist_group=="G3")],paired=TRUE)$p.value
  # [1] 1.458278e-13  # p-value for G2-G3

}
