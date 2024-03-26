boxplot_esc_ratio_per_xci_group <- function(human_xci=NULL) {
  chr="X"
  f1=read.table("/Users/topah/Desktop/hipsci_codes/results/top.table_final_sva_dream_ipsc_mh2_new_FALSE_TRUE.txt")
  f2=read.table("/Users/topah/Desktop/hipsci_codes/results/top.table_final_sva_dream_ipsc_mhe2_new_FALSE_TRUE.txt")
  f3=read.table("/Users/topah/Desktop/hipsci_codes/results/top.table_final_sva_dream_ipsc_ml2_new_FALSE_TRUE.txt")
  xist_ind=which(f1$gene_name=="XIST")
  f1=f1[-xist_ind,]
  f2=f2[-xist_ind,]
  f3=f3[-xist_ind,]
  f1$logFC=(-f1$logFC)
  f2$logFC=(-f2$logFC)
  f3$logFC=(-f3$logFC)

  if (chr=="X") {
    ind_chr=which(f1$chrs=="X")
  } else if (chr=="aut") {
    ind_chr=which(f1$chrs!="X" & f1$chrs!="Y" & f1$chrs!="MT")
  } else {
    ind_chr=(1:dim(f1)[1])
  }
  f1=f1[ind_chr,]
  f2=f2[ind_chr,]
  f3=f3[ind_chr,]

  source(file.path(path,"codes/get_ase_matrix.R"))
  res=get_ase_matrix(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=0)
  source(file.path(path,"codes/plot_heatmap.R"))
  res_10=plot_heatmap(res=res, min_nonna_num=10)
  gene_ids=res_10$gene_info$gene_id
  sig=res_10$sig
  D=res_10$D

  if (is.null(human_xci)) {
  ind_escape=which((f1$class[match(extract_text_before_first_dot(gene_ids),extract_text_before_first_dot(rownames(f1)))]=="Escape") &
                   f1$region[match(extract_text_before_first_dot(gene_ids),extract_text_before_first_dot(rownames(f1)))]!="PAR")

  ind_variable=which((f1$class[match(extract_text_before_first_dot(gene_ids),extract_text_before_first_dot(rownames(f1)))]=="Variable") &
    f1$region[match(extract_text_before_first_dot(gene_ids),extract_text_before_first_dot(rownames(f1)))]!="PAR")

  ind_inactive=which((f1$class[match(extract_text_before_first_dot(gene_ids),extract_text_before_first_dot(rownames(f1)))]=="Inactive") &
                     f1$region[match(extract_text_before_first_dot(gene_ids),extract_text_before_first_dot(rownames(f1)))]!="PAR")

  } else {
    ind_escape=which((human_xci[match(gene_ids,names(human_xci))]=="Escape") &
                       f1$region[match(extract_text_before_first_dot(gene_ids),extract_text_before_first_dot(rownames(f1)))]!="PAR")
    
    ind_variable=which((human_xci[match(gene_ids,names(human_xci))]=="Variable") &
                         f1$region[match(extract_text_before_first_dot(gene_ids),extract_text_before_first_dot(rownames(f1)))]!="PAR")
    
    ind_inactive=which((human_xci[match(gene_ids,names(human_xci))]=="Inactive") &
                         f1$region[match(extract_text_before_first_dot(gene_ids),extract_text_before_first_dot(rownames(f1)))]!="PAR")
    
  }
  escape=rowSums(sig[ind_escape,],na.rm=TRUE)/rowSums(!is.na(sig[ind_escape,]))
  inactive=rowSums(sig[ind_inactive,],na.rm=TRUE)/rowSums(!is.na(sig[ind_inactive,]))
  variable=rowSums(sig[ind_variable,],na.rm=TRUE)/rowSums(!is.na(sig[ind_variable,]))

  df.ie=data.frame(esc_ratio=c(escape,inactive,variable))
  df.ie$status=c(rep("Escape",length(escape)),rep("Inactive",length(inactive)),rep("Variable",length(variable)))
  df.ie$status=as.factor(df.ie$status)

  df.ie$status=factor(df.ie$status,levels=c("Escape","Variable","Inactive"))

  my_comparisons <- list(c("Escape", "Inactive"),  c("Variable","Inactive"))
  library(ggpubr)
  p.xci=ggplot(df.ie,aes(x=status,y=esc_ratio,color=status)) +
    geom_boxplot(outlier.shape = NA,size=1) +
    geom_point(position=position_jitterdodge()) +
    stat_compare_means(comparisons = my_comparisons,method="wilcox.test",method.args = list(alternative = "greater"),size=5) +
    theme_minimal() +
    scale_color_manual(name = "Known XCI status in human tissues", breaks=c(c("Escape","Variable","Inactive")),values=c('red',"darkgoldenrod1",'blue')) +
    theme(legend.position = "none") +
    xlab("Known XCI status in human tissues") +
    ylab("Fraction of cell lines\n where gene ASE > 0.1") +
    ylim(0,1.2) +
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
    theme(text = element_text(size=20))

  ggsave("/Users/topah/Desktop/hipsci_codes/figures/boxplot_esc_ratio_per_xci_group.pdf", p.xci, width=15,height=12,units="cm",limitsize = FALSE)

}