boxplot_esc_ratio_in_female_groups_per_xci_group <- function(human_xci=NULL) {
  
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

  high_escape=rowSums(sig[ind_escape,which(D$xist_group=="High-XIST female")],na.rm=TRUE)/rowSums(!is.na(sig[ind_escape,which(D$xist_group=="High-XIST female")]))
  high_inactive=rowSums(sig[ind_inactive,which(D$xist_group=="High-XIST female")],na.rm=TRUE)/rowSums(!is.na(sig[ind_inactive,which(D$xist_group=="High-XIST female")]))
  low_escape=rowSums(sig[ind_escape,which(D$xist_group=="Low-XIST female")],na.rm=TRUE)/rowSums(!is.na(sig[ind_escape,which(D$xist_group=="Low-XIST female")]))
  low_inactive=rowSums(sig[ind_inactive,which(D$xist_group=="Low-XIST female")],na.rm=TRUE)/rowSums(!is.na(sig[ind_inactive,which(D$xist_group=="Low-XIST female")]))
  high_variable=rowSums(sig[ind_variable,which(D$xist_group=="High-XIST female")],na.rm=TRUE)/rowSums(!is.na(sig[ind_variable,which(D$xist_group=="High-XIST female")]))
  low_variable=rowSums(sig[ind_variable,which(D$xist_group=="Low-XIST female")],na.rm=TRUE)/rowSums(!is.na(sig[ind_variable,which(D$xist_group=="Low-XIST female")]))


  df.ie=data.frame(esc_ratio=c(high_escape,low_escape,high_inactive,low_inactive,high_variable,low_variable))
  df.ie$status=c(rep("Escape",length(high_escape)+length(low_escape)),rep("Inactive",length(high_inactive)+length(low_inactive)),rep("Variable",length(high_variable)+length(low_variable)))
  df.ie$xist_group=c(rep("high",length(high_escape)),rep("low",length(low_escape)),rep("high",length(high_inactive)),rep("low",length(low_inactive)),rep("high",length(high_variable)),rep("low",length(low_variable)))
  df.ie$status=as.factor(df.ie$status)
  df.ie$xist_group=as.factor(df.ie$xist_group)

  df.ie$status=factor(df.ie$status,levels=c("Escape","Variable","Inactive"))
  df.ie$xist_group=factor(df.ie$xist_group,levels=c("low","high"))

  library(ggpubr)

  df.ie2=df.ie
  df.ie2$status[which(df.ie2$status=="Variable")]="Escape"
  df.ie2$status=factor(df.ie2$status,levels=c("Escape","Inactive"))
  my_comparisons <- list(c("low", "high"))
  library(ggpubr)
  p.xci.xist=ggplot(df.ie2,aes(x=status,y=esc_ratio,color=xist_group)) +
    geom_boxplot(outlier.shape = NA,size=1) +
    geom_point(position=position_jitterdodge()) +
    stat_compare_means(aes(group = xist_group),method="wilcox.test",method.args = list(alternative = "greater",paired=TRUE),size=5) + # writes wrong p-values to the plot!
    theme_minimal() +
    scale_color_manual(name = "Female lines", breaks=c(c("low","high")),values=c("#FF99CC","#CC0033"),labels=c("Low-XIST","High-XIST")) +
    xlab("Known XCI status in human tissues") +
    ylab("Fraction of cell lines\n where gene ASE > 0.1") +
    scale_x_discrete(breaks=c("Escape","Inactive"),labels=c("Escape\n+Variable", "Inactive")) +
    ylim(0,1.2) +
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
    theme(text = element_text(size=20)) 
  
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/boxplot_esc_ratio_in_female_groups_per_xci_group_with_wrong_pvals.pdf", p.xci.xist, width=24,height=19,units="cm",limitsize = FALSE)

  AA=df.ie2$esc_ratio[(df.ie2$xist_group=="high" & (df.ie2$status=="Escape"))]
  BB=df.ie2$esc_ratio[(df.ie2$xist_group=="low" & (df.ie2$status=="Escape"))]
  print(paste("Paired one-sided Wilcoxon p-value for escape+variable genes: ",wilcox.test(BB,AA,paired=TRUE)$p.value,sep=""))
  
  AA=df.ie2$esc_ratio[(df.ie2$xist_group=="high" & (df.ie2$status=="Inactive"))]
  BB=df.ie2$esc_ratio[(df.ie2$xist_group=="low" & (df.ie2$status=="Inactive"))]
  print(paste("Paired one-sided Wilcoxon p-value for inactive genes: ",wilcox.test(BB,AA,paired=TRUE)$p.value,sep="")) #,alternative = "greater"
  
}

