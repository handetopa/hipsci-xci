violin_abslogfc_femalebiased_in_G1G2G3 <- function() {
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
  ind_de_all=which(f1$adj.P.Val<0.05 & f2$adj.P.Val<0.05 & f3$adj.P.Val<0.05 & f1$logFC>0 & f2$logFC>0 & f3$logFC>0) # For female-biased genes, 115 genes excluding XIST
  ##ind_de_all=which(f1$adj.P.Val<0.05 & f1$logFC>0) # For female-biased genes in Group 1, 123 genes excluding XIST
  #ind_de_all=which(f1$adj.P.Val<0.05 & f1$logFC<0) # For male-biased genes in G1
  df.effects=data.frame(abs_lfc=c(abs(f1$logFC[ind_de_all]), abs(f2$logFC[ind_de_all]), abs(f3$logFC[ind_de_all])),
                        groups=c(rep("Group 1",length(ind_de_all)),rep("Group 2",length(ind_de_all)),rep("Group 3",length(ind_de_all))))
  df.effects$groups=factor(df.effects$groups,levels=c("Group 1","Group 2","Group 3"))
  my_comparisons <- list( c("Group 1", "Group 2"),  c("Group 2", "Group 3") , c("Group 1", "Group 3"))
  library(ggpubr)
  p.effects=ggplot(df.effects,aes(x=groups,y=abs_lfc,color=groups)) +
    geom_violin(size=1) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    stat_compare_means(comparisons = my_comparisons,paired=TRUE,method="wilcox.test") +
    theme_minimal() +
    scale_color_manual(name = "Female\ncell lines", breaks=c("Group 1","Group 2","Group 3"),values=c("#CC0033", "mediumorchid","#FF99CC")) +
    theme(legend.position = "none") +
    xlab("Female cell lines") +
    #stat_compare_means() +
    # stat_compare_means(ref.group = "Group 1", 
    #                    method = "wilcox.test",
    #                    method.args = list(alternative = "greater")) +
    theme(text = element_text(size=20)) 
  if (chr=="X") {
    p.effects = p.effects + ylab(expression("|Log"[2]*"(Fold change)|")) 
    ggsave("/Users/topah/Desktop/hipsci_codes/figures/female_biased_inG1G2G3_abslogfc_chrX_violin.pdf", p.effects, width=15,height=15,units="cm",limitsize = FALSE)
    #ggsave("/Users/topah/Desktop/hipsci_codes/figures/male_biased_inG1_abslogfc_chrX_violin.pdf", p.effects, width=15,height=15,units="cm",limitsize = FALSE) # for male-biased
  }  else if (chr=="aut") {
    p.effects = p.effects + ylab("Absolute logFC in autosomal genes") 
    ggsave("/Users/topah/Desktop/hipsci_codes/figures/abslogfc_aut.pdf", p.effects, width=15,height=15,units="cm",limitsize = FALSE)
  }
  
  print(paste("Number of genes: ",length(ind_de_all),sep=""))
  wilcox.test(df.effects$abs_lfc[which(df.effects$groups=="Group 2")],df.effects$abs_lfc[which(df.effects$groups=="Group 3")],paired=TRUE)$p.value
  #[1] 9.245681e-18
  wilcox.test(df.effects$abs_lfc[which(df.effects$groups=="Group 1")],df.effects$abs_lfc[which(df.effects$groups=="Group 2")],paired=TRUE)$p.value
  #[1] 1.401934e-20
  wilcox.test(df.effects$abs_lfc[which(df.effects$groups=="Group 1")],df.effects$abs_lfc[which(df.effects$groups=="Group 3")],paired=TRUE)$p.value
  #[1] 1.330216e-20
  
}