# -----------------------------------------------------------
# Script Name: esc_ratio_by_xist_groups.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

esc_ratio_by_xist_groups <- function() {
  
  source("codes/get_ase_matrix.R")
  res=get_ase_matrix(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=0)
  
  D=res$D
  sig=res$sig
  df=data.frame(esc_ratio=c(colSums(sig[,which(D$xist_group=="High-XIST female")],na.rm=TRUE) / colSums(!is.na(sig[,which(D$xist_group=="High-XIST female")])),
                            colSums(sig[,which(D$xist_group=="Low-XIST female")],na.rm=TRUE) / colSums(!is.na(sig[,which(D$xist_group=="Low-XIST female")]))),
                cluster=c(rep("High-XIST",sum(D$xist_group=="High-XIST female")),
                          rep("Low-XIST",sum(D$xist_group=="Low-XIST female"))))
  df$cluster=factor(df$cluster,levels=c("High-XIST","Low-XIST"))
  
  my_comparisons <- list(c("High-XIST","Low-XIST"))
  library(ggpubr)
  pp=ggplot(df,aes(x=cluster,y=esc_ratio,color=cluster)) +
    geom_boxplot(outlier.shape = NA,size=1) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    stat_compare_means(comparisons = my_comparisons, method="wilcox.test", size = 6) +
    theme_minimal() +
    scale_color_manual(name = "", breaks=c("High-XIST","Low-XIST"),values=c("#CC0033","#FF99CC")) +
    theme(legend.position = "none") +
    xlab("Female lines") +
    ylab("Fraction of genes with\n ASE > 0.1") +
    ylim(0,1.1) +
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
    theme(text = element_text(size=22))
  p1=t.test(colSums(sig[,which(D$xist_group=="High-XIST female")],na.rm=TRUE) / colSums(!is.na(sig[,which(D$xist_group=="High-XIST female")])),
            colSums(sig[,which(D$xist_group=="Low-XIST female")],na.rm=TRUE) / colSums(!is.na(sig[,which(D$xist_group=="Low-XIST female")])))$p.value
  p2=wilcox.test(colSums(sig[,which(D$xist_group=="High-XIST female")],na.rm=TRUE) / colSums(!is.na(sig[,which(D$xist_group=="High-XIST female")])),
            colSums(sig[,which(D$xist_group=="Low-XIST female")],na.rm=TRUE) / colSums(!is.na(sig[,which(D$xist_group=="Low-XIST female")])))$p.value
  print(paste("t-test p-val: ",p1,sep=""))
  print(paste("Wilcoxon test p-val: ",p2,sep=""))
  ggsave("figures/esc_ratio_by_xist_groups.pdf", pp, width=15,height=15,units="cm",limitsize = FALSE)
  
  t.test(colSums(sig[,which(D$xist_group=="High-XIST female")],na.rm=TRUE)/colSums(!is.na(sig[,which(D$xist_group=="High-XIST female")])))$estimate
  t.test(colSums(sig[,which(D$xist_group=="High-XIST female")],na.rm=TRUE)/colSums(!is.na(sig[,which(D$xist_group=="High-XIST female")])))$conf.int
  t.test(colSums(sig[,which(D$xist_group=="Low-XIST female")],na.rm=TRUE)/colSums(!is.na(sig[,which(D$xist_group=="Low-XIST female")])))$estimate
  t.test(colSums(sig[,which(D$xist_group=="Low-XIST female")],na.rm=TRUE)/colSums(!is.na(sig[,which(D$xist_group=="Low-XIST female")])))$conf.int
  
}