# -----------------------------------------------------------
# Script Name: mean_ase_by_groups.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

mean_ase_by_groups <- function() {

  path = "/Users/topah/Desktop/hipsci_codes"
  source(file.path(path,"codes/get_ase_matrix.R"))
  res=get_ase_matrix(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=0)
  
  D=res$D
  ase_ratios=res$ase_ratios
  clusters=read.table("/Users/topah/Desktop/hipsci_codes/data/female_clusters.txt",header=TRUE,sep="\t")
  D$xist_group=clusters$clusterName[match(D$lines,clusters$lines)]
  D$xist_group[which(D$sex=="male")]="male"
  
  df=data.frame(mean_ase=c(colMeans(ase_ratios[,which(D$xist_group=="high")],na.rm=TRUE),
                           colMeans(ase_ratios[,which(D$xist_group=="highescape")],na.rm=TRUE),
                           colMeans(ase_ratios[,which(D$xist_group=="low")],na.rm=TRUE)),
                cluster=c(rep("high",sum(D$xist_group=="high")),
                          rep("highescape",sum(D$xist_group=="highescape")),
                          rep("low",sum(D$xist_group=="low"))))
  df$cluster=factor(df$cluster,levels=c("high","highescape","low"),labels=c("group 1","group 2","group 3"))
  
  my_comparisons <- list(c("group 1", "group 2"),  c("group 2", "group 3"), c("group 1", "group 3"))
  library(ggpubr)
  pp=ggplot(df,aes(x=cluster,y=mean_ase,color=cluster)) +
    geom_boxplot(outlier.shape = NA,size=1) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    stat_compare_means(comparisons = my_comparisons,method="wilcox.test", size = 5) +
    theme_minimal() +
    scale_color_manual(name = "", breaks=c("group 1","group 2", "group 3"),values=c("#CC0033", "mediumorchid","#FF99CC")) +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Mean ASE") +
    scale_y_continuous(breaks=c(0,0.25,0.5)) +
    theme(text = element_text(size=20))
  p1=wilcox.test(colMeans(ase_ratios[,which(D$xist_group=="high")],na.rm=TRUE),colMeans(ase_ratios[,which(D$xist_group=="highescape")],na.rm=TRUE))$p.value
  p2=wilcox.test(colMeans(ase_ratios[,which(D$xist_group=="highescape")],na.rm=TRUE),colMeans(ase_ratios[,which(D$xist_group=="low")],na.rm=TRUE))$p.value
  p3=wilcox.test(colMeans(ase_ratios[,which(D$xist_group=="high")],na.rm=TRUE),colMeans(ase_ratios[,which(D$xist_group=="low")],na.rm=TRUE))$p.value
  print(paste("Wilcoxon test p-val (G1-G2): ",p1,sep=""))
  print(paste("Wilcoxon test p-val (G2-G3): ",p2,sep=""))
  print(paste("Wilcoxon test p-val (G1-G3): ",p3,sep=""))
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/mean_ase_by_groups.pdf", pp, width=15,height=15,units="cm",limitsize = FALSE)

}