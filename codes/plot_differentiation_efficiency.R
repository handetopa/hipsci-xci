# -----------------------------------------------------------
# Script Name: plot_differentiation_efficiency.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

plot_differentiation_efficiency <- function(path) {
  #library(readxl)
  #library(ggplot2)
  #library(ggpubr)
  hipsci_datadir = file.path(path,"data")
  deff = read_xlsx("data/Supplement_Puigdevall_etal_2023.xlsx",sheet="TableS1")
  deff=as.data.frame(deff)
  cells=deff$`TableS1: Metadata for the human iPSC lines analysed in the study. Related to Figures 1A, 2A-2B and S1A-S1B.`
  eff=deff$...20
  cells=cells[-1]
  eff=eff[-1]

  load(file.path(hipsci_datadir,"data_for_DE_ipsc.RData"))
  clusters=read.table(file.path(hipsci_datadir,"female_clusters.txt"),header=TRUE,sep="\t")
  D_simple$xist_group=clusters$clusterName[match(D_simple$lines,clusters$lines)]
  D_simple$xist_group[which(D_simple$sex=="male")]="male"
  
  jj=match(cells,D_simple$lines)
  eff_df=data.frame(cells=cells,cluster=D_simple$xist_group[jj],eff=eff)

  my_comparisons <- list( c("Male", "Group1"), c("Male", "Group2"), c("Male", "Group3"), c("Group1", "Group2"), c("Group1", "Group3"), c("Group2", "Group3"))
  eff_df$cluster=as.factor(eff_df$cluster)
  eff_df$eff=as.numeric(eff_df$eff)
  eff_df=eff_df[which(!is.na(eff_df$eff) & !is.na(eff_df$cluster)),]
  eff_df$cluster=factor(eff_df$cluster,levels=c("male","G1","G2","G3"),labels=c("Male","Group1","Group2","Group3"))
  n_males=sum(eff_df$cluster=="Male") #34
  n_g1=sum(eff_df$cluster=="Group1") #21
  n_g2=sum(eff_df$cluster=="Group2") #15
  n_g3=sum(eff_df$cluster=="Group3") #15
  n_all=dim(eff_df)[1] # 85 lines
  
  p1=ggplot(eff_df,aes(x=cluster,y=eff,color=cluster)) +
    geom_boxplot(outlier.shape = NA,size=1) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    scale_color_manual(name = "", breaks=c("Male","Group1","Group2","Group3"),values=c("#3399FF","#CC0033", "mediumorchid","#FF99CC")) +
    stat_compare_means(label="Wilcoxon, p", comparisons = my_comparisons,method="wilcox.test",paired=FALSE) +
    ylab("Differentiation efficiency for iPSC-derived endoderm") +
    xlab(paste("Groups\n (N_males=",n_males,", N_group1=",n_g1,", N_group2=",n_g2,", N_group3=",n_g3,")",sep="")) +
    theme_minimal() +
    theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        plot.title = element_text(size=20)) +
    theme(panel.background = element_rect(fill = "white", colour = "white")) 

  ggsave(filename="figures/diff_efficiency.pdf",p1,width=22,height=20,units="cm",limitsize = FALSE)
}
