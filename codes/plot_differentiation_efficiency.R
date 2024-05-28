# -----------------------------------------------------------
# Script Name: diff_eff.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

diff_eff <- function() {
  library(readxl)
  deff = read_xlsx("1-s2.0-S2666979X2300040X-mmc2.xlsx",sheet="TableS1")
  deff=as.data.frame(deff)
  cells=deff$`TableS1: Metadata for the human iPSC lines analysed in the study. Related to Figures 1A, 2A-2B and S1A-S1B.`
  eff=deff$...20
  cells=cells[-1]
  eff=eff[-1]

  load("counts.RData")
  jj=match(cells,lines_info$lines)
  eff_df=data.frame(cells=cells,cluster=lines_info$cluster_group[jj],eff=eff)

  library(ggplot2)
  library(ggpubr)
  my_comparisons <- list( c("Male", "Group1"), c("Male", "Group2"), c("Male", "Group3"), c("Group1", "Group2"), c("Group1", "Group3"), c("Group2", "Group3"))
  eff_df$cluster=as.factor(eff_df$cluster)
  eff_df$eff=as.numeric(eff_df$eff)
  eff_df=eff_df[which(!is.na(eff_df$eff) & !is.na(eff_df$cluster)),]
  eff_df$cluster=factor(eff_df$cluster,levels=c("male","Group1","Group2","Group3"),labels=c("Male","Group1","Group2","Group3"))

  p1=ggplot(eff_df,aes(x=cluster,y=eff,color=cluster)) +
    geom_boxplot(outlier.shape = NA,size=1) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    scale_color_manual(name = "", breaks=c("Male","Group1","Group2","Group3"),values=c("#3399FF","#CC0033", "mediumorchid","#FF99CC")) +
    stat_compare_means(label="Wilcoxon, p", comparisons = my_comparisons,method="wilcox.test",paired=FALSE) +
    ylab("Differentiation efficiency for iPSC-derived endoderm") +
    xlab("Groups\n (N_males=34, N_group1=21, N_group2=15, N_group3=15)") +
    theme_minimal() +
    theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        plot.title = element_text(size=20)) +
    theme(panel.background = element_rect(fill = "white", colour = "white")) 

  ggsave(filename="diff_efficiency.pdf",p1,width=22,height=20,units="cm",limitsize = FALSE)
  dim(eff_df) # 85 lines

  table(eff_df$cluster)
  #male Group1 Group2 Group3 
  #34     21     15     15 
}
