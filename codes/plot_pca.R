# -----------------------------------------------------------
# Script Name: plot_pca.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

plot_pca <- function(chr = "X", path) {
  
  load(file.path(path,"data/data_for_DE_new.RData"))
  clusters=read.table(file.path(path,"data/female_clusters.txt"),header=TRUE,sep="\t")
  D_simple$cluster_group=clusters$clusterName[match(D_simple$lines,clusters$lines)]
  D_simple$cluster_group[which(D_simple$sex == "male")] = "male"
  if (chr == "X") {
    i=which(gene_info$chr=="X")
    plotname =file.path(path,"figures/PCA_chrX.pdf")
  }
  else {
    i=which(gene_info$chr!="X" & gene_info$chr!="Y" & gene_info$chr!="MT")
    plotname =file.path(path,"figures/PCA_aut.pdf")
  }

  p = prcomp(scale(t(cpm_log[i,])))
  pc_scores = p$x
  pc_var = p$sdev^2
  pc_var_perc = pc_var/sum(pc_var)*100 

  df=data.frame(pc1=pc_scores[,1],pc2=pc_scores[,2],group=D_simple$cluster_group)
  df$group=factor(df$group, levels = c("male", "high", "highescape", "low"))
  levels(df$group) <- c("Male","Group 1","Group 2","Group 3")
  df$group = factor(df$group, levels = c("Male","Group 1", "Group 2", "Group 3"))
  
  pp1=ggplot(df,aes(x=pc1,y=pc2,col=group)) +
    geom_point(size=3,alpha=0.7) +
    xlab(paste("PC1 (",round(pc_var_perc[1],2),"%)",sep="")) +
    ylab(paste("PC2 (",round(pc_var_perc[2],2),"%)",sep="")) +
    theme_minimal() +
    scale_color_manual(name = "", values=c("#3399FF","#CC0033","mediumorchid","#FF99CC")) +
    theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22), 
        legend.title = element_blank(),
        legend.text = element_text(size=22),
        plot.title = element_text(size=22))
  if (chr == "X")
    pp1 = pp1 + labs(title="chrX genes")
  else
    pp1 = pp1 + labs(title="Autosomal genes") 

  ggsave(plotname, pp1, height=16,width=18,units="cm",limitsize = FALSE)

}