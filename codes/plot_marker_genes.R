plot_marker_genes <- function(path) {
  
  library(ggplot2)
  library(ggpubr)
  
  load(file.path(path,"data/data_for_DE_new.RData"))
  clusters=read.table(file.path(path,"data/female_clusters.txt"),header=TRUE,sep="\t")
  D_simple$xist_group=clusters$clusterName[match(D_simple$lines,clusters$lines)]
  D_simple$xist_group[which(D_simple$sex=="male")]="male"
  
  nanog_ind=which(gene_info$gene_name=="NANOG")
  sox2_ind=which(gene_info$gene_name=="SOX2")
  oct4_ind=which(gene_info$gene_name=="POU5F1") 
  
  df_nanog=data.frame(value=cpm_log[nanog_ind,],Sex=D_simple$xist_group,lines=D_simple$lines)
  df_sox2=data.frame(value=cpm_log[sox2_ind,],Sex=D_simple$xist_group,lines=D_simple$lines)
  df_oct4=data.frame(value=cpm_log[oct4_ind,],Sex=D_simple$xist_group,lines=D_simple$lines)

  plotname =file.path(path,"figures/NANOG.pdf")
  df=df_nanog
  df$Sex = factor(df$Sex, levels = c("male", "high", "highescape", "low"))
  levels(df$Sex) <- c("Male","Group 1","Group 2","Group 3")
  df$Sex = factor(df$Sex, levels = c("Male","Group 1", "Group 2", "Group 3"))

  my_comparisons <- list( c("Group 1", "Group 2"), c("Group 1", "Group 3"), c("Group 2", "Group 3"), c("Male", "Group 1"), c("Male", "Group 2"), c("Male", "Group 3"))
  p=ggplot(df, aes(x=Sex, y=value,color=Sex)) +
    geom_boxplot(outlier.shape=NA,size=1) +
    geom_point(position = position_jitterdodge(seed = 42),size=0.6) +
    theme_minimal() +
    #geom_jitter(aes(color=Sex),size=0.5,shape=19) +
    #scale_fill_discrete(labels=c("Female", "Male")) +
    ylab("NANOG expression (log-cpm)") +
    xlab("") +
    scale_color_manual(name = "", values=c("#3399FF","#CC0033","mediumorchid","#FF99CC")) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
    theme(text = element_text(size = 22)) +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=22), axis.title=element_text(size=22))
  ggsave(plotname, p, height=16, width=16, units="cm", limitsize = FALSE)

  plotname =file.path(path,"figures/POU5F1.pdf")
  df=df_oct4
  df$Sex = factor(df$Sex, levels = c("male", "high", "highescape", "low"))
  levels(df$Sex) <- c("Male","Group 1","Group 2","Group 3")
  df$Sex = factor(df$Sex, levels = c("Male","Group 1", "Group 2", "Group 3"))
  
  p=ggplot(df, aes(x=Sex, y=value,color=Sex)) +
    geom_boxplot(outlier.shape=NA,size=1) +
    geom_point(position = position_jitterdodge(seed = 42),size=0.6) +
    theme_minimal() +
    #geom_jitter(aes(color=Sex),size=0.5,shape=19) +
    #scale_fill_discrete(labels=c("Female", "Male")) +
    ylab("POU5F1 expression (log-cpm)") +
    xlab("") +
    scale_color_manual(name = "", values=c("#3399FF","#CC0033","mediumorchid","#FF99CC")) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
    theme(text = element_text(size = 22)) +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=22), axis.title=element_text(size=22))
  ggsave(plotname, p, height=16, width=16, units="cm", limitsize = FALSE)

  plotname =file.path(path,"figures/SOX2.pdf")
  df=df_sox2
  df$Sex = factor(df$Sex, levels = c("male", "high", "highescape", "low"))
  levels(df$Sex) <- c("Male","Group 1","Group 2","Group 3")
  df$Sex = factor(df$Sex, levels = c("Male","Group 1", "Group 2", "Group 3"))
  
  p=ggplot(df, aes(x=Sex, y=value,color=Sex)) +
    geom_boxplot(outlier.shape=NA,size=1) +
    geom_point(position = position_jitterdodge(seed = 42),size=0.6) +
    theme_minimal() +
    ylab("SOX2 expression (log-cpm)") +
    xlab("") +
    scale_color_manual(name = "", values=c("#3399FF","#CC0033","mediumorchid","#FF99CC")) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
    theme(text = element_text(size = 22)) +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=22), axis.title=element_text(size=22))
  ggsave(plotname, p, height=16, width=16, units="cm", limitsize = FALSE)
  
}