# -----------------------------------------------------------
# Script Name: ase_plots.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

ase_plots <- function() {  
  load("data/dff1.RData")
  load("data/dff2.RData")
  dff2=df.sig
  path = "hipsci-xci"
  source(file.path(path,"codes/get_ase_matrix.R"))
  res=get_ase_matrix(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=0)
  #source(file.path(path,"codes/plot_heatmap.R"))
  #res_10=plot_heatmap(res=res, min_nonna_num=10)

  gene_info=res$gene_info
  sig=res$sig
  sig_r=sig
  D=res$D
  ase_ratios=res$ase_ratios
  gene_minorCount=res$gene_minorCount
  gene_allCount=res$gene_allCount
  gene_num_variants=res$gene_num_variants
  n_genes=dim(sig)[1]
  n_samples=dim(sig)[2]

  df.ase=as.data.frame(ase_ratios[n_genes:1,1:n_samples])
  df.ase$gene_name=gene_info$gene_name[n_genes:1]
  df.ase$gene_name=factor(df.ase$gene_name, levels = df.ase$gene_name)
  df.ase.1=melt(df.ase)
  df.ase.1$cluster_group=rep(D$xist_group[1:n_samples],each=dim(ase_ratios)[1])
  df.ase.1$cluster_group=factor(df.ase.1$cluster_group, levels=c("Low-XIST female","High-XIST female"),labels=c("Low-XIST females","High-XIST females"))
  df.ase.1$XIST=rep(D$xist_cpm_log[1:n_samples],each=dim(ase_ratios)[1])


  df.sig=as.data.frame(sig[n_genes:1,1:n_samples])
  df.sig$gene_name=gene_info$gene_name[n_genes:1]
  df.sig$gene_name=factor(df.sig$gene_name, levels = df.sig$gene_name)
  df.sig.1=melt(df.sig)
  df.sig.1$cluster_group=rep(D$xist_group[1:n_samples],each=dim(sig)[1])
  df.sig.1$cluster_group=factor(df.sig.1$cluster_group, levels=c("Low-XIST female","High-XIST female"),labels=c("Low-XIST females","High-XIST females"))
  df.sig.1$XIST=rep(D$xist_cpm_log[1:n_samples],each=dim(sig)[1])
  df.sig.1$avg_num_variants_per_gene=rep(colMeans(gene_num_variants,na.rm=TRUE),each=dim(sig)[1])
  df.sig.1$median_num_variants_per_gene=rep(apply(gene_num_variants,2,median,na.rm=TRUE),each=dim(sig)[1])
  df.sig.1$median_ase=rep(apply(ase_ratios,2,median,na.rm=TRUE),each=dim(sig)[1])
  df.sig.1$num_inactive=rep(colSums(sig==FALSE,na.rm=TRUE),each=dim(sig)[1])
  df.sig.1$num_escape=rep(colSums(sig==TRUE,na.rm=TRUE),each=dim(sig)[1])





  df.new.2=data.frame(num_inactive=100*colSums(sig==FALSE,na.rm=TRUE)/(colSums(sig==FALSE,na.rm=TRUE)+colSums(sig==TRUE,na.rm=TRUE)),num_escape=100*colSums(sig==TRUE,na.rm=TRUE)/(colSums(sig==FALSE,na.rm=TRUE)+colSums(sig==TRUE,na.rm=TRUE)))
  df.new.2$lines=as.character(1:dim(sig)[2])
  df.new.2$lines=factor(df.new.2$lines, levels = df.new.2$lines)
  df.new.1=melt(df.new.2)
  df.new.1$variable=factor(df.new.1$variable, levels=c("num_inactive","num_escape"),labels = c("Inactive","Escape"))

  ind_high=which(D$xist_group=="High-XIST female")
  ind_low=which(D$xist_group=="Low-XIST female")


  clusters=read.table("/Users/topah/Desktop/hipsci_codes/data/female_clusters.txt",header=TRUE,sep="\t")



#df=data.frame(mean_ase=apply(ase_ratios,2,mean,na.rm=TRUE),esc_ratio=colSums(sig,na.rm=TRUE)/colSums(!is.na(sig)),xist=D$xist_cpm_log,xist_group=D$xist_group,clusters=clusters$clusterName[match(commonfiles,clusters$lines)])
  df=data.frame(mean_ase=apply(ase_ratios,2,mean,na.rm=TRUE),median_ase=apply(ase_ratios,2,median,na.rm=TRUE),esc_ratio=colSums(sig,na.rm=TRUE)/colSums(!is.na(sig)),xist=D$xist_cpm_log,xist_group=D$xist_group,clusters=clusters$clusterName[match(D$lines,clusters$lines)])
  df$xist_group=as.factor(df$xist_group)
  levels(df$xist_group)[levels(df$xist_group)=="High-XIST female"]="High-XIST"
  levels(df$xist_group)[levels(df$xist_group)=="Low-XIST female"]="Low-XIST"
  df$clusters=as.factor(df$clusters)
  levels(df$clusters)[levels(df$clusters)=="high"]="Group 1"
  levels(df$clusters)[levels(df$clusters)=="highescape"]="Group 2"
  levels(df$clusters)[levels(df$clusters)=="low"]="Group 3"

  library(ggpubr)

  p0=ggplot(df,aes(x=xist_group,y=median_ase,color=xist_group)) +
    geom_boxplot(outlier.shape = NA, size=1) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"), size=5) +
    theme_minimal() +
    scale_color_manual(values=c("#CC0033", "#FF99CC"),name = "Female\ncell lines") +
    xlab("") +
    ylab("Median ASE in chrX genes") +
    ylim(0,0.5) +
    theme(text = element_text(size=20)) 

  p1=ggplot(df,aes(x=xist_group,y=esc_ratio,color=xist_group)) +
    geom_boxplot(outlier.shape = NA,size=1) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"), size=5) +
    theme_minimal() +
    scale_color_manual(values=c("#CC0033", "#FF99CC"),name = "Female\ncell lines") + #,breaks=c("High-XIST female","Low-XIST female"), labels = c("High-XIST", "Low-XIST")) +
    xlab("") +
    ylab("Fraction of genes with\nASE > 0.1") +
    ylim(0,1) +
    theme(text = element_text(size=20)) 

  p2=ggplot(df,aes(x=xist_group,y=mean_ase,color=xist_group)) +
    geom_boxplot(outlier.shape = NA, size=1) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less"), size=5) +
    theme_minimal() +
    scale_color_manual(values=c("#CC0033", "#FF99CC"),name = "Female\ncell lines") +
    xlab("") +
    ylab("Mean ASE in chrX genes") +
    ylim(0,0.5) +
    theme(text = element_text(size=20)) 

  my_comparisons <- list( c("Group 1", "Group 2"),  c("Group 2", "Group 3") , c("Group 1", "Group 3"))
  p3=ggplot(df,aes(x=clusters,y=esc_ratio,color=clusters)) +
    geom_boxplot(outlier.shape = NA, size=1) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    stat_compare_means(comparisons = my_comparisons,method="wilcox.test", size=5) +
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c("0","0.25","0.5","0.75","1"), limits=c(0,1.25)) +
    theme_minimal() +
    scale_color_manual(name = "Female\ncell lines", values=c("#CC0033", "mediumorchid","#FF99CC")) + #,labels=c("Group 1","Group 2", "Group 3"), breaks=c("high", "highescape", "low")) +
    xlab("") +
    ylab("Fraction of genes with\nASE > 0.1") +
    theme(text = element_text(size=20)) 

  median(df$esc_ratio[which(df$clusters=="Group 1")])
  #[1] 0.1785131
  median(df$esc_ratio[which(df$clusters=="Group 2")])
  #[1] 0.3924051
  median(df$esc_ratio[which(df$clusters=="Group 3")])
  #[1] 0.4128024

  wilcox.test(df$esc_ratio[which(df$clusters=="Group 1")],df$esc_ratio[which(df$clusters=="Group 2")])$p.value
  #[1] 1.419761e-19
  wilcox.test(df$esc_ratio[which(df$clusters=="Group 1")],df$esc_ratio[which(df$clusters=="Group 3")])$p.value
  #[1] 7.927614e-20
  wilcox.test(df$esc_ratio[which(df$clusters=="Group 2")],df$esc_ratio[which(df$clusters=="Group 3")])$p.value
  #[1] 0.3171439

  #df=data.frame(mean_ase=apply(ase_ratios,2,mean,na.rm=TRUE),xist_group=D$xist_group,clusters=clusters$clusterName[match(commonfiles,clusters$lines)])
  my_comparisons <- list( c("Group 1", "Group 2"),  c("Group 2", "Group 3") , c("Group 1", "Group 3"))
  p4=ggplot(df,aes(x=clusters,y=mean_ase,color=clusters)) +
    geom_boxplot(outlier.shape = NA, size=1) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    stat_compare_means(comparisons = my_comparisons, size=5) +
    scale_y_continuous(breaks=c(0,0.25,0.5), labels=c("0","0.25","0.5"), limits=c(0,0.6)) +
    theme_minimal() +
    scale_color_manual(name = "Female\ncell lines", values=c("#CC0033", "mediumorchid","#FF99CC")) + # ,labels=c("high","highescape", "low"), breaks=c("high", "highescape", "low")) +
    xlab("") +
    ylab("Mean ASE") +
    theme(text = element_text(size=20)) 

  wilcox.test(df$mean_ase[which(df$clusters=="Group 1")],df$mean_ase[which(df$clusters=="Group 2")])$p.value
  #[1] 4.959605e-19
  wilcox.test(df$mean_ase[which(df$clusters=="Group 1")],df$mean_ase[which(df$clusters=="Group 3")])$p.value
  #[1] 1.590799e-19
  wilcox.test(df$mean_ase[which(df$clusters=="Group 2")],df$mean_ase[which(df$clusters=="Group 3")])$p.value
  #[1] 0.01889344

  ggsave("/Users/topah/Desktop/deneme/median_ase_xist_groups.pdf", p0, width=15,height=10,units="cm",limitsize = FALSE)
  ggsave("/Users/topah/Desktop/deneme/esc_ratio_xist_groups.pdf", p1, width=15,height=10,units="cm",limitsize = FALSE)
  ggsave("/Users/topah/Desktop/deneme/mean_ase_xist_groups.pdf", p2, width=15,height=10,units="cm",limitsize = FALSE)
  ggsave("/Users/topah/Desktop/deneme/esc_ratio_clusters.pdf", p3, width=15,height=10,units="cm",limitsize = FALSE)
  ggsave("/Users/topah/Desktop/deneme/mean_ase_clusters.pdf", p4, width=15,height=10,units="cm",limitsize = FALSE)

  p5=ggplot(df, aes(x=xist,y=esc_ratio)) +
    geom_point(aes(colour=clusters)) +
    theme_minimal() +
    xlab(expression("XIST expression (log"[2]*"(CPM+0.5))")) + 
    ylab("Fraction of genes with\nASE > 0.1") +
    geom_smooth(method = "loess",formula = y ~ x, color="gray8") +
    geom_vline(xintercept=1.5,linetype=2,color="gray8") +
    scale_color_manual(name = "Female\ncell lines", values=c("#CC0033", "mediumorchid","#FF99CC")) +
    ylim(0,1) +
    theme(text = element_text(size=20)) 

  ggsave("/Users/topah/Desktop/deneme/kmeans_cluster.pdf", p5, height=13,width=18,units="cm",limitsize = FALSE)

  #res_X=get_ase_info(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=10,gene_list=NULL,known_type=NULL) 
  min_nonna_num=10
  res_X=res
  ase_ratios=res_X$ase_ratios
  sig=res_X$sig
  gene_allCount=res_X$gene_allCount
  D=res_X$D
  gene_info=res_X$gene_info
  gene_ids=gene_info$gene_id
  n_genes=dim(sig)[1]
  n_samples=dim(sig)[2]

  ##
  DF.rowwise=data.frame(num_cells_escape=rowSums(sig,na.rm=TRUE),num_cells_nonescape=rowSums(!sig,na.rm=TRUE),num_cells_all=rowSums(!is.na(sig)),mean_gene_counts=apply(gene_allCount/gene_num_variants,1,mean,na.rm=TRUE))
  DF.rowwise$loggene_readCounts <- log10(DF.rowwise$mean_gene_counts)
  DF.rowwise=DF.rowwise[dim(DF.rowwise)[1]:1,]
  DF.rowwise$genes=rownames(DF.rowwise)
  DF.rowwise=DF.rowwise[c(n_genes:1),]
  DF.rowwise$genes=factor(DF.rowwise$genes, levels = DF.rowwise$genes)
  dff2$gene_id==rownames(DF.rowwise)
  dff2$avg_reads_per_gene=DF.rowwise$mean_gene_counts
  dff2$log10_avg_reads_per_gene=DF.rowwise$loggene_readCounts
  library(writexl)
  write_xlsx(dff2,"results/ASE_summary_X-chr_female_genes.xlsx")
  ##

  ttt=which((rowSums(!is.na(ase_ratios[,which(D$xist_group=="High-XIST female")])))>=min_nonna_num & (rowSums(!is.na(ase_ratios[,which(D$xist_group=="Low-XIST female")])))>=min_nonna_num)
  xist_ind=which(gene_ids==gene_info$gene_id[match("XIST",gene_info$gene_name)])
  if (!(xist_ind %in% ttt)) {
    ttt=append(ttt,xist_ind,after=(which(ttt>xist_ind)[1]-1))
  }

  tsig=sig[ttt,]
  tgene_allCount=gene_allCount[ttt,]
  tgene_num_variants=gene_num_variants[ttt,]
  DF.rowwise=data.frame(num_cells_escape=rowSums(tsig,na.rm=TRUE),num_cells_nonescape=rowSums(!tsig,na.rm=TRUE),num_cells_all=rowSums(!is.na(tsig)),mean_gene_counts=apply(tgene_allCount/tgene_num_variants,1,mean,na.rm=TRUE))
  DF.rowwise$loggene_readCounts <- log10(DF.rowwise$mean_gene_counts)
  DF.rowwise=DF.rowwise[dim(DF.rowwise)[1]:1,]
  DF.rowwise$genes=rownames(DF.rowwise)
  DF.rowwise$genes=factor(DF.rowwise$genes, levels = DF.rowwise$genes)

  p.1=ggplot(data=DF.rowwise, aes(y=genes, x=loggene_readCounts)) +
    geom_bar(stat="identity",fill="azure4") +
    theme_minimal() +
    xlab("Average total number of\nreads per gene (log10 scale)") +
    scale_x_continuous(limits = c(0,4), breaks = c(0,1,2,3,4) , labels = as.character(10^(c(0,1,2,3,4)))) +
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    theme(text = element_text(size=20)) 
  ggsave("figures/avg_tot_num_reads_per_gene_rows.pdf", p.1, height=44,width=10,units="cm",limitsize = FALSE)


  p.2=ggplot(data=DF.rowwise, aes(y=genes, x=num_cells_all)) +
    geom_bar(position='stack',stat="identity",fill="darkslategray4") +
    theme_minimal() +
    xlab("Number of cell lines observed") +
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    theme(text = element_text(size=20)) 
  ggsave("figures/num_cell_lines_rows.pdf", p.2, height=44,width=10,units="cm",limitsize = FALSE)


  DF.rowwise=data.frame(num_cells_escape=100*(rowSums(tsig,na.rm=TRUE)/rowSums(!is.na(tsig))),num_cells_nonescape=100*(rowSums(!tsig,na.rm=TRUE)/rowSums(!is.na(tsig))))
  DF.rowwise=DF.rowwise[dim(DF.rowwise)[1]:1,]
  DF.rowwise$genes=rownames(DF.rowwise)
  DF.rowwise$genes=factor(DF.rowwise$genes, levels = DF.rowwise$genes)

  DF.rowwise.1=melt(DF.rowwise)
  DF.rowwise.1$variable=factor(DF.rowwise.1$variable, levels=c("num_cells_nonescape","num_cells_escape"),labels = c("Inactive","Escape"))
  DF.rowwise.1$genes=factor(DF.rowwise.1$genes, levels = unique(DF.rowwise.1$genes))
  cols=c('blue', 'red')
  p.3=ggplot(data=DF.rowwise.1, aes(y=genes, x=value,fill=variable)) +
    geom_bar(position='stack',stat="identity")  +
    theme_minimal() +
    scale_fill_manual(values=cols) +
    xlab("Percent of cell lines") +
    guides(fill="none") +
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    theme(text = element_text(size=20)) 
  ggsave("figures/escape_percent_cell_lines_rows.pdf", p.3, height=44,width=10,units="cm",limitsize = FALSE)

  # For filtered #
  tsig=sig[ttt,]
  tgene_allCount=gene_allCount[ttt,]
  tgene_num_variants=gene_num_variants[ttt,]
  tsig=tsig[,dim(tsig)[2]:1]
  tgene_allCount=tgene_allCount[,dim(tsig)[2]:1]
  tgene_num_variants=tgene_num_variants[,dim(tsig)[2]:1]
  # For filtered end #

  # For not filtered #
  #tsig=sig[,dim(sig)[2]:1]
  #tgene_allCount=gene_allCount[,dim(sig)[2]:1]
  # For not filtered end #
  DF.columnwise=data.frame(num_genes_escape=colSums(tsig,na.rm=TRUE),num_genes_nonescape=colSums(!tsig,na.rm=TRUE),num_genes_all=colSums(!is.na(tsig)),mean_gene_counts=apply(tgene_allCount/tgene_num_variants,2,mean,na.rm=TRUE))
  DF.columnwise$cells=rownames(DF.columnwise)
  DF.columnwise$cells=factor(DF.columnwise$cells, levels = DF.columnwise$cells)

  dff1$avg_reads_per_gene_min10lines_filtered=DF.columnwise$mean_gene_counts
  library(writexl)
  write_xlsx(dff1,"results/ASE_summary_X-chr_female_iPSC_lines.xlsx")

  p.4=ggplot(data=DF.columnwise, aes(x=cells, y=mean_gene_counts)) +
    geom_bar(stat="identity",fill="azure4") +
    xlab("") +
    theme_minimal() +
    ylab("Average total number of\nreads per gene") +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
    theme(text = element_text(size=20)) 
  ggsave("/Users/topah/Desktop/deneme/avg_tot_num_reads_per_gene_columns.pdf", p.4, height=10,width=52,units="cm",limitsize = FALSE)


  p.5=ggplot(data=DF.columnwise, aes(x=cells, y=num_genes_all)) +
    geom_bar(position='stack',stat="identity",fill="darkslategray4") +
    theme_minimal() +
    ylab("Number of genes observed") +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
    theme(text = element_text(size=20)) 
  ggsave("figures/num_genes_columns.pdf", p.5, height=10,width=52,units="cm",limitsize = FALSE)


  DF.columnwise=data.frame(num_genes_escape=100*(colSums(tsig,na.rm=TRUE)/colSums(!is.na(tsig))),num_genes_nonescape=100*(colSums(!tsig,na.rm=TRUE)/colSums(!is.na(tsig))))
  DF.columnwise$cells=rownames(DF.columnwise)
  DF.columnwise$cells=factor(DF.columnwise$cells, levels = DF.columnwise$cells)
  DF.columnwise.1=melt(DF.columnwise)
  DF.columnwise.1$variable=factor(DF.columnwise.1$variable, levels=c("num_genes_nonescape","num_genes_escape"),labels = c("Inactive","Escape"))
  DF.columnwise.1$cells=factor(DF.columnwise.1$cells, levels = unique(DF.columnwise.1$cells))
  cols=c('blue', 'red')
  p.6=ggplot(data=DF.columnwise.1, aes(x=cells, y=value,fill=variable)) +
    geom_bar(position='stack',stat="identity")  +
    theme_minimal() +
    scale_fill_manual(values=cols) +
    ylab("Percent of genes") +
    guides(fill="none") +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
    theme(text = element_text(size=20)) 
  ggsave("figures/escape_percent_genes_columns.pdf", p.6, height=10,width=52,units="cm",limitsize = FALSE)

  cor.test(dff1$num_observed_genes_min10lines_filtered, dff1$XIST_logCPM, method="spearman", use = "complete.obs")
  cor.test(dff1$avg_reads_per_gene_min10lines_filtered, dff1$XIST_logCPM, method="spearman", use = "complete.obs") 


  # overall_p <- function(my_model) {
  #   f <- summary(my_model)$fstatistic
  #   p <- pf(f[1],f[2],f[3],lower.tail=F)
  #   attributes(p) <- NULL
  #   return(p)
  # }
  # overall_p(lm(dff1$num_observed_genes_min10lines_filtered~dff1$XIST_logCPM))
  # overall_p(lm(dff1$avg_reads_per_gene_min10lines_filtered~dff1$XIST_logCPM))

  #cor.test(dff1$num_observed_genes_min10lines_filtered, dff1$XIST_logCPM, method=c("pearson"))
  #cor.test(dff1$avg_reads_per_gene_min10lines_filtered, dff1$XIST_logCPM, method=c("pearson")) 

  }
