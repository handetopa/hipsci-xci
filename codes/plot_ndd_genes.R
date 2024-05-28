# -----------------------------------------------------------
# Script Name: plot_ndd_genes.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

select_genes <- function(df=mh,ndd_genes) {
  logFC=(-df$logFC[match(ndd_genes,df$gene_name)])
  logFC_SE=(df$logFC_SE[match(ndd_genes,df$gene_name)])
  adj.p.val=(df$adj.P.Val[match(ndd_genes,df$gene_name)])
  chr=(df$chrs[match(ndd_genes,df$gene_name)])
  xci_class=(df$class[match(ndd_genes,df$gene_name)])
  logFC_lower=logFC-(1.96*logFC_SE)
  logFC_upper=logFC+(1.96*logFC_SE)
  iii=which((pmax(abs(logFC_lower),abs(logFC_upper))>=0.59) & (adj.p.val<0.05))
  return(iii)
}

ndd_genes_ipsc <- function(path) {
  hipsci_datadir = file.path(path,"data")
  hipsci_resultsdir = file.path(path,"results")

  load(file.path(hipsci_datadir,"data_for_DE_ipsc.RData"))
  clusters=read.table(file.path(hipsci_datadir,"female_clusters.txt"),header=TRUE,sep="\t")
  D_simple$xist_group=clusters$clusterName[match(D_simple$lines,clusters$lines)]
  D_simple$xist_group[which(D_simple$sex=="male")]="male"
  
  mG1=read.table(file.path(hipsci_resultsdir,"top.table_final_sva_dream_ipsc_MG1_FALSE_TRUE.txt"))
  mG2=read.table(file.path(hipsci_resultsdir,"top.table_final_sva_dream_ipsc_MG2_FALSE_TRUE.txt"))
  mG3=read.table(file.path(hipsci_resultsdir,"top.table_final_sva_dream_ipsc_MG3_FALSE_TRUE.txt"))

  ndd_genes=read.table(file.path(hipsci_datadir,"NDD_genes.csv"),sep=",",header=TRUE)
  ndd_genes=unique(ndd_genes$gene)

  iii=select_genes(mG1,ndd_genes)
  iii=append(iii,select_genes(mG2,ndd_genes))
  iii=append(iii,select_genes(mG3,ndd_genes))
  iii=unique(iii)
  
  genes=ndd_genes[iii]

  ind=match(genes,gene_info$gene_name)
  cpm_log=cpm_log[ind,]
  gene_info=gene_info[ind,]
  xist_group=D_simple$xist_group
  
  n_samples=dim(cpm_log)[2]
  n_genes=dim(cpm_log)[1]
  df1=data.frame(expr=c(cpm_log),gene_name=rep(genes,n_samples),xist_group=rep(xist_group,each=n_genes))
  df1$xist_group=factor(df1$xist_group,levels=c("male","G1","G2","G3"),labels=c("Males","Group 1","Group 2", "Group 3"))

  #library(patchwork)
  #library(ggplot2)
  num_plots=0
  for (i in 1:length(genes)) {
    df2=df1[which(df1$gene_name==genes[i]),]
    pp1=ggplot(df2,aes(x=xist_group,y=expr,fill=xist_group)) +
      geom_boxplot(outlier.shape = NA) +
      theme_minimal() +
      ggtitle(paste(genes[i],", Chr: ",gene_info$chr[match(genes[i],gene_info$gene_name)],sep=""))  +
      ylab("") +
      xlab("") +
      #ylab("Gene expression") +
      #xlab("Sex groups") +
      geom_jitter(shape=16, size=1.5, position=position_jitter(0.2))  +
      scale_fill_manual(name = "", breaks=c("Males","Group 1","Group 2","Group 3"),values=c("#3399FF","#CC0033","mediumorchid","#FF99CC")) +
      theme(text = element_text(size=25), axis.text.x = element_text(size=25),axis.text.y = element_text(size=25))
    num_plots=num_plots+1
    if (num_plots==1) {
      pp=pp1
    } else {
      pp=pp+pp1
    }
  }
  p_final = pp + plot_layout(ncol = 4)
  ggsave("figures/ndd_genes_expr_ipsc_minFC1point5.pdf",p_final,width=100,height=160,units="cm",limitsize = FALSE)
  return(genes)
}

ndd_genes_neuron <- function(path) {
  hipsci_datadir = file.path(path,"data")
  hipsci_resultsdir = file.path(path,"results")
  
  load(file.path(hipsci_datadir,"data_for_DE_difneuron.RData"))
  mh=read.table(file.path(hipsci_resultsdir,"top.table_final_sva_dream_neuron_MH_TRUE_TRUE.txt"))
  ml=read.table(file.path(hipsci_resultsdir,"top.table_final_sva_dream_neuron_ML_TRUE_TRUE.txt"))
 
  ndd_genes=read.table(file.path(hipsci_datadir,"NDD_genes.csv"),sep=",",header=TRUE)
  ndd_genes=unique(ndd_genes$gene)
  
  iii=select_genes(mh,ndd_genes)
  iii=append(iii,select_genes(ml,ndd_genes))
  iii=unique(iii)
  
  genes=ndd_genes[iii]
  
  ind=match(genes,gene_info$gene_name)
  cpm_log=cpm_log[ind,]
  gene_info=gene_info[ind,]
  xist_group=D_simple$xist_group
  
  n_samples=dim(cpm_log)[2]
  n_genes=dim(cpm_log)[1]
  df1=data.frame(expr=c(cpm_log),gene_name=rep(genes,n_samples),xist_group=rep(xist_group,each=n_genes))
  df1$xist_group=factor(df1$xist_group,levels=c("male","highXIST","lowXIST"),labels=c("Males","High-XIST\nfemales","Low-XIST\nfemales"))

  num_plots=0
  for (i in 1:length(genes)) {
    df2=df1[which(df1$gene_name==genes[i]),]
    pp1=ggplot(df2,aes(x=xist_group,y=expr,fill=xist_group)) +
      geom_boxplot(outlier.shape = NA) +
      theme_minimal() +
      ggtitle(paste(genes[i],", Chr: ",gene_info$chr[match(genes[i],gene_info$gene_name)],sep=""))  +
      ylab("") +
      xlab("") +
      #ylab("Gene expression") +
      #xlab("Sex groups") +
      geom_jitter(shape=16, size=2.5, position=position_jitter(0.2))  +
      scale_fill_manual(name = "", breaks=c("Males","High-XIST\nfemales","Low-XIST\nfemales"),values=c("#3399FF","#CC0033","#FF99CC")) +
      theme(text = element_text(size=25), axis.text.x = element_text(size=25),axis.text.y = element_text(size=25))
    num_plots=num_plots+1
    if (num_plots==1) {
      pp=pp1
    } else {
      pp=pp+pp1
    }
  }
  p_final = pp + plot_layout(ncol = 4)
  ggsave("figures/ndd_genes_expr_neuron_minFC1point5.pdf",p_final,width=90,height=80,units="cm",limitsize = FALSE)
  return(genes)
}

ndd_genes_ipsc_ase <- function(path, ipsc_genes) {
  #source(file.path(path,"codes/get_ase_matrix.R"))
  res=get_ase_matrix(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=0)
  ase_ratios=res$ase_ratios
  D=res$D
  
  ase_ndd_genes=c(na.omit(res$gene_info$gene_name[match(ipsc_genes,res$gene_info$gene_name)]))
  ind=match(ase_ndd_genes,res$gene_info$gene_name)
  
  clusters=read.table(file.path(hipsci_datadir,"female_clusters.txt"),header=TRUE,sep="\t")
  D$xist_group=clusters$clusterName[match(D$lines,clusters$lines)]
  
  df1=data.frame(gene_name=rep(ase_ndd_genes,dim(D)[1]),ase=c(ase_ratios[ind,]),xist_group=rep(D$xist_group,each=length(ind)))
  df1$xist_group=factor(df1$xist_group,levels=c("G1","G2","G3"),labels=c("Group 1","Group 2", "Group 3"))
  
  ase_ndd_genes=intersect(intersect(df1$gene_name[(!is.na(df1$ase) & df1$xist_group=="Group 1")],
         df1$gene_name[(!is.na(df1$ase) & df1$xist_group=="Group 2")]),
         df1$gene_name[(!is.na(df1$ase) & df1$xist_group=="Group 3")])
  
  num_plots=0
  for (i in 1:length(ase_ndd_genes)) {
    df2=df1[which(df1$gene_name==ase_ndd_genes[i]),]
    pp1=ggplot(df2,aes(x=xist_group,y=ase,fill=xist_group)) +
      geom_boxplot(outlier.shape = NA) +
      theme_minimal() +
      ggtitle(paste(ase_ndd_genes[i],", Chr: X",sep=""))  +
      ylab("") +
      xlab("") +
      ylim(0,0.5) +
      #ylab("ASE") +
      #xlab("Sex groups") +
      geom_jitter(shape=16, size=2.5, position=position_jitter(0.2)) +
      scale_fill_manual(name = "", breaks=c("Group 1","Group 2","Group 3"),values=c("#CC0033","mediumorchid","#FF99CC")) +
      theme(text = element_text(size=25), axis.text.x = element_text(size=25),axis.text.y = element_text(size=25))
    num_plots=num_plots+1
    if (num_plots==1) {
      pp=pp1
    } else {
      pp=pp+pp1
    }
  }
  p_final = pp + plot_layout(ncol = 3)
  ggsave("figures/ndd_genes_ase_ipsc_minFC1point5.pdf",p_final,width=75,height=40,units="cm",limitsize = FALSE)
}

ipsc_genes=ndd_genes_ipsc(path)
ndd_genes_neuron(path)
ndd_genes_ipsc_ase(path, ipsc_genes)
