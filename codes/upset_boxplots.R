# -----------------------------------------------------------
# Script Name: upset_boxplots.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

upset_boxplots <- function(path,chr,plot.pvals=FALSE) {
  # set chr="X" or chr="aut"
  
  #library(ggpubr)
  #library(ggplot2)
  #library(dplyr)
  
  hipsci_datadir = file.path(path,"data")
  hipsci_resultsdir = file.path(path,"results")
  hipsci_figuresdir = file.path(path,"figures")
  if (plot.pvals==TRUE) {
    plotname=file.path(hipsci_figuresdir,paste("ipsc_boxplot_",chr,"_w_pvals.pdf",sep=""))
  } else {
    plotname=file.path(hipsci_figuresdir,paste("ipsc_boxplot_",chr,".pdf",sep=""))
  }
  load(file.path(hipsci_datadir,"data_for_DE_ipsc.RData"))
  clusters=read.table(file.path(hipsci_datadir,"female_clusters.txt"),header=TRUE,sep="\t")
  D_simple$xist_group=clusters$clusterName[match(D_simple$lines,clusters$lines)]
  D_simple$xist_group[which(D_simple$sex=="male")]="male"
  D_ipsc=D_simple
  ind_male=which(D_ipsc$xist_group=="male")
  ind_G1=which(D_ipsc$xist_group=="G1")
  ind_G2=which(D_ipsc$xist_group=="G2")
  ind_G3=which(D_ipsc$xist_group=="G3")
  inds=which(gene_info$gene_name %in% c("DNMT3A","DNMT3B"))
  df=data.frame(expr=c(cpm_log[inds[1],ind_male],cpm_log[inds[1],ind_G1],cpm_log[inds[1],ind_G2],cpm_log[inds[1],ind_G3],
                     cpm_log[inds[2],ind_male],cpm_log[inds[2],ind_G1],cpm_log[inds[2],ind_G2],cpm_log[inds[2],ind_G3]),
                gene=c(rep("DNMT3A",dim(D_simple)[1]),rep("DNMT3B",dim(D_simple)[1])),
                group=rep(c(rep("M",length(ind_male)),rep("G1",length(ind_G1)),rep("G2",length(ind_G2)),rep("G3",length(ind_G3))),2))
  my_comparisons <- list( c("M", "G1"), c("M", "G2"), c("M", "G3"), c("G1","G2"), c("G1","G3"), c("G2","G3") )
  p_a=ggplot(df[which(df$gene=="DNMT3A"),],aes(x=group,y=expr,col=group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=1.5) +
    #ggtitle("DNMT3A") +
    stat_compare_means(comparisons = my_comparisons,method = "t.test") +
    ylab("DNMT3A gene expression (log2 CPM)") +
    scale_color_manual(breaks=c("M","G1","G2","G3"),values=c("#3399FF","#CC0033","mediumorchid","#FF99CC")) +
    theme_minimal() +
    theme(legend.position="none") +
    xlab("") +
    theme(text = element_text(size=20))
  p_b=ggplot(df[which(df$gene=="DNMT3B"),],aes(x=group,y=expr,col=group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=1.5) +
    #ggtitle("DNMT3B") +
    stat_compare_means(comparisons = my_comparisons,method = "t.test") +
    ylab("DNMT3B gene expression (log2 CPM)") +
    scale_color_manual(breaks=c("M","G1","G2","G3"),values=c("#3399FF","#CC0033","mediumorchid","#FF99CC")) +
    theme_minimal() +
    theme(legend.position="none") +
    xlab("") +
    theme(text = element_text(size=20))
  ggsave("figures/dnmt3a_boxplot.pdf", p_a, width=15,height=15,units="cm",limitsize = FALSE)
  ggsave("figures/dnmt3b_boxplot.pdf", p_b, width=15,height=15,units="cm",limitsize = FALSE)
  
  f1=read.table(file.path(hipsci_resultsdir,"top.table_final_sva_dream_ipsc_MG1_FALSE_TRUE.txt"))
  f2=read.table(file.path(hipsci_resultsdir,"top.table_final_sva_dream_ipsc_MG2_FALSE_TRUE.txt"))
  f3=read.table(file.path(hipsci_resultsdir,"top.table_final_sva_dream_ipsc_MG3_FALSE_TRUE.txt"))
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
  cpm_log=cpm_log[ind_chr,]
  gene_info2=gene_info[ind_chr,]
  G=list(f1,f2,f3)
  M=matrix(FALSE,dim(f1)[1],(length(G)*2))
  colnames(M)=c("G1up","G1down","G2up","G2down","G3up","G3down")
  for (i in 1:length(G)) {
    inds_up=which((G[[i]]$adj.P.Val<0.05 & G[[i]]$logFC>0)==TRUE)
    inds_down=which((G[[i]]$adj.P.Val<0.05 & G[[i]]$logFC<0)==TRUE)
    M[inds_up,((i-1)*2+1)]=TRUE
    M[inds_down,((i-1)*2+2)]=TRUE
  }
  df=as.data.frame(M)
  categories= df %>% distinct()
  rm.ind.allfalse=which(apply(!categories,1,all))
  categories=categories[-rm.ind.allfalse,]
  n_categories=dim(categories)[1]
  print(paste("Number of categories: ",n_categories,sep=""))
  gene_cat_inds=matrix(0,dim(f1)[1], n_categories)
  for (i in 1:nrow(categories)) {
    gene_cat_inds[,i]=apply(M,1,function(x) {all(x==categories[i,])})
  }
  n_genes=matrix(0,n_categories,1)
  for (i in 1:n_categories) {
    n_genes[i]=length(which(gene_cat_inds[,i]==1))
  }
  oo=order(n_genes,decreasing=TRUE)
  if (chr=="X") {
    oo[8:9]=oo[9:8]
    oo[13:14]=oo[14:13]
  }
  categories=categories[oo,]
  n_genes=n_genes[oo]
  gene_cat_inds=gene_cat_inds[,oo]
  cpm_log_stdzd=t(scale(t(cpm_log)))
  pp=vector("list", n_categories)
  my_comparisons <- list( c("M", "G1"), c("M", "G2"), c("M", "G3"), c("G1","G2"), c("G1","G3"), c("G2","G3") )
  ttests_X=vector("list", n_categories)
  for (i in 1:n_categories) {
    ind=which(gene_cat_inds[,i]==1)
    if (length(ind)<2) {
      mean_gene_expr_males=mean(cpm_log_stdzd[ind,ind_male])
      mean_gene_expr_G1=mean(cpm_log_stdzd[ind,ind_G1])
      mean_gene_expr_G2=mean(cpm_log_stdzd[ind,ind_G2])
      mean_gene_expr_G3=mean(cpm_log_stdzd[ind,ind_G3])
    } else {
      mean_gene_expr_males=apply(cpm_log_stdzd[ind,ind_male],1,mean)
      mean_gene_expr_G1=apply(cpm_log_stdzd[ind,ind_G1],1,mean)
      mean_gene_expr_G2=apply(cpm_log_stdzd[ind,ind_G2],1,mean)
      mean_gene_expr_G3=apply(cpm_log_stdzd[ind,ind_G3],1,mean)
    }
    df1=data.frame(expr=c(mean_gene_expr_males,mean_gene_expr_G1,mean_gene_expr_G2,mean_gene_expr_G3),
                   group=rep(c("M","G1","G2","G3"),each=length(ind)))
    df1$group=factor(df1$group,levels=c("M","G1","G2","G3"))
    if (length(ind)>2) {
      ttests_X[[i]]=list(M_G1=t.test(df1$expr[which(df1$group=="M")],df1$expr[which(df1$group=="G1")],alternative="two.sided",paired=TRUE),
                         M_G2=t.test(df1$expr[which(df1$group=="M")],df1$expr[which(df1$group=="G2")],alternative="two.sided",paired=TRUE),
                         M_G3=t.test(df1$expr[which(df1$group=="M")],df1$expr[which(df1$group=="G3")],alternative="two.sided",paired=TRUE),
                         G1_G2=t.test(df1$expr[which(df1$group=="G1")],df1$expr[which(df1$group=="G2")],alternative="two.sided",paired=TRUE),
                         G1_G3=t.test(df1$expr[which(df1$group=="G1")],df1$expr[which(df1$group=="G3")],alternative="two.sided",paired=TRUE),
                         G2_G3=t.test(df1$expr[which(df1$group=="G2")],df1$expr[which(df1$group=="G3")],alternative="two.sided",paired=TRUE))
    }
    p=ggplot(df1, aes(x=group,y=expr,col=group)) +
      geom_boxplot(lwd=1) +
      theme_minimal() +
      scale_color_manual(breaks=c("M","G1","G2","G3"),values=c("#3399FF","#CC0033","mediumorchid","#FF99CC")) +
      xlab("") +
      ylab("") +
      theme(legend.position="none") +
      #ylab("Mean standardized gene expression") +
      theme(text = element_text(size=20)) 
    if (i==1) {
      if (plot.pvals==TRUE) {
        p = p + stat_compare_means(comparisons = my_comparisons,method = "t.test",paired=TRUE) + scale_y_continuous("", breaks=c(-1,-0.5,0,0.5,1), labels=c("-1","-0.5","0","0.5","1"))
      } else {
        p= p + scale_y_continuous("", breaks=c(-1,-0.5,0,0.5,1), labels=c("-1","-0.5","0","0.5","1"), limits=c(-1,1))
      }
    } else {
      if (plot.pvals==TRUE) {
        p = p + stat_compare_means(comparisons = my_comparisons,method = "t.test",paired=TRUE) + scale_y_continuous("", breaks=c(-1,-0.5,0,0.5,1), labels=c("","","","",""))
      } else {
        p= p + scale_y_continuous("", breaks=c(-1,-0.5,0,0.5,1), labels=c("","","","",""), limits=c(-1,1))
      }
    }
    print(p)
    pp[[i]]=p
  }
  pp1=ggarrange(plotlist=pp, ncol=n_categories, nrow = 1)
  ggsave(plotname, pp1, width=5*n_categories,height=12,units="cm",limitsize = FALSE)
}

path="~/Documents/Git/github-work/hipsci-xci-final2"
upset_boxplots(path,chr="X")
upset_boxplots(path,chr="aut")
