# -----------------------------------------------------------
# Script Name: get_ase_summary_data.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

get_ase_summary_data <- function(path) {
  source(file.path(path,"codes/get_ase_matrix.R"))
  res=get_ase_matrix(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=0)
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

  library("readxl")
  escape_genes=as.data.frame(read_excel("data/escape_genes.xlsx",col_names = TRUE, skip = 1))
  human_xci=escape_genes$`Combined XCI status`[match(gene_info$gene_id,escape_genes$`Gene ID`)]
  human_xci[which(is.na(human_xci))]="Unknown"
  gene_info$human_xci=human_xci

  sig_r=as.data.frame(sig_r)
  sig_r$gene_biotype=gene_info$gene_biotype
  sig_r$end_position=gene_info$end_position
  sig_r$start_position=gene_info$end_position
  sig_r$gene_name=gene_info$gene_name
  sig_r$gene_id=gene_info$gene_id
  sig_r$chr="X"
  sig_r=sig_r[,dim(sig_r)[2]:1]

  library(writexl)
  write_xlsx(sig_r,paste("results/escape_Xchr",".xlsx",sep=""))

  clusters=read.table("data/female_clusters.txt",header=TRUE,sep="\t")
  clusts=clusters$clusterName[match(D$lines,clusters$lines)]

  ddf=data.frame(num_escape_all=rowSums(sig,na.rm=TRUE),
               num_observed_all=rowSums(!is.na(sig)),
               num_lines_all=dim(D)[1],
               num_escape_highXIST=rowSums(sig[,D$xist_group=="High-XIST female"],na.rm=TRUE),
               num_observed_highXIST=rowSums(!is.na(sig[,D$xist_group=="High-XIST female"])),
               num_lines_highXIST=sum(D$xist_group=="High-XIST female"),
               num_escape_lowXIST=rowSums(sig[,D$xist_group=="Low-XIST female"],na.rm=TRUE),
               num_observed_lowXIST=rowSums(!is.na(sig[,D$xist_group=="Low-XIST female"])),
               num_lines_lowXIST=sum(D$xist_group=="Low-XIST female"),
               num_escape_Group1=rowSums(sig[,clusts=="G1"],na.rm=TRUE),
               num_observed_Group1=rowSums(!is.na(sig[,clusts=="G1"])),
               num_lines_Group1=sum(clusts=="G1"),
               num_escape_Group2=rowSums(sig[,clusts=="G2"],na.rm=TRUE),
               num_observed_Group2=rowSums(!is.na(sig[,clusts=="G2"])),
               num_lines_Group2=sum(clusts=="G2"),
               num_escape_Group3=rowSums(sig[,clusts=="G3"],na.rm=TRUE),
               num_observed_Group3=rowSums(!is.na(sig[,clusts=="G3"])),
               num_lines_Group3=sum(clusts=="G3"))
  ddf$fraction_escape_all=ddf$num_escape_all/ddf$num_observed_all
  ddf$fraction_escape_highXIST=ddf$num_escape_highXIST/ddf$num_observed_highXIST
  ddf$fraction_escape_lowXIST=ddf$num_escape_lowXIST/ddf$num_observed_lowXIST
  ddf$fraction_escape_Group1=ddf$num_escape_Group1/ddf$num_observed_Group1
  ddf$fraction_escape_Group2=ddf$num_escape_Group2/ddf$num_observed_Group2
  ddf$fraction_escape_Group3=ddf$num_escape_Group3/ddf$num_observed_Group3

  t.test(colSums(sig[,which(D$xist_group=="High-XIST female")],na.rm=TRUE)/colSums(!is.na(sig[,which(D$xist_group=="High-XIST female")]),na.rm=TRUE),
       colSums(sig[,which(D$xist_group=="Low-XIST female")],na.rm=TRUE)/colSums(!is.na(sig[,which(D$xist_group=="Low-XIST female")]),na.rm=TRUE))$p.value # 1.308153e-08
  t.test(colSums(sig[,which(D$xist_group=="High-XIST female")],na.rm=TRUE)/colSums(!is.na(sig[,which(D$xist_group=="High-XIST female")]),na.rm=TRUE))
  t.test(colSums(sig[,which(D$xist_group=="Low-XIST female")],na.rm=TRUE)/colSums(!is.na(sig[,which(D$xist_group=="Low-XIST female")]),na.rm=TRUE))

  ind=c(which(gene_info$gene_name=="XIST"),which(ddf$num_observed_highXIST>=10 & ddf$num_observed_lowXIST>=10))
  ind=sort(ind)
  mean(ddf$fraction_escape_highXIST[ind],na.rm=TRUE)

  ddf.small=ddf[ind,]
  sum(ddf.small$num_escape_highXIST>0 | ddf.small$num_escape_lowXIST>0)

  range(colSums(sig[ind,],na.rm=TRUE)/colSums(!is.na(sig[ind,]),na.rm=TRUE))
  range(rowSums(sig[ind,],na.rm=TRUE)/rowSums(!is.na(sig[ind,]),na.rm=TRUE))

  sum(rowSums(sig[ind,],na.rm=TRUE)/rowSums(!is.na(sig[ind,]),na.rm=TRUE)>0.25)
  sum(rowSums(sig[ind,],na.rm=TRUE)/rowSums(!is.na(sig[ind,]),na.rm=TRUE)>0.5)

  min_nonna_num=10
  mychr="X"
  gene_ids=gene_info$gene_id
  ttt=which((rowSums(!is.na(ase_ratios[,which(D$xist_group=="High-XIST female")])))>=min_nonna_num & (rowSums(!is.na(ase_ratios[,which(D$xist_group=="Low-XIST female")])))>=min_nonna_num)
  if (mychr=="X") {
    xist_ind=which(gene_ids==gene_info$gene_id[match("XIST",gene_info$gene_name)])
    if (!(xist_ind %in% ttt)) {
      ttt=append(ttt,xist_ind,after=(which(ttt>xist_ind)[1]-1))
    }
  }
  print(length(ttt))

  sex_groups=names(table(D$xist_group))
  wilcox.pvalue=matrix(NA,n_genes,1)
  fisher.pvalue=matrix(NA,n_genes,1)
  esc_diff=matrix(NA,n_genes,1)
  for (i in 1:n_genes) {
    M <- as.table(rbind(c(sum(sig[i,which(D$xist_group==sex_groups[1])],na.rm=TRUE),sum(!sig[i,which(D$xist_group==sex_groups[1])],na.rm=TRUE)),
                      c(sum(sig[i,which(D$xist_group==sex_groups[2])],na.rm=TRUE),sum(!sig[i,which(D$xist_group==sex_groups[2])],na.rm=TRUE))))
    dimnames(M) <- list(xist_group = c(sex_groups[1],sex_groups[2]),
                      status = c("Escape","Inactive"))
    esc_diff[i]=(M[,1]/rowSums(M))[2]-(M[,1]/rowSums(M))[1]
    tryCatch({
      fisher.pvalue[i]=fisher.test(M)$p.value
    }, error=function(e){})
    tryCatch({
      wilcox.pvalue[i]=wilcox.test(ase_ratios[i,which(D$xist_group==sex_groups[2])],ase_ratios[i,which(D$xist_group==sex_groups[1])],exact=FALSE)$p.value
    }, error=function(e){})
  }

  ddf$diff_fraction_escape_lowXIST_highXIST=ddf$fraction_escape_lowXIST-ddf$fraction_escape_highXIST
  ddf$fisher.pvalue=fisher.pvalue
  ddf$median_ASE_highXIST=apply(ase_ratios[,which(D$xist_group=="High-XIST female")],1,median,na.rm=TRUE)
  ddf$median_ASE_lowXIST=apply(ase_ratios[,which(D$xist_group=="Low-XIST female")],1,median,na.rm=TRUE)
  ddf$diff_median_ASE_lowXIST_highXIST=ddf$median_ASE_lowXIST-ddf$median_ASE_highXIST
  ddf$wilcox.pvalue=wilcox.pvalue
  ddf$median_ASE_Group1=apply(ase_ratios[,which(clusts=="high")],1,median,na.rm=TRUE)
  ddf$median_ASE_Group2=apply(ase_ratios[,which(clusts=="highescape")],1,median,na.rm=TRUE)
  ddf$median_ASE_Group3=apply(ase_ratios[,which(clusts=="low")],1,median,na.rm=TRUE)

  ddf$fisher.pvalue[setdiff(1:dim(ddf)[1],ttt)]=NA
  ddf$wilcox.pvalue[setdiff(1:dim(ddf)[1],ttt)]=NA

  sum(ddf$fisher.pvalue<0.05,na.rm=TRUE) # 41
  sum(ddf$wilcox.pvalue<0.05,na.rm=TRUE) # 58

  ase_summary_by_genes=cbind(gene_info,ddf)

  names(ase_summary_by_genes)[which(names(ase_summary_by_genes)=="wilcox.pvalue")]="diff_median_ASE_wilcox.pvalue"
  names(ase_summary_by_genes)[which(names(ase_summary_by_genes)=="fisher.pvalue")]="diff_fraction_escape_fisher.pvalue"

  #write_xlsx(ase_summary_by_genes,"ASE_summary_X-chr_genes.xlsx")
  save(ase_summary_by_genes, file = "data/genes_ase_summary_data.RData")
  
  ase_summary_by_lines=data.frame(cell_line=D$lines,XIST_logCPM=D$xist_cpm_log,XIST_group=D$xist_group,XIST_cluster=clusts)
  ase_summary_by_lines$mean_ASE=apply(ase_ratios,2,mean,na.rm=TRUE)
  ase_summary_by_lines$median_ASE=apply(ase_ratios,2,median,na.rm=TRUE)
  ase_summary_by_lines$num_escape_genes=colSums(sig,na.rm=TRUE)
  ase_summary_by_lines$num_observed_genes=colSums(!is.na(sig))
  ase_summary_by_lines$fraction_of_genes_escaping=colSums(sig,na.rm=TRUE)/colSums(!is.na(sig))
  ase_summary_by_lines$XIST_cluster[which(ase_summary_by_lines$XIST_cluster=="G3")]="Group3"
  ase_summary_by_lines$XIST_cluster[which(ase_summary_by_lines$XIST_cluster=="G2")]="Group2"
  ase_summary_by_lines$XIST_cluster[which(ase_summary_by_lines$XIST_cluster=="G1")]="Group1"

  ase_summary_by_lines=ase_summary_by_lines[dim(ase_summary_by_lines)[1]:1,]

  t.test(ase_summary_by_lines$fraction_of_genes_escaping[which(ase_summary_by_lines$XIST_group=="Low-XIST female")],ase_summary_by_lines$fraction_of_genes_escaping[which(ase_summary_by_lines$XIST_group=="High-XIST female")]) # Figure 2A
  wilcox.test(ase_summary_by_lines$fraction_of_genes_escaping[which(ase_summary_by_lines$XIST_cluster=="Group1")],ase_summary_by_lines$fraction_of_genes_escaping[which(ase_summary_by_lines$XIST_cluster=="Group2")])$p.value # Figure 3B

  median(ase_summary_by_lines$median_ASE[which(ase_summary_by_lines$XIST_group=="Low-XIST female")]) # 0.0499343
  median(ase_summary_by_lines$median_ASE[which(ase_summary_by_lines$XIST_group=="High-XIST female")]) # 0.01211148
  median(ase_summary_by_lines$mean_ASE[which(ase_summary_by_lines$XIST_group=="High-XIST female")]) # 0.1014025
  median(ase_summary_by_lines$mean_ASE[which(ase_summary_by_lines$XIST_group=="Low-XIST female")]) # 0.1790307

  xist.ind=which(gene_info$gene_name=="XIST")
  pval.xist=matrix(NA,n_samples,1)
  p_est=matrix(NA,n_samples,1)
  low_xist=matrix(NA,n_samples,1)
  up_xist=matrix(NA,n_samples,1)
  for (j in 1:n_samples) {
  if (!is.na(gene_minorCount[xist.ind,j])) {
      p_est[j]=(gene_minorCount[xist.ind,j]/gene_allCount[xist.ind,j])
      low_xist[j]=p_est[j]-sqrt(p_est[j]*(1-p_est[j])/gene_allCount[xist.ind,j])*qnorm(1-0.05/2)
      up_xist[j]=p_est[j]+sqrt(p_est[j]*(1-p_est[j])/gene_allCount[xist.ind,j])*qnorm(1-0.05/2)
      pval.xist[j]=binom.test(gene_minorCount[xist.ind,j],gene_allCount[xist.ind,j],p=0.05,alternative="less")$p.value
    }
  }

  ase_summary_by_lines=data.frame(cell_line=D$lines,XIST_logCPM=D$xist_cpm_log,XIST_group=D$xist_group,XIST_cluster=clusts)
  ase_summary_by_lines$XIST_ASE=ase_ratios[xist_ind,]
  ase_summary_by_lines$XIST_ASE_CI_lower=low_xist
  ase_summary_by_lines$XIST_ASE_CI_upper=up_xist
  ase_summary_by_lines$XIST_is_monoallelic_binom.pval=pval.xist
  ase_summary_by_lines$mean_ASE=apply(ase_ratios,2,mean,na.rm=TRUE)
  ase_summary_by_lines$median_ASE=apply(ase_ratios,2,median,na.rm=TRUE)
  ase_summary_by_lines$num_escape_genes=colSums(sig,na.rm=TRUE)
  ase_summary_by_lines$num_observed_genes=colSums(!is.na(sig))
  ase_summary_by_lines$fraction_of_genes_escaping=colSums(sig,na.rm=TRUE)/colSums(!is.na(sig))
  ase_summary_by_lines$num_escape_genes_min10lines_filtered=colSums(sig[ttt,],na.rm=TRUE)
  ase_summary_by_lines$num_observed_genes_min10lines_filtered=colSums(!is.na(sig[ttt,]))
  ase_summary_by_lines$fraction_of_genes_escaping_min10lines_filtered=colSums(sig[ttt,],na.rm=TRUE)/colSums(!is.na(sig[ttt,]))
  ase_summary_by_lines$XIST_cluster[which(ase_summary_by_lines$XIST_cluster=="G3")]="Group3"
  ase_summary_by_lines$XIST_cluster[which(ase_summary_by_lines$XIST_cluster=="G2")]="Group2"
  ase_summary_by_lines$XIST_cluster[which(ase_summary_by_lines$XIST_cluster=="G1")]="Group1"

  ase_summary_by_lines=ase_summary_by_lines[dim(ase_summary_by_lines)[1]:1,]
  #write_xlsx(ase_summary_by_lines,"ASE_summary_X-chr_female_iPSC_lines.xlsx")
  sum(ase_summary_by_lines$XIST_is_monoallelic_binom.pval<0.05,na.rm=TRUE)/sum(!is.na(ase_summary_by_lines$XIST_is_monoallelic_binom.pval),na.rm=TRUE) # 0.75
  plot(ase_summary_by_lines$XIST_ASE,ase_summary_by_lines$median_ASE)
  save(ase_summary_by_lines, file = "data/lines_ase_summary_data.RData")
}
