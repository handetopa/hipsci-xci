get_ase_matrix <- function(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=10) {

  library(qvalue)
  library(ggplot2)
  library(reshape)
  library(viridis)
  library(ggpubr)
  
  source("/Users/topah/Desktop/hipsci_codes/codes/position_scripts.R")
  source("/Users/topah/Desktop/hipsci_codes/codes/extract_text_before_first_dot.R")
  ga_path="/Users/topah/Desktop/hipsci_codes/data/ase_files"
  load("/Users/topah/Desktop/hipsci_codes/data/gene_info/ref_hs37d5.RData")
  load("/Users/topah/Desktop/hipsci_codes/data/data_for_DE_new.RData")
  DD=D_simple
  rm(D_simple)
  
  DD$xist_group=DD$sex
  DD$xist_group[which(DD$sex=="female" & DD$xist_cpm_log<xist_lim)]="Low-XIST female"
  DD$xist_group[which(DD$sex=="female" & DD$xist_cpm_log>=xist_lim)]="High-XIST female"
  DD$xist_group[which(DD$sex=="male")]="Male"
  
  print(table(DD$xist_group))
  
  if (include!="all") {
    D=DD[which(DD$xist_status==include),]
  } else {
    D=DD
  }
  
  if (mysex!="all") {
    ind_sex=which(D$sex==mysex)
    D=D[ind_sex,]
    commonfiles=D$lines
  }
  
  if (samples_orderby=="xist") {
    oo=order(D$xist_cpm_log)
    D=D[oo,]
    commonfiles=D$lines
  }
  
  print(table(D$xist_group))
  #write.table(D,file="/Users/topah/Desktop/ase_plots/HipSci_metadata.txt",quote=FALSE,col.names = TRUE, row.names = FALSE, sep="\t")
  
  #load("/Users/topah/Desktop/ipsc_xist.RData")
  #D$xist_cpm_log=xist2[match(D$lines,names(xist2))]
  
  genewise_ASEfiles=paste(commonfiles,"_min20_snpwiseASE_withGenePval.txt",sep="")
  genewise_ASEfiles=file.path(path=ga_path, genewise_ASEfiles)
  
  gw_ase=read.table(genewise_ASEfiles[1],header=TRUE)
  gw_ase=gw_ase[which(gw_ase$chromosome==as.character(mychr)),]
  gene_ids=gw_ase$gene_id
  
  for (i in (2:length(genewise_ASEfiles))) {
    gw_ase=read.table(genewise_ASEfiles[i],header=TRUE)
    gw_ase=gw_ase[which(gw_ase$chromosome==as.character(mychr)),]
    gene_ids=unique(append(gw_ase$gene_id,gene_ids))
  }
  
  # Remove gene biotypes other than protein coding and lincRNA
  gene_biotype=ref$gene_biotype[match(gene_ids,ref$ensembl_gene_id)]
  gene_ids=gene_ids[which(gene_biotype %in% c("protein_coding","lincRNA"))]
  print(length(gene_ids))
  
  gene_info=data.frame(gene_id=gene_ids,
                       gene_name=ref$external_gene_id[match(gene_ids,ref$ensembl_gene_id)],
                       start_position=ref$start_position[match(gene_ids,ref$ensembl_gene_id)],
                       end_position=ref$end_position[match(gene_ids,ref$ensembl_gene_id)],
                       gene_biotype=ref$gene_biotype[match(gene_ids,ref$ensembl_gene_id)])
  
  aa=gene_info$gene_name[which(gene_info$gene_name %in% names(which(table(gene_info$gene_name)!=1)))]
  if (length(aa)>0) {
    sel_aa=setdiff(1:length(aa),match(unique(aa),aa))
    gene_info=gene_info[-which(gene_info$gene_name %in% names(which(table(gene_info$gene_name)!=1)))[sel_aa],] # Remove IDS which don't exist in DE analysis.
    #gene_info$gene_name[which(gene_info$gene_name %in% names(which(table(gene_info$gene_name)!=1)))]=paste(aa," (",as.character(1:length(aa)),")",sep="")
  }
  
  if (genes_orderby=="pos") {
    oo=order(gene_info$start_position)
    gene_info=gene_info[oo,]
    gene_ids=gene_info$gene_id
  }
  
  print ("Before filtering:")
  print(paste("Number of genes: ",length(gene_ids),sep=""))
  print(paste("Number of cell lines: ",length(genewise_ASEfiles),sep=""))
  
  ase_ratios=matrix(0,length(gene_ids),length(genewise_ASEfiles))
  var_snps=matrix(NA,length(gene_ids),length(genewise_ASEfiles))
  var_snps_xist=matrix(NA,length(gene_ids),length(genewise_ASEfiles))
  gene_num_variants=matrix(0,length(gene_ids),length(genewise_ASEfiles))
  gene_minorCount=matrix(0,length(gene_ids),length(genewise_ASEfiles))
  gene_allCount=matrix(0,length(gene_ids),length(genewise_ASEfiles))
  gene.pval.binom=matrix(0,length(gene_ids),length(genewise_ASEfiles))
  gene_ase_ratio_weightedMean=matrix(0,length(gene_ids),length(genewise_ASEfiles))
  gene_ase_ratio_weightedvar=matrix(0,length(gene_ids),length(genewise_ASEfiles))
  
  for (i in (1:length(genewise_ASEfiles))) {
    gw_ase=read.table(genewise_ASEfiles[i],header=TRUE)
    im=match(gene_ids,gw_ase$gene_id)
    im_na=which(is.na(im)==TRUE)
    im_notna=which(is.na(im)==FALSE)
    
    var_snps[im_na,i]=NA
    for (hu in 1:length(gene_ids[im_notna])) {
      #var_snps[im_notna[hu],i]=var(gw_ase[which(gw_ase$gene_id==gene_ids[im_notna][hu]),"snp_ase_ratio"])
      var_snps[im_notna[hu],i]=sum((gw_ase[which(gw_ase$gene_id==gene_ids[im_notna][hu]),"snp_ase_ratio"])>0.1)/length((gw_ase[which(gw_ase$gene_id==gene_ids[im_notna][hu]),"snp_ase_ratio"]))
      #var_snps[im_notna[hu],i]=diff(range(gw_ase[which(gw_ase$gene_id==gene_ids[im_notna][hu]),"snp_ase_ratio"]))
    }
    for (hu in 1:length(gene_ids[im_notna])) {
      #var_snps[im_notna[hu],i]=var(gw_ase[which(gw_ase$gene_id==gene_ids[im_notna][hu]),"snp_ase_ratio"])
      var_snps_xist[im_notna[hu],i]=sum((gw_ase[which(gw_ase$gene_id==gene_ids[im_notna][hu]),"snp_ase_ratio"])<0.05)/length((gw_ase[which(gw_ase$gene_id==gene_ids[im_notna][hu]),"snp_ase_ratio"]))
      #var_snps[im_notna[hu],i]=diff(range(gw_ase[which(gw_ase$gene_id==gene_ids[im_notna][hu]),"snp_ase_ratio"]))
    }
    ase_ratios[im_na,i]=NA
    ase_ratios[im_notna,i]=gw_ase$gene_ase_ratio[im[im_notna]]
    gene_num_variants[im_na,i]=NA
    gene_num_variants[im_notna,i]=gw_ase$gene_num_variants[im[im_notna]]
    gene_minorCount[im_na,i]=NA
    gene_minorCount[im_notna,i]=gw_ase$gene_minorCount[im[im_notna]]
    gene_allCount[im_na,i]=NA
    gene_allCount[im_notna,i]=gw_ase$gene_allCount[im[im_notna]]
    gene.pval.binom[im_na,i]=NA
    gene.pval.binom[im_notna,i]=gw_ase$gene.pval.binom[im[im_notna]]
    gene_ase_ratio_weightedMean[im_na,i]=NA
    gene_ase_ratio_weightedMean[im_notna,i]=gw_ase$gene_ase_ratio_weightedMean[im[im_notna]]
    gene_ase_ratio_weightedvar[im_na,i]=NA
    gene_ase_ratio_weightedvar[im_notna,i]=gw_ase$gene_ase_ratio_weightedvar[im[im_notna]]
  }
  
  ttt=which((rowSums(!is.na(ase_ratios[,which(D$xist_group=="High-XIST female")])))>=min_nonna_num & (rowSums(!is.na(ase_ratios[,which(D$xist_group=="Low-XIST female")])))>=min_nonna_num)
  if (mychr=="X") {
    xist_ind=which(gene_ids==gene_info$gene_id[match("XIST",gene_info$gene_name)])
    if (!(xist_ind %in% ttt)) {
      ttt=append(ttt,xist_ind,after=(which(ttt>xist_ind)[1]-1))
    }
  }
  print(length(ttt)) # Number of genes after min_nonna_num filtering: 139
  
  ase_ratios=ase_ratios[ttt,]
  gene_minorCount=gene_minorCount[ttt,]
  gene_allCount=gene_allCount[ttt,]
  gene_num_variants=gene_num_variants[ttt,]
  var_snps=var_snps[ttt,]
  var_snps_xist=var_snps_xist[ttt,]
  gene_ids=gene_ids[ttt]
  gene_info=gene_info[ttt,]
  
  n_genes=dim(ase_ratios)[1]
  n_samples=dim(ase_ratios)[2]
  
  low_l=matrix(NA,n_genes,n_samples)
  up_l=matrix(NA,n_genes,n_samples)
  p_est=matrix(NA,n_genes,n_samples)
  pval=matrix(NA,n_genes,n_samples)
  for (i in 1:n_genes) {
    for (j in 1:n_samples) {
      if (!is.na(gene_minorCount[i,j])) {
        p_est[i,j]=(gene_minorCount[i,j]/gene_allCount[i,j])
        low_l[i,j]=p_est[i,j]-sqrt(p_est[i,j]*(1-p_est[i,j])/gene_allCount[i,j])*qnorm(1-0.05/2)
        up_l[i,j]=p_est[i,j]+sqrt(p_est[i,j]*(1-p_est[i,j])/gene_allCount[i,j])*qnorm(1-0.05/2)
        pval[i,j]=binom.test(gene_minorCount[i,j],gene_allCount[i,j],p=min_ase_ratio_for_escape,alternative=altern_hypt)$p.value
      }
    }
  }
  
  qval=qvalue(pval)$qvalues
  sig=(qval<0.01)
  
  colnames(ase_ratios)=commonfiles
  rownames(ase_ratios)=gene_ids
  
  colnames(sig)=commonfiles
  rownames(sig)=gene_ids
  
  res=list(ase_ratios=ase_ratios,sig=sig,D=D,gene_info=gene_info,pval=pval,
           gene_minorCount=gene_minorCount,gene_allCount=gene_allCount,gene_num_variants=gene_num_variants)
  return(res)
}
  