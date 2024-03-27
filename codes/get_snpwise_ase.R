# -----------------------------------------------------------
# Script Name: get_snpwise_ase.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

get_snpwise_ase <- function(ase_fileName1,ase_fileName2=NULL,gene=NULL,chr=NULL,write_file=TRUE) {
  
  library(rdist)
  load("/Users/topah/Desktop/hipsci_codes/data/gene_info/ref_hs37d5.RData")
  source("/Users/topah/Desktop/hipsci_codes/codes/position_scripts.R")
  source("/Users/topah/Desktop/hipsci_codes/codes/extract_text_before_first_dot.R")
  if (!file.exists(ase_fileName1)) {
    print(paste("No ASE file found for ",ase_fileName1,sep=""))
  } else {
    if (!is.null(ase_fileName2)) {
      # Needed for combining results from technical reps when using results with WASP-correction
      library(dplyr)
      ase_info1=read.table(ase_fileName1,header=TRUE)
      ase_info1$contig=gsub('chr', '', ase_info1$contig)
      ase_info2=read.table(ase_fileName2,header=TRUE)
      ase_info2$contig=gsub('chr', '', ase_info2$contig)
      ase_info=rbind(ase_info1,ase_info2)
      ase_info3 = ase_info %>% group_by(contig,position,refAllele,altAllele) %>% summarise(refCount=sum(refCount),altCount=sum(altCount),totalCount=sum(totalCount),
                                                                                         lowMAPQDepth=sum(lowMAPQDepth), lowBaseQDepth=sum(lowBaseQDepth), rawDepth=sum(rawDepth), otherBases=sum(otherBases), improperPairs=sum(improperPairs))
      ase_info=as.data.frame(ase_info3)
    } else {
      ase_info=read.table(ase_fileName1,header=TRUE)
    }
    if (!is.null(gene)) {
      j=which(ref$ensembl_gene_id==gene | ref$external_gene_id==gene)
    } else if (!is.null(chr)) {
      j=which(ref$chromosome_name==chr)
    } else {
      j=(1:dim(ref)[1])
    }
    gene_name=ref$external_gene_id[j]
    gene_id=ref$ensembl_gene_id[j]
    chromosome=ref$chromosome_name[j]
    start_position=ref$start_position[j]
    end_position=ref$end_position[j]
    gene_length=ref$end_position[j]-ref$start_position[j]+1
    gene_ase_info=apply(as.matrix(j), 1, function(x) get_ase_ratio2(get_variant_info(ase_info,chr=ref$chromosome_name[x],select_position_vector=unlist(ref[x,4:5]),min_total_count=20,min_read_length=100)))
    
    k=apply(as.matrix(1:length(j)), 1, function(x) gene_ase_info[[x]]$num_variants>0)
    gene_name=gene_name[k]
    gene_id=gene_id[k]
    gene_length=gene_length[k]
    start_position=start_position[k]
    end_position=end_position[k]
    chromosome=chromosome[k]
    gene_ase_ratio_meansnps=apply(as.matrix(which(k==TRUE)), 1, function(x) gene_ase_info[[x]]$ase_ratio_gene_meansnps)
    gene_ase_ratio_weightedMean=apply(as.matrix(which(k==TRUE)), 1, function(x) gene_ase_info[[x]]$ase_ratio_gene_weightedMean)
    gene_ase_ratio_weightedvar=apply(as.matrix(which(k==TRUE)), 1, function(x) gene_ase_info[[x]]$ase_ratio_gene_weightedVar)
    gene_ase_ratio_sumsnps=apply(as.matrix(which(k==TRUE)), 1, function(x) gene_ase_info[[x]]$ase_ratio_gene_sumsnps)
    gene_num_variants=apply(as.matrix(which(k==TRUE)), 1, function(x) gene_ase_info[[x]]$num_variants)
    pval.binom=apply(as.matrix(which(k==TRUE)), 1, function(x) gene_ase_info[[x]]$p.value.binom)
    pval.t=apply(as.matrix(which(k==TRUE)), 1, function(x) gene_ase_info[[x]]$p.value.t)
    gene_minorCount=apply(as.matrix(which(k==TRUE)), 1, function(x) gene_ase_info[[x]]$gene_minorCount)
    gene_allCount=apply(as.matrix(which(k==TRUE)), 1, function(x) gene_ase_info[[x]]$gene_allCount)
    #df=data.frame(gene_id=gene_id,gene_name=gene_name,chromosome=chromosome,gene_ase_ratio=gene_ase_ratio,gene_num_variants=gene_num_variants) # for genewise
    
    ## For snpwise ##
    jj=unlist(apply(as.matrix(which(k==TRUE)), 1, function(x) gene_ase_info[[x]]$snp_ratios))
    snp_position=unlist(apply(as.matrix(which(k==TRUE)), 1, function(x) gene_ase_info[[x]]$variant_info$position))
    major_allele=unlist(apply(as.matrix(which(k==TRUE)), 1, function(x) get_major_allele(gene_ase_info[[x]]$variant_info)))
    minor_allele=unlist(apply(as.matrix(which(k==TRUE)), 1, function(x) get_minor_allele(gene_ase_info[[x]]$variant_info)))
    ref_allele=unlist(apply(as.matrix(which(k==TRUE)), 1, function(x) get_ref_allele(gene_ase_info[[x]]$variant_info)))
    alt_allele=unlist(apply(as.matrix(which(k==TRUE)), 1, function(x) get_alt_allele(gene_ase_info[[x]]$variant_info)))
    snp_minorCount=unlist(apply(as.matrix(which(k==TRUE)), 1, function(x) pmin(gene_ase_info[[x]]$variant_info[,c("refCount","altCount")][,1],gene_ase_info[[x]]$variant_info[,c("refCount","altCount")][,2])))
    snp_allCount=unlist(apply(as.matrix(which(k==TRUE)), 1, function(x) as.matrix(rowSums(gene_ase_info[[x]]$variant_info[,c("refCount","altCount")]))))
    
    df=data.frame(gene_id=rep(gene_id,gene_num_variants),gene_name=rep(gene_name,gene_num_variants),chromosome=rep(chromosome,gene_num_variants),
                  snp_ase_ratio=jj,gene_ase_ratio=rep(gene_ase_ratio_sumsnps,gene_num_variants),
                  gene_num_variants=rep(gene_num_variants,gene_num_variants),start_position=rep(start_position,gene_num_variants),
                  end_position=rep(end_position,gene_num_variants),snp_position=snp_position,major_allele=major_allele,minor_allele=minor_allele,
                  ref_allele=ref_allele,alt_allele=alt_allele,
                  gene_minorCount=rep(gene_minorCount,gene_num_variants),gene_allCount=rep(gene_allCount,gene_num_variants),
                  snp_minorCount=snp_minorCount,snp_allCount=snp_allCount,
                  gene.pval.binom=rep(pval.binom,gene_num_variants),gene.pval.z=rep(pval.t,gene_num_variants),
                  gene_ase_ratio_meansnps=rep(gene_ase_ratio_meansnps,gene_num_variants),gene_ase_ratio_weightedMean=rep(gene_ase_ratio_weightedMean,gene_num_variants),
                  gene_ase_ratio_weightedvar=rep(gene_ase_ratio_weightedvar,gene_num_variants))
    
    if (!write_file) {
      return(df)
    } else {
      write.table(df,file=paste(extract_text_before_first_dot(ase_fileName1),"_min20_snpwiseASE_withGenePval.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
    }
  }
}




