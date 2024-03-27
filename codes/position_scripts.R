# -----------------------------------------------------------
# Script Name: position_scripts.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

library(rdist)

find_overlapping_genes <- function(biomart_object,gene="XIST") {
  
  target_gene=biomart_object[which(biomart_object$external_gene_id==gene),]
  chr=target_gene$chromosome_name
  t_s=target_gene$start_position
  t_e=target_gene$end_position
  
  s=biomart_object$start_position[biomart_object$chromosome_name==chr]
  e=biomart_object$end_position[biomart_object$chromosome_name==chr]
  genes=biomart_object$external_gene_id[biomart_object$chromosome_name==chr]

 overlapping_gene_ids=which((s>=t_s & e<=t_e) | (s<=t_s & e>t_s & e<=t_e) | (s>=t_s & s<=t_e & e>t_e) | (s<=t_s & e>=t_e))
 
 if (length(overlapping_gene_ids)>0) {
   overlapping_genes=genes[overlapping_gene_ids]
 } else {
   overlapping_genes="none"
 }
 
 if (gene %in% overlapping_genes) {overlapping_genes=overlapping_genes[overlapping_genes!=gene]}
 
 return(overlapping_genes)
 
}

find_overlapping_positions <- function(biomart_object,gene1="XIST",gene2="TSIX") {
  
  gene1=biomart_object[which(biomart_object$external_gene_id==gene1),]
  gene2=biomart_object[which(biomart_object$external_gene_id==gene2),]
  
  if (gene1$chromosome_name != gene2$chromosome_name) {
    stop("Error: Genes are not on the same chromosome!")
  }
  
  s_1=gene1$start_position
  e_1=gene1$end_position
  s_2=gene2$start_position
  e_2=gene2$end_position
  
  if (s_2>=s_1 & e_2<=e_1) {
    s_intersect=s_2 
    e_intersect=e_2
  } else if (s_2<=s_1 & e_2>=s_1 & e_2<=e_1) {
    s_intersect=s_1
    e_intersect=e_2
  } else if (s_2>=s_1 & s_2<=e_1 & e_2>=e_1) {
    s_intersect=s_2 
    e_intersect=e_1
  } else if (s_2<=s_1 & e_2>=e_1) {
    s_intersect=s_1
    e_intersect=e_1
  } else {
    s_intersect=NA
    e_intersect=NA
  }

  overlapping_region=c(s_intersect, e_intersect)
  names(overlapping_region)=c("start_position", "end_position")
  
  return(overlapping_region)

}


find_positions <- function(biomart_object,chr=NULL,select_gene=NULL,omit_gene=NULL,ensembl=TRUE) {
  
  if (!is.null(select_gene)) {
    if (ensembl) {
      select_gene=biomart_object[which(biomart_object$ensembl_gene_id==select_gene),]
    } else {
      select_gene=biomart_object[which(biomart_object$external_gene_id==select_gene),]
    }
    r=c(select_gene$start_position,select_gene$end_position)
    names(r)=c("start_position", "end_position")
  }
  if (!is.null(chr)) {
    chr=biomart_object[which(biomart_object$chromosome_name==chr),]
    r=cbind(chr$start_position,chr$end_position)
    colnames(r)=c("start_position", "end_position")
    if (ensembl) {
      rownames(r)=chr$ensembl_gene_id
    } else {
      rownames(r)=chr$external_gene_id
    }
    if (!is.null(omit_gene)) {
      r=r[-which(rownames(r)==omit_gene),]
    }
  }
  
  return(r)

}


get_variant_info <- function(ase, chr, select_position_vector=NULL, omit_position_vector=NULL, min_total_count=20, min_read_length=NULL) {
  
  ase=ase[which(ase$contig==chr & ase$totalCount>=min_total_count),]
  if (!is.null(select_position_vector)) {
    r=ase[which(ase$position>=select_position_vector[1] & ase$position<=select_position_vector[2]),]
  } else if (!is.null(omit_position_vector)) {
    r=ase[-which(ase$position>=omit_position_vector[1] & ase$position<=omit_position_vector[2]),]
  } else {
    r=ase
  }
  
  
  if (((dim(r)[1]>0) & (!is.null(min_read_length)))) {
    p=pdist(r$position)
    p[lower.tri(p,diag=TRUE)]=NA
    w=which(p<min_read_length,arr.ind=TRUE)
    if (dim(w)[1]>0) {
      exc=c()
      for (z in (1:(dim(w)[1]))) {
        m=r$totalCount[w[z,]]
        exc=append(exc,w[z,m<max(m)])
      }
      exc=unique(exc)
      r=r[-exc,]
    }
  }
  
  return(r)
  
}

get_ase_ratio <- function(variant_info,overlapping_info=NULL) {
  
  # if overlapping_info != NULL, omit overlapping positions
  if (!is.null(overlapping_info)) {
    variant_info=variant_info[!(variant_info$position %in% overlapping_info$position),]
  }
  snp_ratios=pmin(variant_info$refCount,variant_info$altCount)/(variant_info$refCount+variant_info$altCount)
  r=mean(pmin(variant_info$refCount,variant_info$altCount)/(variant_info$refCount+variant_info$altCount))
  num_variants=dim(variant_info)[1]
  res=list(ase_ratio=r,num_variants=num_variants,variant_info=variant_info,snp_ratios=snp_ratios)
  return(res)
  
}

get_ase_ratio2 <- function(variant_info,overlapping_info=NULL,alpha=4,beta=4,test.value=0.1) {
  
  # if overlapping_info != NULL, omit overlapping positions
  if (!is.null(overlapping_info)) { 
    variant_info=variant_info[!(variant_info$position %in% overlapping_info$position),]
  }
  num_variants=dim(variant_info)[1]
  AA=variant_info
  y_=pmin(AA$refCount,AA$altCount)
  n_=AA$totalCount
  Ep=(alpha+y_)/(alpha+beta+n_)
  Vp=((alpha+y_)*(beta+n_-y_))/(((alpha+beta+n_)^2)*(alpha+beta+n_+1))
  alpha_post=alpha+y_
  beta_post=beta+n_-y_
  weighted_mean=sum(Ep/Vp)/sum(1/Vp)
  weighted_var=1/sum(1/Vp)
  #loglik_bial=dbeta(0.5, alpha_post, beta_post, ncp = 0, log = TRUE)
  #loglik_monoal=dbeta(0, alpha_post, beta_post, ncp = 0, log = TRUE)
  snp_ratios=pmin(variant_info$refCount,variant_info$altCount)/(variant_info$refCount+variant_info$altCount)
  r1=mean(pmin(variant_info$refCount,variant_info$altCount)/(variant_info$refCount+variant_info$altCount)) # mean of snp ratios
  r2=sum(pmin(variant_info$refCount,variant_info$altCount))/sum(variant_info$refCount+variant_info$altCount) # sum over snps
  if (num_variants>0) {
    gene_minorAllele=sum(pmin(variant_info$refCount,variant_info$altCount))
    gene_allAllele=sum(variant_info$refCount+variant_info$altCount)
    p.val=binom.test(sum(pmin(variant_info$refCount,variant_info$altCount)), sum(variant_info$refCount+variant_info$altCount), test.value, alternative = "less")$p.value 
    #res=list(ase_ratio_gene_sumsnps=r2,p.value.binom=p.val,ase_ratio_gene_meansnps=r1,ase_ratio_gene_weightedMean=weighted_mean,ase_ratio_gene_weightedVar=weighted_var,
    #         t=NA,p.value.t=NA,num_variants=num_variants,variant_info=variant_info,snp_ratios=snp_ratios,snp_ratios_postMean=Ep,snp_ratios_postVar=Vp,alpha_post=alpha_post,beta_post=beta_post)
    #if (num_variants>1) {
      t=(weighted_mean-test.value)/sqrt(weighted_var)
      p.val.t=pnorm(t, mean = 0, sd = 1, lower.tail = TRUE) 
      #p.val.t=pt(q=t, df=num_variants-1, lower.tail=TRUE)
      res=list(ase_ratio_gene_sumsnps=r2,p.value.binom=p.val,ase_ratio_gene_meansnps=r1,ase_ratio_gene_weightedMean=weighted_mean,ase_ratio_gene_weightedVar=weighted_var,
               t=t,p.value.t=p.val.t,num_variants=num_variants,variant_info=variant_info,snp_ratios=snp_ratios,snp_ratios_postMean=Ep,snp_ratios_postVar=Vp,alpha_post=alpha_post,beta_post=beta_post,gene_minorCount=gene_minorAllele,gene_allCount=gene_allAllele)
    #}
  } else {
    res=list(ase_ratio_gene_sumsnps=r2,ase_ratio_gene_meansnps=r1,ase_ratio_gene_weightedMean=weighted_mean,ase_ratio_gene_weightedVar=weighted_var,
             num_variants=num_variants,variant_info=variant_info,snp_ratios=snp_ratios,snp_ratios_postMean=Ep,snp_ratios_postVar=Vp,alpha_post=alpha_post,beta_post=beta_post)
  }
  return(res)
}

get_major_allele <- function(variant_info) {
  major_allele=apply(as.matrix(variant_info[,which(colnames(variant_info)=="refAllele"):(which(colnames(variant_info)=="refAllele")+3)]), 1, function(x) x[as.integer(as.numeric(x["refCount"])<as.numeric(x["altCount"]))+1])
  return(as.vector(major_allele))
}

get_minor_allele <- function(variant_info) {
  minor_allele=apply(as.matrix(variant_info[,which(colnames(variant_info)=="refAllele"):(which(colnames(variant_info)=="refAllele")+3)]), 1, function(x) x[as.integer(as.numeric(x["refCount"])>as.numeric(x["altCount"]))+1])
  return(as.vector(minor_allele))
}

get_ref_allele <- function(variant_info) {
  ref_allele=apply(as.matrix(variant_info[,which(colnames(variant_info)=="refAllele"):(which(colnames(variant_info)=="refAllele")+3)]), 1, function(x) x["refAllele"])
  return(as.vector(ref_allele))
}

get_alt_allele <- function(variant_info) {
  alt_allele=apply(as.matrix(variant_info[,which(colnames(variant_info)=="refAllele"):(which(colnames(variant_info)=="refAllele")+3)]), 1, function(x) x["altAllele"])
  return(as.vector(alt_allele))
}

name2id <- function(biomart_object,name) {

  return(biomart_object$ensembl_gene_id[match(name,biomart_object$external_gene_id)])
  #return(biomart_object$ensembl_gene_id[which(biomart_object$external_gene_id==name)])
  
}
    
id2name <- function(biomart_object,id) {
  
  return(biomart_object$external_gene_id[match(id,biomart_object$ensembl_gene_id)])
  #return(biomart_object$external_gene_id[which(biomart_object$ensembl_gene_id==id)])
  
}