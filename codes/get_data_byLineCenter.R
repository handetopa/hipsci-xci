# -----------------------------------------------------------
# Script Name: get_data_byLineCenter.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

get_data_byLineCenter <- function(path) {
  
  source(file.path(path,"codes/sum_counts_of_technical_reps.R"))
  hipsci_datadir = file.path(path,"data")
  raw_data=file.path(hipsci_datadir,"hipsci_rnaseq_rawData.RData")
  load(raw_data)
  gene_info_orig=gene_info
  counts_new=sum_counts_of_technical_reps(counts,hipsci_info,bycenter=TRUE)
  
  hipsci_info_SC=hipsci_info[which((hipsci_info$lines %in% colnames(counts_new$SC)) & (hipsci_info$center_name=="SC")),]
  unique_lines=colnames(counts_new$SC)
  n_unique_lines=length(unique_lines)
  read_count_summed=matrix(0,n_unique_lines,1)
  base_count_summed=matrix(0,n_unique_lines,1)
  num_of_technical_reps=matrix(0,n_unique_lines,1)
  for (i in (1:n_unique_lines)) {
    num_of_technical_reps[i]=length(which(hipsci_info_SC$lines==unique_lines[i]))
    read_count_summed[i]=sum(hipsci_info_SC$read_count[which(hipsci_info_SC$lines==unique_lines[i])])
    base_count_summed[i]=sum(hipsci_info_SC$base_count[which(hipsci_info_SC$lines==unique_lines[i])])
  }
  hipsci_info_SC=hipsci_info_SC[!duplicated(hipsci_info_SC$lines),]
  hipsci_info_SC=hipsci_info_SC[match(colnames(counts_new$SC),hipsci_info_SC$lines),]
  hipsci_info_SC$read_count=read_count_summed
  hipsci_info_SC$base_count=base_count_summed
  hipsci_info_SC$num_of_technical_reps=num_of_technical_reps
  
  hipsci_info_WT=hipsci_info[which((hipsci_info$lines %in% colnames(counts_new$`WELLCOME TRUST SANGER INSTITUTE`)) & (hipsci_info$center_name=="WELLCOME TRUST SANGER INSTITUTE")),]
  unique_lines=colnames(counts_new$`WELLCOME TRUST SANGER INSTITUTE`)
  n_unique_lines=length(unique_lines)
  read_count_summed=matrix(0,n_unique_lines,1)
  base_count_summed=matrix(0,n_unique_lines,1)
  num_of_technical_reps=matrix(0,n_unique_lines,1)
  for (i in (1:n_unique_lines)) {
    num_of_technical_reps[i]=length(which(hipsci_info_WT$lines==unique_lines[i]))
    read_count_summed[i]=sum(hipsci_info_WT$read_count[which(hipsci_info_WT$lines==unique_lines[i])])
    base_count_summed[i]=sum(hipsci_info_WT$base_count[which(hipsci_info_WT$lines==unique_lines[i])])
  }
  hipsci_info_WT=hipsci_info_WT[!duplicated(hipsci_info_WT$lines),]
  hipsci_info_WT=hipsci_info_WT[match(colnames(counts_new$`WELLCOME TRUST SANGER INSTITUTE`),hipsci_info_WT$lines),]
  hipsci_info_WT$read_count=read_count_summed
  hipsci_info_WT$base_count=base_count_summed
  hipsci_info_WT$num_of_technical_reps=num_of_technical_reps
  
  counts_SC=counts_new$SC
  counts_WT=counts_new$`WELLCOME TRUST SANGER INSTITUTE`
  colnames(counts_WT)=paste(colnames(counts_WT),"_WT",sep="")
  colnames(counts_SC)=paste(colnames(counts_SC),"_SC",sep="")
  hipsci_info_combined=rbind(hipsci_info_SC,hipsci_info_WT)
  counts_combined=cbind(counts_SC,counts_WT)
  hipsci_info_combined$line_center=paste(hipsci_info_combined$lines,"_",hipsci_info_combined$center_name,sep="")
  hipsci_info_combined$donor_center=paste(hipsci_info_combined$donor,"_",hipsci_info_combined$center_name,sep="")
  
  save(hipsci_info_combined,counts_combined,gene_info_orig,file=file.path(hipsci_datadir,"hipsci_rnaseq_byLineCenter.RData"))
  
}