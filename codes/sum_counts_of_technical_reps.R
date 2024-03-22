sum_counts_of_technical_reps <- function(counts,hipsci_info,bycenter=TRUE) {
  
  sum_reps <- function(H,C,n) {
    # H:hipsci_info, C:counts, n:n_genes
    unique_lines=unique(H$lines)
    n_unique_lines=length(unique_lines)
    counts_summed=matrix(0,n,n_unique_lines)
    for (i in (1:n_unique_lines)) {
      accs_num=H$run_accession[which(H$lines==unique_lines[i])]
      inds_cellline=which(colnames(C) %in% accs_num)
      counts_summed[,i]=rowSums(as.matrix(C[,inds_cellline]))
    }
    colnames(counts_summed)=unique_lines
    rownames(counts_summed)=rownames(C)
    return(counts_summed)
  }
  
  n_genes=dim(counts)[1]
  if (bycenter==TRUE) {
    center_names=unique(hipsci_info$center_name)
    n_center_names=length(center_names)
    counts_new=list(NA,n_center_names)
    names(counts_new)=center_names
    for (i in (1:n_center_names)) {
      CC=counts[,which(hipsci_info$center_name==center_names[i])]
      HH=hipsci_info[which(hipsci_info$center_name==center_names[i]),]
      counts_new[[i]]=sum_reps(HH,CC,n_genes)
    }
  } else {
    counts_new=sum_reps(hipsci_info,counts,n_genes)
  }
  
  return(counts_new)
}


