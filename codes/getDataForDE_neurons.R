# -----------------------------------------------------------
# Script Name: getDataForDE_neurons.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

getDataForDE_neurons <- function(path) {
  
  source(file.path(path,"codes/sum_counts_of_technical_reps.R"))
  hipsci_datadir = file.path(path,"data")
  countsData=file.path(hipsci_datadir,"counts_neurons.RData")
  load(countsData)
  report=read.table(file.path(hipsci_datadir,"filereport_read_run_PRJEB18630_tsv.txt"),header=TRUE,sep="\t")
  report=report[which(report$library_strategy=="RNA-Seq"),]
  
  report$lines=paste("H",gsub("\\,.*","",gsub(".*H","",report$sample_title)),sep="")
  report$donor=gsub("\\_.*", "", report$lines)
  
  counts_new=sum_counts_of_technical_reps(counts,report,bycenter=TRUE)
  counts_neurons=counts_new$SC
  
  hipsci_info_neurons=report[which((report$lines %in% colnames(counts_new$SC)) & (report$center_name=="SC")),]
  hipsci_info_neurons=hipsci_info_neurons[!duplicated(hipsci_info_neurons$lines),]
  hipsci_info_neurons=hipsci_info_neurons[match(colnames(counts_new$SC),hipsci_info_neurons$lines),]
  colnames(counts_neurons)=paste(colnames(counts_neurons),"_neuron",sep="")

  hipsci_info_neurons=hipsci_info_neurons[,-which(names(hipsci_info_neurons) %in% c("experiment_accession","run_accession","read_count","base_count","sample_title","study_title","fastq_ftp"))]
  
  #dim(hipsci_info_neurons)==length(unique(hipsci_info_neurons$donor))
  
  load(file.path(hipsci_datadir,"data_for_DE_ipsc.RData"))
  D_simple_ipsc=D_simple
  
  male_donors=unique(D_simple$donor[D_simple$sex=="male"]) 
  female_donors=unique(D_simple$donor[D_simple$sex=="female"]) 
  unique_donors=c(male_donors,female_donors)
  
  counts_simple=d0$counts
  inds=which(colnames(D_simple) %in% c("lines", "disease_status", "sex", "age", "passage_number","ethnicity", "center_name", "predicted_population", "donor", "number_of_lines_from_same_donor"))
  D_simple=D_simple[,inds]
  
  rm.ind=which(!(hipsci_info_neurons$donor %in% D_simple$donor))
  hipsci_info_neurons$donor[rm.ind]
  colnames(counts_neurons)[rm.ind]
  
  hipsci_info_neurons=hipsci_info_neurons[-rm.ind,]
  counts_neurons=counts_neurons[,-rm.ind]
  counts_neurons=counts_neurons[match(gene_info$gene_id,rownames(counts_neurons)),]
  hipsci_info_neurons$number_of_lines_from_same_donor=1
  
  hipsci_info_neurons$disease_status=D_simple$disease_status[match(hipsci_info_neurons$donor,D_simple$donor)]
  hipsci_info_neurons$sex=D_simple$sex[match(hipsci_info_neurons$donor,D_simple$donor)]
  hipsci_info_neurons$age=D_simple$age[match(hipsci_info_neurons$donor,D_simple$donor)]
  hipsci_info_neurons$ethnicity=D_simple$ethnicity[match(hipsci_info_neurons$donor,D_simple$donor)]
  hipsci_info_neurons$predicted_population=D_simple$predicted_population[match(hipsci_info_neurons$donor,D_simple$donor)]
  hipsci_info_neurons$passage_number=D_simple$passage_number[match(hipsci_info_neurons$lines,D_simple$lines)]
  hipsci_info_neurons=hipsci_info_neurons[,match(names(D_simple),names(hipsci_info_neurons))]
  hipsci_info_neurons$lines=as.character(hipsci_info_neurons$lines)
  hipsci_info_neurons$donor=as.character(hipsci_info_neurons$donor)
  donor_of_neurons=unique(hipsci_info_neurons$donor)
  
  D_simple=hipsci_info_neurons
  counts_simple=counts_neurons
  groups=D_simple$sex
  xist_ind=which(rownames(counts_simple)=="ENSG00000229807.10")
  keep=filterByExpr(counts_simple,group=groups,min.count=10)
  keep[xist_ind]=TRUE
  counts_simple=counts_simple[keep,]
  gene_info=gene_info[keep,]
  d0=DGEList(counts_simple)
  d0=calcNormFactors(d0,method="TMM")
  cpm_log=cpm(d0,log=TRUE,prior.count = 0.5)
  
  D_simple$xist_cpm_log=cpm_log[which(rownames(cpm_log)=="ENSG00000229807.10"),]
  D_simple$xist_group=NA  # neuron_xist_group
  D_simple$xist_group[which(D_simple$xist_cpm_log<1.5 & D_simple$sex=="female")]="lowXIST"
  D_simple$xist_group[which(D_simple$xist_cpm_log>=1.5 & D_simple$sex=="female")]="highXIST"
  D_simple$xist_group[which(D_simple$sex=="male")]="male"
  D_simple$ipsc_xist_group=D_simple_ipsc$xist_group[match(D_simple$lines,D_simple_ipsc$lines)]

  clusters=read.table(file.path(hipsci_datadir,"female_clusters.txt"),header=TRUE,sep="\t")
  D_simple$ipsc_xist_cluster=clusters$clusterName[match(D_simple$lines,clusters$lines)]
  D_simple$ipsc_xist_cluster[which(D_simple$sex=="male")]="male"
  D_simple$type="neuron"
  
  save(D_simple,cpm_log,gene_info,d0,file=file.path(hipsci_datadir,"data_for_DE_difneuron.RData"))

}
