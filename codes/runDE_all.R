# -----------------------------------------------------------
# Script Name: runDE_all.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

runDE_all <- function(path) {
  
  source(file.path(path,"codes/runDE.R"))
  
  # ipscs
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="MH",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=FALSE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="ML",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=FALSE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="HL",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=FALSE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="MHL",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=FALSE)
  
  # ipscs with clusters
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="hl",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="mh",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="mhe",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="ml",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="hhe",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="lhe",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="lheh",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="lhhe",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="mlhe",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="mhhe",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="mlheh",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_new.RData"),comparison="mlh",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  
  # diff neurons
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_difneuron.RData"),comparison="MH",celltype="neuron",treatasind=TRUE,usesva=TRUE,by_cluster=FALSE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_difneuron.RData"),comparison="ML",celltype="neuron",treatasind=TRUE,usesva=TRUE,by_cluster=FALSE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_difneuron.RData"),comparison="HL",celltype="neuron",treatasind=TRUE,usesva=TRUE,by_cluster=FALSE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_difneuron.RData"),comparison="MHL",celltype="neuron",treatasind=TRUE,usesva=TRUE,by_cluster=FALSE)

}