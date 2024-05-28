# -----------------------------------------------------------
# Script Name: runDE_all.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

runDE_all <- function(path) {
  
  source(file.path(path,"codes/runDE.R"))
  
  # ipscs with 2 female XIST groups
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_ipsc.RData"),comparison="HL",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=FALSE)

  # ipscs with 3 female clusters
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_ipsc.RData"),comparison="G1G2",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_ipsc.RData"),comparison="MG1",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_ipsc.RData"),comparison="MG2",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_ipsc.RData"),comparison="MG3",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_ipsc.RData"),comparison="G1G3",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_ipsc.RData"),comparison="G3G2",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE)

  # diff neurons
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_difneuron.RData"),comparison="MH",celltype="neuron",treatasind=TRUE,usesva=TRUE,by_cluster=FALSE)
  runDE(path,dataForDE=file.path(path,"data/data_for_DE_difneuron.RData"),comparison="ML",celltype="neuron",treatasind=TRUE,usesva=TRUE,by_cluster=FALSE)
  
}