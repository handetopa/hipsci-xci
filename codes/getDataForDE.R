# -----------------------------------------------------------
# Script Name: getDataForDE.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

getDataForDE <- function(path) {
  
  hipsci_datadir = file.path(path,"data")
  countsData=file.path(hipsci_datadir,"hipsci_rnaseq_byLineCenter.RData")
  load(countsData)
 
  counts_simple=counts_combined
  D_simple=hipsci_info_combined
  gene_info=gene_info_orig
  
  library(edgeR)
  
  #male_donors=unique(D_simple$donor[D_simple$sex=="male"]) ## 77 male donors
  #female_donors=unique(D_simple$donor[D_simple$sex=="female"]) ## 110 female donors
  unique_donors=unique(D_simple$donor)
  D_simple$number_of_lines_from_same_donor=NA
  for (i in 1:length(unique_donors)) {
    D_simple$number_of_lines_from_same_donor[which(D_simple$donor==unique_donors[i])]=length(unique(D_simple$lines[which(D_simple$donor==unique_donors[i])]))
  }
  
  D_simple$number_of_centers_sequenced_the_line=table(D_simple$lines)[match(D_simple$lines,names(table(D_simple$lines)))]
  
  xist_ind=which(rownames(counts_simple)=="ENSG00000229807.10")
  hist(counts_simple[xist_ind,which(D_simple$sex=="male")])
  hist(counts_simple[xist_ind,which(D_simple$sex=="female")])
  
  groups=D_simple$sex
  keep=filterByExpr(counts_simple,group=groups,min.count=10)
  keep[xist_ind]=TRUE
  counts_simple=counts_simple[keep,]
  gene_info=gene_info[keep,]
  d0=DGEList(counts_simple)
  d0=calcNormFactors(d0,method="TMM")
  cpm_log=cpm(d0,log=TRUE,prior.count = 0.5)
  
  D_simple$xist_cpm_log=cpm_log[which(rownames(cpm_log)=="ENSG00000229807.10"),]
  
  library(reshape)
  library(ggplot2)
  # Check the expression of pluripotency marker genes:
  nanog_ind=which(gene_info$gene_name=="NANOG") 
  sox2_ind=which(gene_info$gene_name=="SOX2") 
  oct4_ind=which(gene_info$gene_name=="POU5F1") 
  
  df=data.frame(NANOG=cpm_log[nanog_ind,],SOX2=cpm_log[sox2_ind,],PO5F1=cpm_log[oct4_ind,],Sex=D_simple$sex,lines=D_simple$lines)
  df1=melt(df)
  
  lab1=sub("_[^_]+$", "", colnames(cpm_log)[which(cpm_log[nanog_ind,]<0)])
  lab2=sub("_[^_]+$", "", colnames(cpm_log)[which(cpm_log[sox2_ind,]<0)])
  lab3=sub("_[^_]+$", "", colnames(cpm_log)[which(cpm_log[oct4_ind,]<0)])
  
  pdf(file.path(path,"figures/pluripotency_genes.pdf"),width=8, height=8)
  par(oma = c(0,0,0,0),ps=20,cex=2)
  p=ggplot(df1, aes(x=variable, y=value,color=Sex)) + 
    geom_boxplot(outlier.shape=NA) +
    geom_point(position = position_jitterdodge(seed = 42),size=0.6) +
    scale_fill_discrete(labels=c("Female", "Male")) +
    ylab("Expression (log-cpm)") +
    xlab("Pluripotency marker genes") + 
    theme(text = element_text(size = 20)) +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20))
  print(p)
  dev.off()
  
  D_simple$xist_group=NA
  D_simple$xist_group[which(D_simple$xist_cpm_log<1.5 & D_simple$sex=="female")]="lowXIST"
  D_simple$xist_group[which(D_simple$xist_cpm_log>=1.5 & D_simple$sex=="female")]="highXIST"
  D_simple$xist_group[which(D_simple$sex=="male")]="male"
  
  discordant_donors=intersect(D_simple$donor[which(D_simple$xist_group=="lowXIST")],D_simple$donor[which(D_simple$xist_group=="highXIST")])
  discordant_lines=D_simple$lines[which(D_simple$donor %in% discordant_donors)]
  D_simple$concordant=TRUE
  D_simple$concordant[which(D_simple$lines %in% discordant_lines)]=FALSE
  
  D_simple$type="ipsc"
  #D_simple$center_name[which(D_simple$center_name=="WELLCOME TRUST SANGER INSTITUTE")]=gsub(" ", "_", "WELLCOME TRUST SANGER INSTITUTE")
  save(D_simple,cpm_log,gene_info,d0,file=file.path(hipsci_datadir,"data_for_DE_new.RData"))
  
}