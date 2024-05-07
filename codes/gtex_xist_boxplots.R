# -----------------------------------------------------------
# Script Name: gtex_xist_boxplots.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

# Download GTEx Analysis v8 RNA-seq gene read counts by tissue (https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression)
# and place them in files folder. Also, download GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt from Metadata in GTEx portal.
gtex_counts_files=list.files(path="codes/gtex_gene_counts/files")
setwd("codes/gtex_gene_counts/files")
pp=read.table("../GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",header=TRUE,sep="\t")

library(edgeR)

i=1
print(i)
dd=read.table(gtex_counts_files[i],skip=2,header=TRUE)
sex=pp$SEX[match(paste("GTEX-",sub('.*GTEX.','',sub( "(^[^.]+[.][^.]+)(.+$)", "\\1", colnames(dd)[-c(1,2,3)])),sep=""),pp$SUBJID)]
counts_simple=dd[,-c(1,2,3)]
rownames(counts_simple)=dd$Name
groups=sex
keep=filterByExpr(counts_simple,group=groups,min.count=10) 
keep[which(rownames(counts_simple)=="ENSG00000229807.10")]=TRUE
counts_simple=counts_simple[keep,]
d0=DGEList(counts_simple)
d0=calcNormFactors(d0,method="TMM")
cpm_log=cpm(d0,log=TRUE,prior.count = 0.5)
xist_cpm_log=cpm_log[which(rownames(cpm_log)=="ENSG00000229807.10"),]

df=data.frame(xist=xist_cpm_log,sex=sex,tissue=as.character(i))
df$sex=as.factor(df$sex)
DF=df[which(df$sex=="2"),]


for (i in 2:length(gtex_counts_files)) {
  print(i)
  dd=read.table(gtex_counts_files[i],skip=2,header=TRUE)
  sex=pp$SEX[match(paste("GTEX-",sub('.*GTEX.','',sub( "(^[^.]+[.][^.]+)(.+$)", "\\1", colnames(dd)[-c(1,2,3)])),sep=""),pp$SUBJID)]
  counts_simple=dd[,-c(1,2,3)]
  rownames(counts_simple)=dd$Name
  groups=sex
  keep=filterByExpr(counts_simple,group=groups,min.count=10) # All is kept in ipscs, not in neurons.
  keep[which(rownames(counts_simple)=="ENSG00000229807.10")]=TRUE
  counts_simple=counts_simple[keep,]
  d0=DGEList(counts_simple)
  d0=calcNormFactors(d0,method="TMM")
  cpm_log=cpm(d0,log=TRUE,prior.count = 0.5)
  xist_cpm_log=cpm_log[which(rownames(cpm_log)=="ENSG00000229807.10"),]
  
  df=data.frame(xist=xist_cpm_log,sex=sex,tissue=as.character(i))
  df$sex=as.factor(df$sex)
  
  DF=rbind(DF,df[which(df$sex=="2"),])
  
}

DF$tissue2=factor(DF$tissue, labels=sub('.*\\v8_', '', gsub("\\..*", "", gtex_counts_files[1:54])), levels=as.character(1:54))
DF$tissue2=paste("GTEx_",DF$tissue2,sep="")
load("../../../data/data_for_DE_new.RData")
df_ipsc=data.frame(xist=D_simple$xist_cpm_log[which(D_simple$sex=="female")],sex=rep("2",sum(D_simple$sex=="female")),tissue=rep("HipSci_iPSC",sum(D_simple$sex=="female")),tissue2=rep("HipSci_iPSC",sum(D_simple$sex=="female")))
DF=rbind(DF,df_ipsc)

library(ggplot2)
DF$tissue2=as.factor(DF$tissue2)
pp=ggplot(DF,aes(x=tissue2,y=xist,fill=tissue2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="gray34", size=0.4, alpha=0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "none") +
  ylab(expression("XIST expression (log"[2]*"(CPM+0.5))")) +
  xlab("") +
  geom_hline(yintercept=1.5,linetype=2,color="gray8") +
  theme(text = element_text(size=20)) 
ggsave("../../../figures/gtex_xist_boxplots.pdf", pp,width=50,height=25,units="cm",limitsize = FALSE)


load("../../../data/data_for_DE_new.RData")
dim(D_simple)
df=data.frame(xist=D_simple$xist_cpm_log,xist_group=D_simple$xist_group)
df$xist_group=as.factor(df$xist_group)

pp=ggplot(df,aes(x=xist_group,y=xist,color=xist_group)) +
  geom_boxplot(outlier.shape = NA,size=1) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_minimal() +
  scale_color_manual(name = "Female\ncell lines", breaks=c("highXIST","lowXIST","male"),values=c("#CC0033","#FF99CC","#3399FF")) +
  scale_x_discrete(breaks=c("highXIST","lowXIST","male"),labels=c("high-XIST\nfemales","low-XIST\nfemales","males")) +
  theme(legend.position = "none") +
  ylab(expression("XIST expression (log"[2]*"(CPM+0.5))")) +
  xlab("") +
  theme(text = element_text(size=20)) 
ggsave("../../../figures/xist_levels_by_groups.pdf", pp, width=15,height=15,units="cm",limitsize = FALSE)
