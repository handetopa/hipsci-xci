get_hipsci_counts <- function(part=1,hipsci_info="/scratch/project_2004435/hipsci/rnaseq-pipeline/INFO/hipsci_info_3parts.Rdata",rnaseqcDir="/scratch/project_2004435/hipsci/rnaseq-pipeline/rnaseqc",countsDir="/scratch/project_200443\
5/hipsci/rnaseq-pipeline/INFO") {
  
  load(hipsci_info)
  if (part==1) {
    accession=d1$accession
  } else if (part==2) {
    accession=d2$accession
  } else if (part==3) {
    accession=d3$accession
  }
  
  n_samples=length(accession)
  i=1
  fileName=file.path(rnaseqcDir,paste(accession[i],".gene_reads.gct",sep=""))
  rnaseqc=read.table(fileName,skip=3)
  n_genes=dim(rnaseqc)[1]
  gene_id=rnaseqc$V1
  gene_name=rnaseqc$V2
  counts=matrix(0,n_genes,n_samples)
  counts[,i]=rnaseqc$V3
  
  for (i in 2:n_samples) {
    fileName=file.path(rnaseqcDir,paste(accession[i],".gene_reads.gct",sep=""))
    rnaseqc=read.table(fileName,skip=3)
    counts[,i]=rnaseqc$V3
  }
  
  colnames(counts)=accession
  rownames(counts)=gene_id
  countsFile=file.path(countsDir,paste("counts_part",part,".Rdata",sep=""))
  save(counts,file=countsFile)
  
  #  gene_info=data.frame(gene_id=gene_id,gene_name=gene_name)
  #  gene_infoFile=file.path(countsDir,"geneInfo.Rdata")
  #  save(gene_info,file=gene_infoFile)
  
}

