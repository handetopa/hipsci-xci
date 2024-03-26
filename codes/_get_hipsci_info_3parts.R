get_hipsci_info_3parts <- function() {
  
  lines=read.table("hipsci_lines.tsv",header=TRUE,sep="\t")
  files=read.delim("hipsci_files.tsv", header = TRUE, stringsAsFactors = FALSE, quote = "", sep = "\t")
  report=read.table("filereport_read_run_PRJEB7388_tsv.txt",header=TRUE,sep="\t")
  
  lines_orig=lines
  files_orig=files
  report_orig=report
  
  
  ind=match(lines$Name,files$Cell.line)
  accession=files$Accession[ind]
  passage_number=files$Cell.line.passage.number.for.this.assay[ind]
  
  ind2=match(accession,report$run_accession)
  fastqs=report$fastq_ftp[ind2]
  read_count=report$read_count[ind2]
  base_count=report$base_count[ind2]
  
  
  
  d1=data.frame(lines=lines$Name,disease_status=lines$Disease.Status,sex=lines$Sex,age=lines$Age,
                pluripot_score=lines$Pluritest.pluripotency.score,novelty_score=lines$Pluritest.novelty.score,
                CNV.num.different.regions=lines$CNV.num.different.regions,
                CNV.length.of.different.regions...Mbp=lines$CNV.length.of.different.regions...Mbp,culture=lines$Culture,
                passage_number=passage_number,read_count=read_count,base_count=base_count,accession=accession,fastq_ftp=fastqs)
  
  
  files=files[-ind,]
  report=report[-ind2,]
  ind=match(lines$Name,files$Cell.line)
  accession=files$Accession[ind]
  passage_number=files$Cell.line.passage.number.for.this.assay[ind]
  
  ind2=match(accession,report$run_accession)
  fastqs=report$fastq_ftp[ind2]
  read_count=report$read_count[ind2]
  base_count=report$base_count[ind2]
  
  d2=data.frame(lines=lines$Name,disease_status=lines$Disease.Status,sex=lines$Sex,age=lines$Age,
                pluripot_score=lines$Pluritest.pluripotency.score,novelty_score=lines$Pluritest.novelty.score,
                CNV.num.different.regions=lines$CNV.num.different.regions,
                CNV.length.of.different.regions...Mbp=lines$CNV.length.of.different.regions...Mbp,culture=lines$Culture,
                passage_number=passage_number,read_count=read_count,base_count=base_count,accession=accession,fastq_ftp=fastqs)
  
  
  files=files[-ind,]
  report=report[-ind2,]
  ind=match(lines$Name,files$Cell.line)
  accession=files$Accession[ind]
  passage_number=files$Cell.line.passage.number.for.this.assay[ind]
  
  ind2=match(accession,report$run_accession)
  fastqs=report$fastq_ftp[ind2]
  read_count=report$read_count[ind2]
  base_count=report$base_count[ind2]
  
  d3=data.frame(lines=lines$Name,disease_status=lines$Disease.Status,sex=lines$Sex,age=lines$Age,
                pluripot_score=lines$Pluritest.pluripotency.score,novelty_score=lines$Pluritest.novelty.score,
                CNV.num.different.regions=lines$CNV.num.different.regions,
                CNV.length.of.different.regions...Mbp=lines$CNV.length.of.different.regions...Mbp,culture=lines$Culture,
                passage_number=passage_number,read_count=read_count,base_count=base_count,accession=accession,fastq_ftp=fastqs)
  
  filter_in=which(d3$accession!="NA")
  d3=d3[filter_in,]
  
  save(d1,d2,d3,file="hipsci_info_3parts.Rdata")
  
}

###

