get_rawData <- function(path) {
  
  source(file.path(path,"codes/extract_text_before_first_dot.R"))
  hipsci_datadir = file.path(path,"data")
  load(file.path(hipsci_datadir,"counts.Rdata")) # counts
  
  lines=read.table(file.path(hipsci_datadir,"hipsci_lines.tsv"),header=TRUE,sep="\t")
  files=read.delim(file.path(hipsci_datadir,"hipsci_files.tsv"), header = TRUE, stringsAsFactors = FALSE, quote = "", sep = "\t")
  report=read.table(file.path(hipsci_datadir,"filereport_read_run_PRJEB7388_tsv.txt"),header=TRUE,sep="\t")
  dim(files)
  dim(report)
  ind=match(colnames(counts),files$Accession)
  files=files[ind,]
  ind=match(colnames(counts),report$run_accession)
  report=report[ind,]
  
  lines_ind=match(files$Cell.line,lines$Name)
  
  hipsci_info=data.frame(lines=files$Cell.line,run_accession=files$Accession, accession_type=files$Accession.type,
                         disease_status=files$Disease.status,sex=files$Sex,age=lines$Age[lines_ind],ethnicity=lines$Ethnicity[lines_ind],source_material=lines$Source.Material[lines_ind],
                         culture=files$Cell.line.growing.conditions.for.this.assay,cell_type=files$Cell.type,instrument=files$Instrument,assay=files$Assay,
                         passage_number=files$Cell.line.passage.number.for.this.assay, center_name=report$center_name,read_count=report$read_count,
                         tissue_provider=lines$Tissue.Provider[lines_ind],predicted_population=lines$Predicted.Population[lines_ind],method_of_derivation=lines$Method.of.derivation[lines_ind],
                         base_count=report$base_count, pluripot_score=lines$Pluritest.pluripotency.score[lines_ind],novelty_score=lines$Pluritest.novelty.score[lines_ind],
                         CNV.num.different.regions=lines$CNV.num.different.regions[lines_ind], CNV.length.of.different.regions...Mbp=lines$CNV.length.of.different.regions...Mbp[lines_ind],
                         fastq_ftp=report$fastq_ftp)
  
  hipsci_info$donor=gsub("\\_.*", "", hipsci_info$lines)
  
  # Remove feeder-dependent:
  counts=counts[,which(hipsci_info$culture=="Feeder-free")]
  hipsci_info=hipsci_info[which(hipsci_info$culture=="Feeder-free"),]
  
  table(hipsci_info$disease_status)
  table(hipsci_info$sex)
  table(hipsci_info$ethnicity)
  table(hipsci_info$source_material)
  table(hipsci_info$culture)
  table(hipsci_info$cell_type)
  table(hipsci_info$instrument)
  table(hipsci_info$assay)
  table(hipsci_info$center_name)
  table(hipsci_info$tissue_provider)
  table(hipsci_info$predicted_population)
  table(hipsci_info$method_of_derivation)
  
  table(hipsci_info$sex[which(hipsci_info$center_name=="SC")])
  table(hipsci_info$sex[which(hipsci_info$center_name=="WELLCOME TRUST SANGER INSTITUTE")])
  
  table(hipsci_info$age[which(hipsci_info$center_name=="SC")])
  table(hipsci_info$age[which(hipsci_info$center_name=="WELLCOME TRUST SANGER INSTITUTE")])
  
  table(hipsci_info$ethnicity[which(hipsci_info$center_name=="SC")])
  table(hipsci_info$ethnicity[which(hipsci_info$center_name=="WELLCOME TRUST SANGER INSTITUTE")])
  
  table(hipsci_info$method_of_derivation[which(hipsci_info$center_name=="SC")])
  table(hipsci_info$method_of_derivation[which(hipsci_info$center_name=="WELLCOME TRUST SANGER INSTITUTE")])
  
  hist(hipsci_info$passage_number[which(hipsci_info$center_name=="SC")])
  hist(hipsci_info$passage_number[which(hipsci_info$center_name=="WELLCOME TRUST SANGER INSTITUTE")])
  
  gtex_v8_info=read.table(file.path(hipsci_datadir,"gencode_genes.txt"),header=FALSE)
  ii=match(rownames(counts),gtex_v8_info[,2])
  gene_info=data.frame(gene_id=rownames(counts),gene_name=gtex_v8_info[ii,4],chr=gtex_v8_info[ii,1],gene_biotype=gtex_v8_info[ii,3])
  gene_info$chr=gsub('chr','',gene_info$chr)
  gene_info$chr[which(gene_info$chr=="M")]="MT"
  
  ind_mt_genes=which(gene_info$chr=="MT")
  
  # Remove PAR_Y genes:
  ind_par_y=which(grepl("\\PAR_Y$", gene_info$gene_id)==TRUE)
  par_genes=gsub("\\_.*", "", gene_info$gene_id[ind_par_y])
  sum(counts[ind_par_y,])
  counts=counts[-ind_par_y,]
  gene_info=gene_info[-ind_par_y,]
  
  sum(colnames(counts)==hipsci_info$run_accession)==dim(counts)[2]
  
  gene_info$region=NA
  gene_info$region[match(par_genes,gene_info$gene_id)]="PAR"
  gene_info$class=NA
  gene_info$expressed_allele=NA
  
  library("readxl")
  xci_status=read_excel(file.path(hipsci_datadir,"combined_status.xlsx"),col_names = TRUE)
  
  imprint_status=read_excel(file.path(hipsci_datadir,"geneimprint_22092021.xlsx"),col_names = TRUE)
  gene_names=c()
  status=c()
  expressed_allele=c()
  for (i in 1:dim(imprint_status)[1]) {
    gene_name=unique(c(unlist(strsplit(imprint_status$Gene[i],",",fixed=TRUE)),unlist(strsplit(imprint_status$Aliases[i],",",fixed=TRUE))))
    gene_names=append(gene_names,gene_name)
    status=append(status,rep(imprint_status$Status[i],length(gene_name)))
    expressed_allele=append(expressed_allele,rep(imprint_status$`Expressed Allele`[i],length(gene_name)))
  }
  imprint_status=data.frame(Gene=gene_names,Status=status,Expressed_allele=expressed_allele)
  
  mm=match(imprint_status$Gene,gene_info$gene_name)
  gene_info$class[mm[!is.na(mm)]]=imprint_status$Status[!is.na(mm)]
  gene_info$expressed_allele[mm[!is.na(mm)]]=imprint_status$Expressed_allele[!is.na(mm)]
  rr=which(!is.na(gene_info$expressed_allele[which(gene_info$chr=="X")]))
  ## Remove XIST and TSIX imprint status so that they don't mix with escape status; imprint status is then taken only for autosomal genes
  gene_info[which(gene_info$chr=="X")[rr],]$class=NA
  gene_info[which(gene_info$chr=="X")[rr],]$expressed_allele=NA

  mm3=match(extract_text_before_first_dot(xci_status$`Gene ID`) ,extract_text_before_first_dot(gene_info$gene_id))
  mm3[which(is.na(mm3)==TRUE)]=match(xci_status$`Gene name`[which(is.na(mm3)==TRUE)],gene_info$gene_name)
  gene_info$class[mm3[!is.na(mm3)]]=xci_status$`Reported XCI status`[!is.na(mm3)]
  gene_info$region[mm3[!is.na(mm3)]]=xci_status$Region[!is.na(mm3)]
  
  gene_info$region[which(gene_info$chr!="X" & gene_info$chr!="Y" & gene_info$chr!="MT")]="AUT"
  gene_info$region[which(gene_info$chr=="Y")]="Y"
  gene_info$region[which(gene_info$chr=="MT")]="MT"
  
  # Filter genes to only lincRNA and protein_coding:
  sel_ind=which(((gene_info$gene_biotype=="protein_coding") | (gene_info$gene_biotype =="lincRNA")))
  gene_info=gene_info[sel_ind,]
  counts=counts[sel_ind,]
  
  save(hipsci_info,counts,gene_info,file=file.path(hipsci_datadir,"hipsci_rnaseq_rawData.RData"))
  
}

