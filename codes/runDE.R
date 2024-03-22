runDE <- function(path,dataForDE,comparison="MHL",celltype="ipsc",treatasind=FALSE,usesva=TRUE,by_cluster=TRUE) {

  library(edgeR)
  library(sva)
  library(BiocParallel)
  library(variancePartition)
  library(ggpubr)
  library(ggplot2)
  
  load(dataForDE)
  D_simple$center_name[which(D_simple$center_name=="WELLCOME TRUST SANGER INSTITUTE")]=gsub(" ", "_", "WELLCOME TRUST SANGER INSTITUTE")
  
  if (by_cluster==TRUE) {
    clusters=read.table(file.path(path,"data/female_clusters.txt"),header=TRUE,sep="\t")
    D_simple$xist_group=clusters$clusterName[match(D_simple$lines,clusters$lines)]
    D_simple$xist_group[which(D_simple$sex=="male")]="male"
  }

  table(D_simple$xist_group)
  D_simple_orig=D_simple
  d0_orig=d0
  cpm_log_orig=cpm_log
  gene_info_orig=gene_info
  
  df=data.frame(pass_no=D_simple$passage_number[which(D_simple$sex=="female")],xist=D_simple$xist_cpm_log[which(D_simple$sex=="female")],
                center_name=D_simple$center_name[which(D_simple$sex=="female")],type=D_simple$type[which(D_simple$sex=="female")],
                group=D_simple$xist_group[which(D_simple$sex=="female")])
  if (by_cluster==TRUE) {
    df$group=factor(factor(df$group),levels=c("low","highescape","high"))
  } else {
    df$group=factor(factor(df$group),levels=c("lowXIST","highXIST"))
  }

  if (length(unique(df$center_name))>1) {
    library(ggpubr)
    p=ggplot(df,aes(x=group,y=pass_no)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(aes(color=center_name),size=2,shape=19) +
      xlab("Female XIST group") +
      ylab("Passage number") +
      scale_color_discrete(name = "Center name", labels = c("SC", "Wellcome Trust\nSanger Institute", "C")) +
      if (by_cluster==TRUE) {
        stat_compare_means(comparisons = list(c("low","highescape"),c("low","high"),c("high","highescape")),method = "t.test")
      } else {
        stat_compare_means(method = "t.test",label.x = 1.5, label.y = 42)
      }
    ggsave(filename=file.path(path,paste("figures/boxplot_",celltype,"_",by_cluster,".pdf",sep="")), plot=p, width = 15, height = 12, units = "cm")  
  } else {
    p=ggplot(df,aes(x=group,y=pass_no)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(aes(color=center_name),size=2,shape=17) +
      xlab("Female XIST group") +
      ylab("Passage number") +
      theme(legend.position = "none") +
      stat_compare_means(method = "t.test",label.x = 1.5, label.y = 42)
    ggsave(filename=file.path(path,paste("figures/boxplot_",celltype,".pdf",sep="")), plot=p, width = 10, height = 10, units = "cm")  
  }
  
  D_simple$sex=relevel(as.factor(D_simple$sex), ref="male")
  D_simple$xist_group=relevel(as.factor(D_simple$xist_group), ref="male")
  
  if (comparison=="MF") {
    metadata=data.frame(sex_group=as.factor(D_simple$sex))
    mod <- model.matrix(~0 + sex_group, data=metadata)
  } else {
    metadata=data.frame(sex_group=as.factor(D_simple$xist_group))
    mod <- model.matrix(~0 + sex_group , data=metadata)
  }
  
  rownames(metadata)=colnames(cpm_log)
  param = SnowParam(4, "SOCK", progressbar=TRUE)
  register(param)
  
  mod0 <- model.matrix(~ 1, data=metadata)
  if (usesva==TRUE) {
    svobj <- svaseq(cpm(d0), mod, mod0) 
    sv=as.data.frame(svobj$sv)
    colnames(sv)=paste("sv",as.character(c(1:dim(svobj$sv)[2])),sep="")
    sv=cbind(metadata,sv)
  } else {
    sv=metadata
  }
  sv$donor=as.factor(D_simple$donor)
  metadata=sv
  rownames(metadata)=colnames(cpm_log)
  
  form1=paste(c(0,head(colnames(metadata),-1)),collapse="+")
  if (celltype=="ipsc" & treatasind==FALSE) {
    form1=paste("~",form1,"+ (1|donor)",sep="")
  } else {
    form1=paste("~",form1,sep="")
  }
  form=parse(text = form1)[[1]]
  
  # estimate weights using linear mixed model of dream
  vobjDream = voomWithDreamWeights(d0, form, metadata, BPPARAM=param ) # same as limma:voom() except that it allows random efffects in the formula.
  
  if (comparison=="ML") {
    L = makeContrastsDream(form, metadata, contrasts = c(ML="sex_groupmale - sex_grouplowXIST"))
  } else if (comparison=="MH") {
    L = makeContrastsDream(form, metadata, contrasts = c(MH="sex_groupmale - sex_grouphighXIST"))
  } else if (comparison=="HL") {
    L = makeContrastsDream(form, metadata, contrasts = c(HL="sex_grouphighXIST - sex_grouplowXIST"))
  } else if (comparison=="MF") {
    L = makeContrastsDream(form, metadata, contrasts = c(MF="sex_groupmale - sex_groupfemale"))
  } else if (comparison=="MHL") {
    L = makeContrastsDream(form, metadata, contrasts = c(MHL="sex_groupmale - ((sex_grouphighXIST+sex_grouplowXIST)/2)"))
  } else if (comparison=="MLMH") {
    L = makeContrastsDream(form, metadata, contrasts = c(MLMH="(sex_groupmale-sex_grouplowXIST) - (sex_groupmale-sex_grouphighXIST)"))
  } else if (comparison=="ml") {
    L = makeContrastsDream(form, metadata, contrasts = c(ml="sex_groupmale - sex_grouplow"))
  } else if (comparison=="mh") {
    L = makeContrastsDream(form, metadata, contrasts = c(mh="sex_groupmale - sex_grouphigh"))
  } else if (comparison=="mhe") {
    L = makeContrastsDream(form, metadata, contrasts = c(mhe="sex_groupmale - sex_grouphighescape"))
  } else if (comparison=="mlhe") {
    L = makeContrastsDream(form, metadata, contrasts = c(mlhe="sex_groupmale - (sex_grouplow + sex_grouphighescape)/2"))
  } else if (comparison=="mhhe") {
    L = makeContrastsDream(form, metadata, contrasts = c(mhhe="sex_groupmale - (sex_grouphigh + sex_grouphighescape)/2"))
  } else if (comparison=="mlheh") {
    L = makeContrastsDream(form, metadata, contrasts = c(mlheh="sex_groupmale - (sex_grouplow + sex_grouphighescape + sex_grouphigh)/3"))
  } else if (comparison=="mlh") {
    L = makeContrastsDream(form, metadata, contrasts = c(mlh="sex_groupmale - (sex_grouplow + sex_grouphigh)/2"))
  } else if (comparison=="hhe") {
    L = makeContrastsDream(form, metadata, contrasts = c(hhe="sex_grouphigh - sex_grouphighescape"))
  } else if (comparison=="hl") {
    L = makeContrastsDream(form, metadata, contrasts = c(hl="sex_grouphigh - sex_grouplow"))
  } else if (comparison=="lhe") {
    L = makeContrastsDream(form, metadata, contrasts = c(lhe="sex_grouplow - sex_grouphighescape"))
  } else if (comparison=="lheh") {
    L = makeContrastsDream(form, metadata, contrasts = c(lheh="(sex_grouplow + sex_grouphighescape)/2 - sex_grouphigh"))
  } else if (comparison=="lhhe") {
    L = makeContrastsDream(form, metadata, contrasts = c(lhhe="sex_grouplow - (sex_grouphigh + sex_grouphighescape)/2"))
  } 

  # Visualize contrast matrix
  plotContrasts(L) 
  
  # fit dream model with contrasts
  fit = dream( vobjDream, form, metadata, L, BPPARAM = param)
  fit1 = eBayes(fit)
  
  # get names of available coefficients and contrasts for testing
  colnames(fit1)
  
  levels=c(as.character(unique(metadata$sex_group)),"all")
  nbrs=matrix(0,1,length(levels))
  nbr_colnames=c()
  for (lev in 1:length(levels)) {
    nbr_colnames=append(nbr_colnames,paste(eval(parse(text = "levels[lev]")),"Nbr",sep=""))
    if (levels[lev]!="all") {
      nbrs[lev]=sum(metadata$sex_group==levels[lev])
    } else {
      nbrs[lev]=length(metadata$sex_group)
    }
  }
  nbrs=as.data.frame(matrix(rep(nbrs, each = dim(fit1$coefficients)[1]),nrow=dim(fit1$coefficients)[1]))
  colnames(nbrs)=nbr_colnames
  if (celltype=="ipsc") { 
    ind1=c(which(grepl(comparison, colnames(fit1$coefficients), fixed = TRUE)),which(grepl("sex_group", colnames(fit1$coefficients), fixed = TRUE)))
    metadata_df=as.data.frame(fit1$coefficients[,ind1]) # Is this same as mean(E*weights) ?
    metadata_df_colnames=c(paste(sub('.*sex_group', '', colnames(fit1$coefficients)[ind1]),"Mean",sep=""),
                           paste(sub('.*sex_group', '', colnames(fit1$coefficients)[ind1]),"Var",sep=""))
    metadata_var=lapply(fit1$cov.coefficients.list, function (x) diag((x[ind1,ind1])))
    metadata_var=matrix(unlist(metadata_var), ncol = length(ind1))
    metadata_df=cbind(metadata_df,metadata_var)
    colnames(metadata_df)=metadata_df_colnames
    metadata_df=cbind(metadata_df,nbrs)
    metadata_df$dof=fit1$df.total[,which(colnames(fit1$df.total)==comparison)]
  } else {
    ind1=c(which(grepl(comparison, colnames(fit1$coefficients), fixed = TRUE)),which(grepl("sex_group", colnames(fit1$coefficients), fixed = TRUE)))
    metadata_df=as.data.frame(fit1$coefficients[,ind1]) # Is this same as mean(E*weights) ?
    metadata_df_colnames=c(paste(sub('.*sex_group', '', colnames(fit1$coefficients)[ind1]),"Mean",sep=""))
    colnames(metadata_df)=metadata_df_colnames
    metadata_df=cbind(metadata_df,nbrs)
    metadata_df$dof=fit1$df.total
  }

  top.table <- topTable(fit1,coef=comparison,sort.by = "none", n = Inf, confint=TRUE)
  top.table$logFC_SE=top.table$logFC/top.table$t
  top.table$gene_name=gene_info$gene_name[match(rownames(top.table),gene_info$gene_id)]
  top.table$chrs=gene_info$chr[match(rownames(top.table),gene_info$gene_id)]
  top.table$gene_biotype=gene_info$gene_biotype[match(rownames(top.table),gene_info$gene_id)]
  top.table$region=gene_info$region[match(rownames(top.table),gene_info$gene_id)]
  top.table$class=gene_info$class[match(rownames(top.table),gene_info$gene_id)]
  top.table=cbind(top.table,metadata_df)
  top.table$allMean=rowMeans(vobjDream$E)
  top.table$allVar=rowVars(vobjDream$E)
  
  library(ashr)
  top.table$lfsr=ash(top.table$logFC,top.table$logFC_SE,df=median(top.table$dof))$result$lfsr
  if (by_cluster==FALSE) {
    write.table(top.table,file=file.path(path,paste("results/top.table_final_sva_dream_",celltype,"_",comparison,"_new_",treatasind,"_",usesva,".txt",sep="")),quote=FALSE,sep="\t")
  } else {
    write.table(top.table,file=file.path(path,paste("results/top.table_final_sva_dream_",celltype,"_",comparison,"2_new_",treatasind,"_",usesva,".txt",sep="")),quote=FALSE,sep="\t")
  }
  
}


