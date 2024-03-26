density_aut_chrX <- function() {
  
  source(file.path(path,"codes/get_ase_matrix.R"))
  chr=1
  res_aut=get_ase_matrix(include="all",mychr=as.character(chr),mysex="all",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",altern_hypt="less",min_nonna_num=0)
  DD=res_aut$D
  ase_aut=res_aut$ase_ratios
  for (chr in 2:22) {
    print(chr)
    res_aut=get_ase_matrix(include="all",mychr=as.character(chr),mysex="all",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",altern_hypt="less",min_nonna_num=0)
    ase_aut=rbind(ase_aut,res_aut$ase_ratios)
  }
  dim(ase_aut)[1] # Number of autosomal genes with ASE info: 12500
  sum(rowSums(!is.na(ase_aut[,DD$xist_group=="Male"]))>0)  # Number of autosomal genes with ASE info in males: 11744
  sum(rowSums(!is.na(ase_aut[,DD$xist_group=="Low-XIST female"]))>0)  # Number of autosomal genes with ASE info in low-XIST: 11039
  sum(rowSums(!is.na(ase_aut[,DD$xist_group=="High-XIST female"]))>0)  # Number of autosomal genes with ASE info in high-XIST: 11807
  ##
  common_ind=intersect(which(rowSums(!is.na(ase_aut[,DD$xist_group=="Low-XIST female"]))>0),
                       intersect(which(rowSums(!is.na(ase_aut[,DD$xist_group=="Male"]))>0),
                                 which(rowSums(!is.na(ase_aut[,DD$xist_group=="High-XIST female"]))>0)))
  length(common_ind) # Number of autosomal genes common in all: 10516
  DF_common=data.frame(ase_ratios=c(c(ase_aut[common_ind,which(DD$xist_group=="Male")]),c(ase_aut[common_ind,which(DD$xist_group=="High-XIST female")]),
                                    c(ase_aut[common_ind,which(DD$xist_group=="Low-XIST female")])),
                       xist_group=c(rep("Male",length(common_ind)*sum(DD$xist_group=="Male")),rep("High-XIST female",length(common_ind)*sum(DD$xist_group=="High-XIST female")),rep("Low-XIST female",length(common_ind)*sum(DD$xist_group=="Low-XIST female"))))
  ##
  NN_male=dim(ase_aut)[1]*sum(DD$xist_group=="Male")
  NN_low=dim(ase_aut)[1]*sum(DD$xist_group=="Low-XIST female")
  NN_high=dim(ase_aut)[1]*sum(DD$xist_group=="High-XIST female")
  
  DF=data.frame(ase_ratios=c(c(ase_aut[,which(DD$xist_group=="Male")]),c(ase_aut[,which(DD$xist_group=="High-XIST female")]),
                           c(ase_aut[,which(DD$xist_group=="Low-XIST female")])),
              xist_group=c(rep("Male",NN_male),rep("High-XIST female",NN_high),rep("Low-XIST female",NN_low)))
  
  ks.test(c(ase_aut[,which(DD$xist_group=="Male")]),c(ase_aut[,which(DD$xist_group=="High-XIST female")])) #Kolmogorov-Smirnov p.value for equality of distributions: 0.2196
  ks.test(c(ase_aut[,which(DD$xist_group=="Male")]),c(ase_aut[,which(DD$xist_group=="Low-XIST female")])) #Kolmogorov-Smirnov p.value for equality of distributions: 0.747
  ks.test(c(ase_aut[,which(DD$xist_group=="Low-XIST female")]),c(ase_aut[,which(DD$xist_group=="High-XIST female")])) #Kolmogorov-Smirnov p.value for equality of distributions: 0.1964
  
  p_A=ggplot(DF, aes(x=ase_ratios, fill=xist_group)) +
    geom_density(alpha=0.5) +
    scale_fill_manual(name="Female\ncell lines", values=c("#CC0033", "#FF99CC", "#3399FF"),labels=c("High-XIST females", "Low-XIST females", "Males"), breaks=c("High-XIST female", "Low-XIST female", "Male")) +
    xlab("ASE per gene per cell line\n(Autosomal chromosomes)") +
    ylab("Density") +
    theme_minimal() +
    theme(text = element_text(size=20))
  
  res=get_ase_matrix(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=0)
  ase_ratios_X=res$ase_ratios
  DD=res$D
  NN_low=dim(ase_ratios_X)[1]*sum(DD$xist_group=="Low-XIST female")
  NN_high=dim(ase_ratios_X)[1]*sum(DD$xist_group=="High-XIST female")
  DF_X=data.frame(ase_ratios=c(c(ase_ratios_X[,which(DD$xist_group=="High-XIST female")]),
                               c(ase_ratios_X[,which(DD$xist_group=="Low-XIST female")])),
                  xist_group=c(rep("High-XIST female",NN_high),rep("Low-XIST female",NN_low)))
  
  dim(ase_ratios_X)[1] # Number of chrX genes with ASE info: 411
  sum(rowSums(!is.na(ase_ratios_X[,DD$xist_group=="Low-XIST female"]))>0)  # Number of chrX genes with ASE info in low-XIST: 344
  sum(rowSums(!is.na(ase_ratios_X[,DD$xist_group=="High-XIST female"]))>0)  # Number of chrX genes with ASE info in high-XIST: 390
  ks.test(c(ase_ratios_X[,which(DD$xist_group=="High-XIST female")]),
          c(ase_ratios_X[,which(DD$xist_group=="Low-XIST female")])) #Kolmogorov-Smirnov p.value for equality of distributions: <2.2e-16
  ##
  common_ind=intersect(which(rowSums(!is.na(ase_ratios_X[,DD$xist_group=="Low-XIST female"]))>0),
                                 which(rowSums(!is.na(ase_ratios_X[,DD$xist_group=="High-XIST female"]))>0))
  length(common_ind) # Number of autosomal genes common in all: 323
  DF_X_common=data.frame(ase_ratios=c(c(ase_ratios_X[common_ind,which(DD$xist_group=="High-XIST female")]),
                                    c(ase_ratios_X[common_ind,which(DD$xist_group=="Low-XIST female")])),
                       xist_group=c(rep("High-XIST female",length(common_ind)*sum(DD$xist_group=="High-XIST female")),rep("Low-XIST female",length(common_ind)*sum(DD$xist_group=="Low-XIST female"))))
  ##
  
  p_X=ggplot(DF_X, aes(x=ase_ratios, fill=xist_group)) +
    geom_density(alpha=0.8) +
    scale_fill_manual(name="Female\ncell lines", values=c("#CC0033", "#FF99CC"),labels=c("High-XIST", "Low-XIST"), breaks=c("High-XIST female", "Low-XIST female")) +
    xlab("ASE per gene per cell line\n(X chromosome)") +
    ylab("Density") +
    theme_minimal() +
    theme(text = element_text(size=20))
  
  
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/density_autosomes.pdf", p_A, width=16,height=12,units="cm",limitsize = FALSE)
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/density_X.pdf", p_X, width=15,height=12,units="cm",limitsize = FALSE)

}