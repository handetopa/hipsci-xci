# -----------------------------------------------------------
# Script Name: plot_xist_ase.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

plot_xist_ase <- function() {
  
  source(file.path(path,"codes/get_ase_matrix.R"))
  res=get_ase_matrix(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=0)
  source(file.path(path,"codes/get_ase_matrix.R"))
  res_xist=get_ase_matrix(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.05,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="less",min_nonna_num=0,gene="XIST")
  
  med_ase_x=apply(res$ase_ratios,2,median,na.rm=TRUE)
  df_xist=data.frame(xist=res_xist$ase_ratios,l=pmax(0,res_xist$low_l),u=pmin(1,res_xist$up_l),med_ase_x=med_ase_x,xist_group=D$xist_group)
  df_xist$l[is.na(df_xist$xist)]=(-0.1)
  df_xist$u[is.na(df_xist$xist)]=(-0.1)
  df_xist$xist[is.na(df_xist$xist)]=(-0.1)
  
  p_xist_ase=ggplot(df_xist,aes(x=xist,y=med_ase_x,xmin=l,xmax=u,col=xist_group)) +
    geom_point(alpha=0.8,size=1.6) +
    geom_errorbarh(aes(xmin=l, xmax=u)) + 
    ylab("Median ASE in chrX\ngenes per line") + 
    xlab("XIST ASE (minor/total)") +
    scale_x_continuous(breaks=c(-0.1,0,0.1,0.2,0.3,0.4,0.5),labels=c("NA","0.0","0.1","0.2","0.3","0.4","0.5")) +
    scale_color_manual(name="Female\ncell lines", values=c("#CC0033", "#FF99CC"),labels=c("High-XIST", "Low-XIST"), breaks=c("High-XIST female", "Low-XIST female")) +
    ylim(0,0.5) +
    theme_minimal() +
    theme(text = element_text(size=20))
  
  p_x_ase_density=ggplot(df_xist, aes(y = med_ase_x, fill = xist_group)) +
    geom_density(alpha = 0.8, color=NA) +
    scale_fill_manual(name="Female\ncell lines", values=c("#CC0033", "#FF99CC"),labels=c("High-XIST", "Low-XIST"), breaks=c("High-XIST female", "Low-XIST female")) +
    ylab("") +
    xlab("Density") +
    ylim(0,0.5) +
    theme_minimal() +
    theme(text = element_text(size=20))
  
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/xist_ase.pdf", p_xist_ase, width=15,height=10,units="cm",limitsize = FALSE)
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/x_ase_density.pdf", p_x_ase_density, width=10,height=10,units="cm",limitsize = FALSE)
  
  ase_ratios_clonal=ase_ratios[,which(res_xist$pval<0.05)]
  ase_ratios_nonclonal=ase_ratios[,which(res_xist$pval>=0.05)]
  sig_clonal=sig[,which(res_xist$pval<0.05)]
  sig_nonclonal=sig[,which(res_xist$pval>=0.05)]
  dim(ase_ratios_clonal)
  dim(ase_ratios_nonclonal)
  
  esc_ratio_clonal=colSums(sig_clonal,na.rm=TRUE)/colSums(!is.na(sig_clonal))
  esc_ratio_nonclonal=colSums(sig_nonclonal,na.rm=TRUE)/colSums(!is.na(sig_nonclonal))
  
  num_genes_clonal=colSums(!is.na(sig_clonal))
  num_genes_nonclonal=colSums(!is.na(sig_nonclonal))
  
  med_ase_clonal=apply(ase_ratios_clonal,2,median,na.rm=TRUE)
  med_ase_nonclonal=apply(ase_ratios_nonclonal,2,median,na.rm=TRUE)
  
  df=data.frame(num_genes=c(num_genes_clonal,num_genes_nonclonal),median_ase=c(med_ase_clonal,med_ase_nonclonal),
              esc_ratio=c(esc_ratio_clonal,esc_ratio_nonclonal),group=c(rep("clonal",length(med_ase_clonal)),rep("nonclonal",length(med_ase_nonclonal))))
  
  p1=ggplot(df,aes(x=group,y=esc_ratio)) +
    geom_boxplot(outlier.shape = NA) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1,binwidth=0.02,fill="#CC0033") +
    theme_minimal() +
    ylim(0,1) +
    ylab("Fraction of genes with\nASE > 0.1") +
    scale_x_discrete(name="", breaks=c("clonal","nonclonal"), labels=c("XIST ASE < 0.05", expression("XIST ASE ">=" 0.05"))) +
    stat_compare_means(size=6) +
    theme(text = element_text(size=20)) 
  
  p2=ggplot(df,aes(x=group,y=median_ase)) +
    geom_boxplot(outlier.shape = NA) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1,binwidth=0.01,fill="#CC0033") +
    theme_minimal() +
    ylim(0,0.5) +
    ylab("Median ASE of chrX\nper line") +
    scale_x_discrete(name="", breaks=c("clonal","nonclonal"), labels=c("XIST ASE < 0.05", expression("XIST ASE ">=" 0.05"))) +
    stat_compare_means(size=6) +
    theme(text = element_text(size=20)) 
  
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/clonality_xchr_escape_ratio_unicolor.pdf", p1, width=20,height=15,units="cm",limitsize = FALSE)
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/clonality_xchr_median_ase_unicolor.pdf", p2, width=20,height=15,units="cm",limitsize = FALSE)

}