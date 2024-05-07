# -----------------------------------------------------------
# Script Name: run_all.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

run_all <- function() {
  
  #BiocManager::install('PCAtools')
  library(PCAtools)
  library(ggplot2)
  
  path = "hipsci-xci"
  setwd(path)
  source(file.path(path,"codes/get_rawData.R"))
  source(file.path(path,"codes/get_data_byLineCenter.R"))
  source(file.path(path,"codes/getDataForDE.R"))
  source(file.path(path,"codes/plot_pca.R"))
  source(file.path(path,"codes/plot_marker_genes.R"))
  source(file.path(path,"codes/runDE_all.R"))
  
  get_rawData(path = path)
  get_data_byLineCenter(path = path)
  getDataForDE(path = path)
  
  plot_pca(chr = "X", path = path)
  plot_pca(chr = "aut", path = path)
  plot_marker_genes(path = path)
  
  runDE_all(path = path)
  
  source(file.path(path,"codes/get_snpwise_ase.R"))
  ase_files=list.files(path="data/ase_files/",pattern="*.table",full.names=TRUE)
  for (i in seq(1,length(ase_files))) {
    print(i)
    ase_fileName1=ase_files[i]
    get_snpwise_ase(ase_fileName1)
  }
  source(file.path(path,"codes/get_ase_matrix.R"))
  res=get_ase_matrix(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=0)
  source(file.path(path,"codes/plot_heatmap.R"))
  res_10=plot_heatmap(res=res, min_nonna_num=10)
  source(file.path(path,"codes/plot_heatmap_sig.R"))
  res_10=plot_heatmap_sig(res=res, min_nonna_num=10)
  source(file.path(path,"codes/plot_heatmap_expr.R"))
  res_10=plot_heatmap_sig(res=res, min_nonna_num=10)
  plot_heatmap_expr(res_10=res_10)
  
  source(file.path(path,"codes/get_ase_matrix.R"))
  res=get_ase_matrix(include="all",mychr="X",mysex="female",min_ase_ratio_for_escape=0.1,xist_lim=1.5,include.na=FALSE,genes_orderby="pos",samples_orderby="xist",altern_hypt="greater",min_nonna_num=0)
  ttX=as.data.frame(read_xlsx("data/escape_genes.xlsx",skip=1))
  human_xci=ttX$`Combined XCI status`
  names(human_xci)=esc$`Gene ID`
  human_xci[human_xci=="escape"]="Escape"
  human_xci[human_xci=="variable"]="Variable"
  human_xci[human_xci=="inactive"]="Inactive"
  table(human_xci)
  clusters=read.table("data/female_clusters.txt",header=TRUE,sep="\t")
  clusts=clusters$clusterName[match(res$D$lines,clusters$lines)]
  ddf=data.frame(gene_id=res$gene_info$gene_id,
                 num_escape_low=rowSums(res$sig[,clusts=="low"],na.rm=TRUE),
                 num_observed_low=rowSums(!is.na(res$sig[,clusts=="low"])),
                 num_lines_low=sum(clusts=="low"),
                 num_escape_highescape=rowSums(res$sig[,clusts=="highescape"],na.rm=TRUE),
                 num_observed_highescape=rowSums(!is.na(res$sig[,clusts=="highescape"])),
                 num_lines_highescape=sum(clusts=="highescape"),
                 num_escape_high=rowSums(res$sig[,clusts=="high"],na.rm=TRUE),
                 num_observed_high=rowSums(!is.na(res$sig[,clusts=="high"])),
                 num_lines_high=sum(clusts=="high"))
  
  source(file.path(path,"codes/get_escape_cell_percentage_plots.R"))
  source(file.path(path,"codes/human_tissues_vs_hipsci_comparison.R"))
  get_escape_cell_percentage_plots(ddf)
  human_tissues_vs_hipsci_comparison(human_xci,ddf,th_inactive=0,th_esc=0.8)
  
  source(file.path(path,"codes/violin_abslogfc_femalebiased_in_G1G2G3.R"))
  source(file.path(path,"codes/violin_abslogfc_malebiased_in_G1.R"))
  violin_abslogfc_femalebiased_in_G1G2G3()
  violin_abslogfc_malebiased_in_G1()
  
  source(file.path(path,"codes/boxplot_esc_ratio_per_xci_group.R"))
  boxplot_esc_ratio_per_xci_group(human_xci)
  source(file.path(path,"codes/boxplot_esc_ratio_in_female_groups_per_xci_group.R"))
  boxplot_esc_ratio_in_female_groups_per_xci_group(human_xci)
  
  source(file.path(path,"codes/std_gene_expr_by_groups.R"))
  std_gene_expr_by_groups()
  
  source(file.path(path,"codes/mean_ase_by_groups.R"))
  mean_ase_by_groups()
  source(file.path(path,"codes/esc_ratio_by_xist_groups.R"))
  esc_ratio_by_xist_groups()
  
  source(file.path(path,"codes/plot_xist_ase.R"))
  plot_xist_ase()
  source(file.path(path,"codes/density_aut_chrX.R"))
  density_aut_chrX()
}