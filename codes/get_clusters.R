# -----------------------------------------------------------
# Script Name: get_clusters.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

get_clusters <- function() {
  load("data/lines_ase_summary_data.RData")
  df=data.frame(lines=ase_summary_by_lines$cell_line,xist_group=ase_summary_by_lines$XIST_group,
                xist_expr=ase_summary_by_lines$XIST_logCPM,num_escape_genes=ase_summary_by_lines$num_escape_genes,
                num_observed_genes=ase_summary_by_lines$num_observed_genes)
  # k-means clustering:
  T=as.matrix(cbind(df$xist_expr,df$num_escape_genes/df$num_observed_genes),ncol=2)
  colnames(T) <- c( "xist_expr","esc_ratio")
  T=as.matrix(scale(T))
  df55=as.data.frame(T)
  ##
  # Optimal number of clusters:
  #library(factoextra)
  fviz_nbclust(df55, kmeans, method = "gap_stat", nboot=1000 ,print.summary = TRUE) # Gap-statistics
  ##
  numc=3
  df55$cluster <- as.factor(kmeans(df55, centers=numc,iter.max=10000,nstart=100)$clust)

  df$cluster=df55$cluster
  df$clusterName=as.character(df$cluster)
  df$clusterName[df$cluster==1]="G3"
  df$clusterName[df$cluster==2]="G2"
  df$clusterName[df$cluster==3]="G1"

  write.table(df,file="data/female_clusters.txt",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
}
