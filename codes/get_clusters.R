# -----------------------------------------------------------
# Script Name: get_clusters.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

get_clusters <- function() {
  load("data/dff1.RData")
  df=data.frame(lines=dff1$cell_line,xist_group=dff1$XIST_group, xist_expr=dff1$XIST_logCPM,num_escape_genes=dff1$num_escape_genes,num_observed_genes=dff1$num_observed_genes)
  # k-means clustering:
  T=as.matrix(cbind(df$xist_expr,df$num_escape_genes/df$num_observed_genes),ncol=2)
  colnames(T) <- c( "xist_expr","esc_ratio")
  T=as.matrix(scale(T))
  df55=as.data.frame(T)
  ##
  # Optimal number of clusters:
  library(factoextra)
  fviz_nbclust(df55, kmeans, method = "gap_stat", nboot=1000 ,print.summary = TRUE) # Gap-statistics
  ##
  numc=3
  df55$cluster <- as.factor(kmeans(df55, centers=numc,iter.max=10000,nstart=100)$clust)

  df$cluster=df55$cluster
  df$clusterName=as.character(df$cluster)
  df$clusterName[df$cluster==1]="low"
  df$clusterName[df$cluster==2]="highescape"
  df$clusterName[df$cluster==3]="high"

  write.table(df,file="data/female_clusters.txt",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
}
