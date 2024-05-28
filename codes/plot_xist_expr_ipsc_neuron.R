# -----------------------------------------------------------
# Script Name: plot_xist_expr_ipsc_neuron.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

plot_xist_expr_ipsc_neuron <- function(path) {
  
  hipsci_datadir = file.path(path,"data")
  load(file.path(hipsci_datadir,"data_for_DE_ipsc.RData"))
  clusters=read.table(file.path(hipsci_datadir,"female_clusters.txt"),header=TRUE,sep="\t")
  D_simple$xist_group=clusters$clusterName[match(D_simple$lines,clusters$lines)]
  D_simple$xist_group[which(D_simple$sex=="male")]="male"
  D_ipsc=D_simple
  
  load(file.path(hipsci_datadir,"data_for_DE_difneuron.RData"))
  D_neuron=D_simple
  
  df=data.frame(ipsc_xist=D_ipsc$xist_cpm_log[match(D_neuron$lines,D_ipsc$lines)],
              neuron_xist=D_neuron$xist_cpm_log,
              ipsc_groups=D_ipsc$xist_group[match(D_neuron$lines,D_ipsc$lines)],neuron_groups=D_neuron$xist_group)
  df=df[-which(df$ipsc_groups=="male"),]
  df$ipsc_xist_group=df$ipsc_groups
  df$ipsc_xist_group[which(df$ipsc_xist<1.5)]="Low-XIST"
  df$ipsc_xist_group[which(df$ipsc_xist>=1.5)]="High-XIST"
  df$ipsc_xist_group=factor(df$ipsc_xist_group)
  df$ipsc_groups=factor(df$ipsc_groups)
  df$neuron_groups=factor(df$neuron_groups)
  
  p=ggplot(df,aes(x=ipsc_xist,y=neuron_xist,col=ipsc_xist_group)) +
    geom_point(size=3) +
    geom_hline(yintercept=1.5,linetype="dashed") +
    geom_vline(xintercept=1.5,linetype="dashed") +
    xlab(expression("iPSCs - XIST expression (log"[2]*"(CPM+0.5))")) +
    ylab(expression("Diff. neurons - XIST expression (log"[2]*"(CPM+0.5))")) +
    scale_color_manual(name="iPSC groups", values=c("#CC0033","#FF99CC"),breaks=c("High-XIST","Low-XIST"),labels=c("High XIST","Low XIST")) +
    xlim(-8,8) +
    ylim(-8,8) +
    theme_minimal() +
    theme(text = element_text(size=20))
  
  ggsave("figures/ipsc_neuron_xist_expr.pdf", p, width=20,height=18,units="cm",limitsize = FALSE) 

}
