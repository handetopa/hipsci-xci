# -----------------------------------------------------------
# Script Name: get_escape_cell_percentage_plots.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

get_escape_cell_percentage_plots <- function(ddf,human_tissue_inactive_ratio=0.69,human_tissue_variable_ratio=0.16) {
  
  res=ddf
  high_ind=which(res$num_observed_high>0)
  low_ind=which(res$num_observed_low>0)
  highescape_ind=which(res$num_observed_highescape>0)
  
  res$num_escape_all=res$num_escape_high+res$num_escape_low+res$num_escape_highescape
  res$num_observed_all=res$num_observed_high+res$num_observed_low+res$num_observed_highescape
  all_ind=which(res$num_observed_all>0)
  
  esc_high=res$num_escape_high[high_ind]/res$num_observed_high[high_ind]
  esc_low=res$num_escape_low[low_ind]/res$num_observed_low[low_ind]
  esc_highescape=res$num_escape_highescape[highescape_ind]/res$num_observed_highescape[highescape_ind]
  esc_all=res$num_escape_all[all_ind]/res$num_observed_all[all_ind]
  
  df_esc=data.frame(escape_ratio=c(esc_high,esc_highescape,esc_low),xist_group=c(rep("high",length(esc_high)),rep("highescape",length(esc_highescape)),rep("low",length(esc_low))))
  
  p_density=ggplot(df_esc, aes(x = escape_ratio, fill = xist_group)) +
    geom_density(alpha = 0.7, color=NA) +
    theme_minimal() +
    scale_fill_manual(name="Female\ncell lines", values=c("#CC0033", "mediumorchid","#FF99CC"),labels=c("Group 1","Group 2", "Group 3"), breaks=c("high", "highescape", "low")) +
    xlab("Fraction of cell lines where gene ASE > 0.1") +
    ylab("Density") +
    theme(text = element_text(size=20)) 
  
  p_cumulative=ggplot(df_esc, aes(x=escape_ratio, colour = xist_group)) +  
    stat_ecdf(geom="step",linewidth=1.2) +
    theme_minimal() +
    scale_color_manual(name="Female\ncell lines", values=c("#CC0033", "mediumorchid","#FF99CC"),labels=c("Group 1","Group 2", "Group 3"), breaks=c("high", "highescape", "low")) +
    xlab("Fraction of cell lines where gene ASE > 0.1") +
    ylab("Cumulative percentage of genes") +
    #geom_hline(yintercept = human_tissue_inactive_ratio, linetype="dashed", color="gray8",size=1) + # Reference to what is seen as inactive ratio in human tissues
    #geom_hline(yintercept = human_tissue_inactive_ratio+human_tissue_variable_ratio, linetype="dashed", color="gray8",size=1) + # Reference to what is seen as variable ratio in human tissues
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c(0,25,50,75,100)) +
    theme(text = element_text(size=20)) 
  
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/escape_ratio_density_clusters.pdf", p_density, width=18,height=10,units="cm",limitsize = FALSE)
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/escape_ratio_cumulative_with_ref_human_tissues_clusters.pdf", p_cumulative, width=18,height=15,units="cm",limitsize = FALSE)
  
}