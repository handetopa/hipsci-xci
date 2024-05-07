# -----------------------------------------------------------
# Script Name: plot_heatmap.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

plot_heatmap <- function(res, min_nonna_num=10) {
  
  matr=res$sig
  ttt=which((rowSums(!is.na(matr[,which(res$D$xist_group=="High-XIST female")])))>=min_nonna_num & (rowSums(!is.na(matr[,which(res$D$xist_group=="Low-XIST female")])))>=min_nonna_num)
  xist_ind=which(res$gene_info$gene_name=="XIST")
  if (!(xist_ind %in% ttt)) {
    ttt=append(ttt,xist_ind,after=(which(ttt>xist_ind)[1]-1))
  }
  print(length(ttt)) 
  
  res_10=res
  res_10$gene_info=res$gene_info[ttt,]
  res_10$ase_ratios=res$ase_ratios[ttt,]
  res_10$sig=res$sig[ttt,]
  res_10$pval=res$pval[ttt,]
  res_10$gene_minorCount=res$gene_minorCount[ttt,]
  res_10$gene_allCount=res$gene_allCount[ttt,]
  res_10$gene_num_variants=res$gene_num_variants[ttt,]

  matr_filtered=as.matrix(res_10$ase_ratios)
  n_genes=dim(matr_filtered)[1]
  n_samples=dim(matr_filtered)[2]
  gene_ids_filtered=rownames(res_10$ase_ratios)
  
  df.ase=as.data.frame(res_10$ase_ratios[n_genes:1,n_samples:1])
  df.ase$gene_name=res_10$gene_info$gene_name[match(rev(gene_ids_filtered),res_10$gene_info$gene_id)]
  df.ase$gene_name=factor(df.ase$gene_name, levels = df.ase$gene_name)
  df.ase.1=melt(df.ase)
  df.ase.1$cluster_group=rep(res_10$D$xist_group[n_samples:1],each=n_genes)
  df.ase.1$cluster_group=factor(df.ase.1$cluster_group, levels=c("High-XIST female","Low-XIST female"),labels=c("High-XIST females","Low-XIST females"))
  df.ase.1$XIST=rep(res_10$D$xist_cpm_log[n_samples:1],each=n_genes)
  df=df.ase.1
  
  library(ggside)
  library(paletteer)
  colramp=(as.vector(paletteer_c("ggthemes::Red-Blue-White Diverging",  (25-1)))) 
  
  labs <- levels(df$gene_name)
  
  p_ase=ggplot(df, aes(y = as.numeric(gene_name), x = variable)) +
    geom_tile(aes(fill = value)) +
    scale_fill_viridis(name ="", discrete=FALSE, breaks=c(0,0.25,0.5)) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("") +
    ylab("") +
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25,face="bold")) +
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=35), #change legend title font size
          legend.text = element_text(size=30)) +
    guides(fill=guide_colourbar(title="ASE ratio")) +
    scale_x_discrete(expand=c(0,0)) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
  
  p_ase=p_ase + scale_y_continuous(expand=c(0,0),breaks = as.numeric(df$gene_name[seq(1,length(labs),2)]), labels = labs[seq(1,length(labs),2)], sec.axis = sec_axis(~., breaks = as.numeric(df$gene_name[seq(2,length(labs),2)]), labels = labs[seq(2,length(labs),2)]))
  
  p_ase2 = p_ase + geom_xsidetile(data=df[match(unique(df$variable),df$variable),], aes(y=as.numeric(1), xfill = XIST, height=0.1)) +  
    scale_xsidey_continuous(breaks = c(1), labels = "XIST",expand = expansion(add = c(0, 0.3))) +
    scale_xfill_gradientn(colours = colramp,name=expression("XIST / log"[2]*"(CPM)"),trans = 'reverse') + 
    theme(legend.position="bottom", legend.direction = "horizontal",legend.title=element_text(vjust=1)) +
    theme(panel.background = element_rect(fill = "white", colour = "white")) +
    theme(axis.text.x = element_blank()) 
  
  df$cluster_group = factor(df$cluster_group, levels=c("High-XIST females","Low-XIST females"))
  myp= p_ase2 + facet_grid(.~cluster_group, 
                           scales = "free_x",space="free_x", # Let the x axis vary across facets
                           switch = "x") +
    theme(strip.text.y = element_text(size=28,angle = 90)) +
    theme(strip.text.x = element_text(size=35)) +
    theme(plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 
  
  library("readxl")
  escape_genes=as.data.frame(read_excel("data/escape_genes.xlsx",col_names = TRUE, skip = 1))
  human_xci=escape_genes$`Combined XCI status`[match(res_10$gene_info$gene_id,escape_genes$`Gene ID`)]
  human_xci[which(is.na(human_xci))]="Unknown"
  human_xci[which(human_xci=="escape")]="Escape"
  human_xci[which(human_xci=="inactive")]="Inactive"
  human_xci[which(human_xci=="variable")]="Variable"
  xci=as.matrix(human_xci,ncol=1)
  table(xci)
  
  df.xci=data.frame(xci_status=xci[n_genes:1,])
  df.xci$gene_name=res_10$gene_info$gene_name[match(rev(gene_ids_filtered),res_10$gene_info$gene_id)]
  df.xci$gene_name=factor(df.xci$gene_name, levels = df.xci$gene_name)
  df.xci$xci_status=factor(df.xci$xci_status, levels = unique(df.xci$xci_status))
  df.xci$cell="human tissues"
  df.xci.1=melt(df.xci)
  cols=c('red',"darkgoldenrod1",'blue',"black")
  
  p_xci=ggplot(df.xci.1, aes(y = as.numeric(gene_name), x = cell)) +
    geom_tile(aes(fill = xci_status)) +
    scale_fill_manual(name ="Reported\nXCI status", breaks=c("Escape","Variable","Inactive","Unknown"),values=cols) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("") +
    ylab("") +
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25,face="bold")) +
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=35), #change legend title font size
          legend.text = element_text(size=30)) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    theme(legend.position="left")
  
  library("cowplot")
  ppp=ggdraw() +
    draw_plot(myp, x = 0.15, y = 0, width = .85, height = 1) +
    draw_plot(p_xci, x = 0, y = 0.0331, width = .15, height = 0.92) 
  
  ggsave("figures/x_chr_heatmap_ase_ratios.pdf", ppp, width=65,height=70,units="cm",limitsize = FALSE)
  
  return(res_10)
}
