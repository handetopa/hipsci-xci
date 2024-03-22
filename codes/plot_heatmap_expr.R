plot_heatmap_expr <- function(res_10) {
  
  library(paletteer)
  source("/Users/topah/Desktop/hipsci_codes/codes/extract_text_before_first_dot.R")
  load("/Users/topah/Desktop/hipsci_codes/data/data_for_DE_new.RData")
  gene_ids_filtered2=rownames(res_10$sig)
  cpm_log_filtered=cpm_log[match(extract_text_before_first_dot(gene_ids_filtered2),extract_text_before_first_dot(rownames(cpm_log))),match(res_10$D$lines,sub('_[^_]*$', '', colnames(cpm_log)))]
  rm.ind=which(rowSums(is.na(cpm_log_filtered))>0)
  print(paste("Missing expression data for gene: ",res_10$gene_info$gene_id[match(gene_ids_filtered2[rm.ind],res_10$gene_info$gene_id)],
            " (",res_10$gene_info$gene_name[match(gene_ids_filtered2[rm.ind],res_10$gene_info$gene_id)],")",sep=""))
  cpm_log_filtered=cpm_log_filtered[-rm.ind,]
  n_genes=dim(cpm_log_filtered)[1]
  n_samples=dim(cpm_log_filtered)[2]
  
  cpm_log_filtered=t(scale(t(cpm_log_filtered)))
  df.ase=as.data.frame(cpm_log_filtered[n_genes:1,n_samples:1])
  df.ase$gene_name=res_10$gene_info$gene_name[match(rev(gene_ids_filtered2[-rm.ind]),res_10$gene_info$gene_id)]
  df.ase$gene_name=factor(df.ase$gene_name, levels = df.ase$gene_name)
  df.ase.1=melt(df.ase)
  df.ase.1$cluster_group=rep(res_10$D$xist_group[n_samples:1],each=n_genes)
  df.ase.1$cluster_group=factor(df.ase.1$cluster_group, levels=c("High-XIST female","Low-XIST female"),labels=c("High-XIST females","Low-XIST females"))
  df.ase.1$XIST=rep(res_10$D$xist_cpm_log[n_samples:1],each=n_genes)
  df=df.ase.1
  
  colramp=(as.vector(paletteer_c("ggthemes::Red-Blue-White Diverging",  (25-1))))
  labs <- levels(df$gene_name)
  
  p_ase=ggplot(df, aes(y = as.numeric(gene_name), x = variable)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradientn(colours = rev(colramp),limits=c(-5,5),breaks=c(-5,0,5),labels=c("<-5","0",">5")) +
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
    guides(fill=guide_colourbar(title="Standardized gene\nexpression")) +
    scale_x_discrete(expand=c(0,0)) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
  
  p_ase=p_ase + scale_y_continuous(expand=c(0,0),breaks = as.numeric(df$gene_name[seq(2,length(labs),2)]), labels = labs[seq(2,length(labs),2)], sec.axis = sec_axis(~., breaks = as.numeric(df$gene_name[seq(1,length(labs),2)]), labels = labs[seq(1,length(labs),2)]))
  
  p_ase2 = p_ase + geom_xsidetile(data=df[match(unique(df$variable),df$variable),], aes(y=as.numeric(1), xfill = XIST, height=0.1)) +
    scale_xsidey_continuous(breaks = c(1), labels = "XIST",expand = expansion(add = c(0, 0.3))) +
    scale_xfill_gradientn(colours = colramp,name=expression("XIST / log"[2]*"(CPM)"),trans="reverse") +
    theme(legend.position="bottom", legend.direction = "horizontal",legend.title=element_text(vjust=1)) +
    theme(panel.background = element_rect(fill = "white", colour = "white")) +
    theme(axis.text.x = element_blank()) 
  myp= p_ase2 + facet_grid(.~cluster_group, 
                         scales = "free_x",space="free_x", # Let the x axis vary across facets
                         switch = "x") +
    theme(strip.text.y = element_text(size=28,angle = 90)) +
    theme(strip.text.x = element_text(size=35)) +
    theme(plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/x_chr_heatmap_gene_expression.pdf", myp, width=55,height=70,units="cm",limitsize = FALSE)
  
  
  top.table=read.table("/Users/topah/Desktop/hipsci_codes/results/top.table_final_sva_dream_ipsc_HL_new_FALSE_TRUE.txt")
  adj.p.val=top.table$adj.P.Val[match(extract_text_before_first_dot(gene_ids_filtered2[-rm.ind]),extract_text_before_first_dot(rownames(top.table)))]
  logFC=top.table$logFC[match(extract_text_before_first_dot(gene_ids_filtered2[-rm.ind]),extract_text_before_first_dot(rownames(top.table)))]
  
  logFC_original=top.table$logFC[match(extract_text_before_first_dot(gene_ids_filtered2),extract_text_before_first_dot(rownames(top.table)))]
  adj.p.val_original=top.table$adj.P.Val[match(extract_text_before_first_dot(gene_ids_filtered2),extract_text_before_first_dot(rownames(top.table)))]
  ase_diff=apply(res_10$ase_ratios[,which(res_10$D$xist_group=="Low-XIST female")],1,median,na.rm=TRUE)-apply(res_10$ase_ratios[,which(res_10$D$xist_group=="High-XIST female")],1,median,na.rm=TRUE)
   
  sex_groups=names(table(res_10$D$xist_group))
  wilcox.pvalue=matrix(NA,n_genes+1,1)
  fisher.pvalue=matrix(NA,n_genes+1,1)
  esc_diff=matrix(NA,n_genes+1,1)
  for (i in 1:(n_genes+1)) {
    M <- as.table(rbind(c(sum(res_10$sig[i,which(res_10$D$xist_group==sex_groups[1])],na.rm=TRUE),sum(!res_10$sig[i,which(res_10$D$xist_group==sex_groups[1])],na.rm=TRUE)),
                        c(sum(res_10$sig[i,which(res_10$D$xist_group==sex_groups[2])],na.rm=TRUE),sum(!res_10$sig[i,which(res_10$D$xist_group==sex_groups[2])],na.rm=TRUE))))
    dimnames(M) <- list(xist_group = c(sex_groups[1],sex_groups[2]),
                        status = c("Escape","Inactive"))
    esc_diff[i]=(M[,1]/rowSums(M))[2]-(M[,1]/rowSums(M))[1]
    tryCatch({
      fisher.pvalue[i]=fisher.test(M)$p.value
    }, error=function(e){})
    tryCatch({
      wilcox.pvalue[i]=wilcox.test(res_10$ase_ratios[i,which(res_10$D$xist_group==sex_groups[2])],res_10$ase_ratios[i,which(res_10$D$xist_group==sex_groups[1])],exact=FALSE)$p.value
    }, error=function(e){})
  }
  wilcox.pvalue_orig=wilcox.pvalue
  wilcox.pvalue=wilcox.pvalue[-rm.ind]
  
  fisher.pvalue_orig=fisher.pvalue
  fisher.pvalue=fisher.pvalue[-rm.ind]
  
  df.diff=data.frame(logFC=(-logFC_original),adj.p.val=adj.p.val_original,diff_median_ase=ase_diff,
                     wilcox.pvalue=wilcox.pvalue_orig,esc_diff=esc_diff,fisher.pvalue=fisher.pvalue_orig,
                     genes=res_10$gene_info$gene_name[match(gene_ids_filtered2,res_10$gene_info$gene_id)])
  
  df.diff$ase_sig=(df.diff$wilcox.pvalue<0.05)
  df.diff$esc_sig=(df.diff$fisher.pvalue<0.05)
  df.diff$logFC[which(is.na(df.diff$logFC)==TRUE)]=0
  df.diff$diff_median_ase[which(is.na(df.diff$diff_median_ase)==TRUE)]=0 #XIST
  df.diff$esc_diff[which(is.na(df.diff$esc_diff)==TRUE)]=0 #XIST
  df.diff$expr_sig=(df.diff$adj.p.val<0.05)
  
  df.median_ASE_diff=df.diff
  df.median_ASE_diff$ase_sign=""
  df.median_ASE_diff$ase_sign[df.median_ASE_diff$ase_sig==TRUE]="*"
  df.median_ASE_diff$genes_ase=paste(df.median_ASE_diff$ase_sign, df.median_ASE_diff$genes,sep=" ")
  df.median_ASE_diff$genes_ase <- factor(df.median_ASE_diff$genes_ase, levels = df.median_ASE_diff$genes_ase)
  df.median_ASE_diff$ase_sign[df.median_ASE_diff$ase_sig==TRUE]="#333333"
  df.median_ASE_diff$ase_sign[df.median_ASE_diff$ase_sig==FALSE]="#CCCCCC"
  
  df.median_ASE_diff$expr_sign=""
  df.median_ASE_diff$expr_sign[df.median_ASE_diff$expr_sig==TRUE]="*"
  df.median_ASE_diff$genes_expr=paste(df.median_ASE_diff$expr_sign, df.median_ASE_diff$genes,sep=" ")
  df.median_ASE_diff$genes_expr <- factor(df.median_ASE_diff$genes_expr, levels = df.median_ASE_diff$genes_expr)
  df.median_ASE_diff$expr_sign[df.median_ASE_diff$expr_sig==TRUE]="#333333"
  df.median_ASE_diff$expr_sign[df.median_ASE_diff$expr_sig==FALSE]="#CCCCCC"
  
  df.median_ASE_diff$esc_sign=""
  df.median_ASE_diff$esc_sign[df.median_ASE_diff$esc_sig==TRUE]="*"
  df.median_ASE_diff$genes_esc=paste(df.median_ASE_diff$esc_sign, df.median_ASE_diff$genes,sep=" ")
  df.median_ASE_diff$genes_esc <- factor(df.median_ASE_diff$genes_esc, levels = df.median_ASE_diff$genes_esc)
  df.median_ASE_diff$esc_sign[df.median_ASE_diff$esc_sig==TRUE]="#333333"
  df.median_ASE_diff$esc_sign[df.median_ASE_diff$esc_sig==FALSE]="#CCCCCC"
  
  df.median_ASE_diff=df.median_ASE_diff[dim(df.median_ASE_diff)[1]:1,]
  df.median_ASE_diff$genes_esc <- factor(df.median_ASE_diff$genes_esc, levels = df.median_ASE_diff$genes_esc)
  df.median_ASE_diff$genes_expr <- factor(df.median_ASE_diff$genes_expr, levels = df.median_ASE_diff$genes_expr)
  df.median_ASE_diff$genes_ase <- factor(df.median_ASE_diff$genes_ase, levels = df.median_ASE_diff$genes_ase)
  df.median_ASE_diff$genes <- factor(df.median_ASE_diff$genes, levels = df.median_ASE_diff$genes)
  
  pp_ase_diff=ggplot(data=df.median_ASE_diff, aes(y=genes_ase, x=diff_median_ase,fill=ase_sig,na.rm=FALSE)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("#CCCCCC","#333333"),breaks=c("FALSE","TRUE"),labels=c(expression(p >= 0.05), "p < 0.05")) +
    scale_x_continuous(position="top") +
    ylab("") +
    xlab("Difference in median ASE") +
    guides(fill=guide_legend(title="")) +
    theme_minimal() +
    theme(axis.text.y = element_blank()) +
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25),
          legend.title = element_text(size=25),
          legend.text = element_text(size=25)) +
    theme(panel.background = element_rect(fill = "white", colour = "white")) 
  
  print(paste("XIST logFC (",df.median_ASE_diff$logFC[which(df.median_ASE_diff$logFC<(-0.5))],") cropped to -0.5.",sep=""))
  df.median_ASE_diff$logFC[which(df.median_ASE_diff$logFC<(-0.5))]=(-0.5) # Crop XIST logFC to -0.5
  df.median_ASE_diff$genes_expr <- factor(df.median_ASE_diff$genes_expr, levels = df.median_ASE_diff$genes_expr)
  df.median_ASE_diff$genes_ase <- factor(df.median_ASE_diff$genes_ase, levels = df.median_ASE_diff$genes_ase)
  df.median_ASE_diff$genes_esc <- factor(df.median_ASE_diff$genes_esc, levels = df.median_ASE_diff$genes_esc)
  df.median_ASE_diff$genes <- factor(df.median_ASE_diff$genes, levels = df.median_ASE_diff$genes)
  
  pp_esc_diff=ggplot(data=df.median_ASE_diff, aes(y=genes_esc, x=esc_diff,fill=esc_sig,na.rm=FALSE)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("#CCCCCC","#333333"),breaks=c("FALSE","TRUE"),labels=c(expression(p >= 0.05), "p < 0.05")) +
    scale_x_continuous(position="top") +
    ylab("") +
    xlab("Difference in fraction of\ngenes with ASE ratio > 0.1") +
    guides(fill=guide_legend(title="")) +
    theme_minimal() +
    theme(axis.text.y = element_blank()) +
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25),
          legend.title = element_text(size=25),
          legend.text = element_text(size=25)) +
    theme(panel.background = element_rect(fill = "white", colour = "white")) 
  
  pp_expr_diff=ggplot(data=df.median_ASE_diff[-seq(dim(df.median_ASE_diff)[1],1)[rm.ind],], aes(y=genes_expr, x=logFC,fill=expr_sig,na.rm=FALSE)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("#CCCCCC","#333333"),breaks=c("FALSE","TRUE"),labels=c(expression(adj.p >= 0.05), "adj.p < 0.05")) +
    ylab("") +
    xlab(expression(log[2]*"(Fold change)")) +
    scale_x_continuous(breaks=c(-0.5,-0.25,0,0.25,0.5), labels=c("< -0.5","-0.25","0","0.25","0.5"), limits=c(-0.5,0.5),position="top") +
    guides(fill=guide_legend(title="")) +
    theme_minimal() +
    theme(axis.text.y = element_blank()) +
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25),
          legend.title = element_text(size=25),
          legend.text = element_text(size=25)) +
    theme(panel.background = element_rect(fill = "white", colour = "white")) 
  
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/diff_median_ase.pdf", pp_ase_diff, height=70,width=20,units="cm",limitsize = FALSE)
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/diff_expr.pdf", pp_expr_diff, height=70,width=20,units="cm",limitsize = FALSE)
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/diff_esc_ratio.pdf", pp_esc_diff, height=70,width=20,units="cm",limitsize = FALSE)
  
  print(paste("Number of significant genes according to difference btw ASE: ",sum(df.median_ASE_diff$ase_sig,na.rm=TRUE),sep=""))
  print(paste("Number of significant genes according to difference btw expression: ",sum(df.median_ASE_diff[-seq(dim(df.median_ASE_diff)[1],1)[rm.ind],]$expr_sig),sep=""))
  print(paste("Number of significant genes according to difference btw escape ratio: ",sum(df.median_ASE_diff$esc_sig),sep=""))
  print(paste("Number of significant genes according to difference btw escape ratio and expression: ",sum(df.median_ASE_diff$esc_sig & df.median_ASE_diff$expr_sig,na.rm=TRUE),sep=""))
  table1=as.matrix(table(df.median_ASE_diff$esc_sig, df.median_ASE_diff$expr_sig))
  table1[2,1]=table1[2,1]+1
  table1
  print(paste("Fisher's exact test p-value: ",fisher.test(table1)$p.value,sep=""))
  
  overall_p <- function(my_model) {
    f <- summary(my_model)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  BB=df.median_ASE_diff[-seq(dim(df.median_ASE_diff)[1],1)[rm.ind],]
  BB=BB[-which(BB$genes=="XIST"),]
  model1=lm(BB$logFC~BB$esc_diff)
  
  library(ggpmisc)
  p1=ggplot(BB,aes(x=esc_diff,y=logFC)) +
    stat_poly_line() +
    stat_poly_eq(use_label(c("eq", "adj.R2", "p"))) +
    geom_point() +
    ylim(-0.2,0.6) +
    theme_minimal() +
    xlab("Difference in fraction of lines\nwhere gene ASE > 0.1") +
    ylab(expression("log"[2]*"(Fold changes)")) +
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18), 
          legend.title = element_text(size=18),
          legend.text = element_text(size=18),
          plot.title = element_text(size=18))
  
  print(paste("p-value for regression model: ",overall_p(model1),sep=""))
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/esc_ratio_diff_logFC_chrX.pdf", p1, width=15,height=15,units="cm",limitsize = FALSE)
  
}
