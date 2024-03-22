human_tissues_vs_hipsci_comparison <- function(human_xci,ddf,th_inactive=0,th_esc=0.25) {
  
  # Take the genes where info is available in both human adult tissues and ipscs
  common_genes=intersect(extract_text_before_first_dot(names(human_xci)),extract_text_before_first_dot(ddf$gene_id))
  
  ddf=ddf[match(common_genes,extract_text_before_first_dot(ddf$gene_id)),]
  human_xci=human_xci[match(common_genes,extract_text_before_first_dot(names(human_xci)))]
  table(human_xci)/length(human_xci)
  
  # Include only genes where it is observed in at least one cell line
  ind_a=which(ddf$num_observed_high>0)
  ind_b=which(ddf$num_observed_highescape>0)
  ind_c=which(ddf$num_observed_low>0)
  
  a=ddf$num_escape_high[ind_a]/ddf$num_observed_high[ind_a]
  b=ddf$num_escape_highescape[ind_b]/ddf$num_observed_highescape[ind_b]
  c=ddf$num_escape_low[ind_c]/ddf$num_observed_low[ind_c]
  
  ths=seq(0,1,by=0.05)
  ths=ths[which(ths>th_inactive)]
  p.vals=matrix(0,length(ths),3)
  for (i in 1:length(ths)) {
    th=ths[i]
    p.vals[i,1]=chisq.test(rbind(c(sum(a>=th),sum(a<=th_inactive),sum(a<th & a>th_inactive)),table(human_xci[ind_a])))$p.value
    
    p.vals[i,2]=chisq.test(rbind(c(sum(b>=th),sum(b<=th_inactive),sum(b<th & b>th_inactive)),table(human_xci[ind_b])))$p.value
    
    p.vals[i,3]=chisq.test(rbind(c(sum(c>=th),sum(c<=th_inactive),sum(c<th & c>th_inactive)),table(human_xci[ind_c])))$p.value
  }
  
  print("Chi-squared test p-val for Group 1 vs Human XCI status proportions:")
  print(p.vals[which(ths==th_esc),])
  
  df.pvals=as.data.frame(p.vals)
  names(df.pvals)=c("Group 1","Group 2","Group 3")
  df.pv=melt(df.pvals)
  df.pv$escape_threshold=rep(ths,3)
  
  p.escthr=ggplot(df.pv,aes(x=escape_threshold,y=value,color=variable)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept=0.05) +
    xlab("Minimum fraction of cell lines\nwhere gene ASE > 0.1") +
    ylab("Chi-squared test p-value\n(group XCI vs human XCI proportions)") +
    scale_x_continuous(breaks=ths) +
    scale_y_continuous(breaks=c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
    scale_color_manual(name="Female groups", values=c("#CC0033", "mediumorchid","#FF99CC"),labels=c("Group 1","Group 2", "Group 3"), breaks=c("Group 1","Group 2", "Group 3")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
    theme(text = element_text(size=20)) 
  
  ggsave("/Users/topah/Desktop/hipsci_codes/figures/esc_threshold.pdf", p.escthr, width=20,height=18,units="cm",limitsize = FALSE)
  
  group1_xci=c("",length(a))
  group1_xci[which(a>=th_esc)]="Escape"
  group1_xci[which(a<th_esc & a>th_inactive)]="Variable"
  group1_xci[which(a<=th_inactive)]="Inactive"
  
  group2_xci=c("",length(b))
  group2_xci[which(b>=th_esc)]="Escape"
  group2_xci[which(b<th_esc & b>th_inactive)]="Variable"
  group2_xci[which(b<=th_inactive)]="Inactive"
  
  group3_xci=c("",length(c))
  group3_xci[which(c>=th_esc)]="Escape"
  group3_xci[which(c<th_esc & c>th_inactive)]="Variable"
  group3_xci[which(c<=th_inactive)]="Inactive"
  
  df_g1=data.frame(nums=c(table(group1_xci,human_xci[ind_a])), human_xci=c(rep("Escape",3),rep("Inactive",3),rep("Variable",3)),group_xci=rep(c("Escape","Inactive","Variable"),3))
  df_g1$tot=length(a)
  df_g1$group="Group1"
  
  df_g2=data.frame(nums=c(table(group2_xci,human_xci[ind_b])), human_xci=c(rep("Escape",3),rep("Inactive",3),rep("Variable",3)),group_xci=rep(c("Escape","Inactive","Variable"),3))
  df_g2$tot=length(b)
  df_g2$group="Group2"
  
  df_g3=data.frame(nums=c(table(group3_xci,human_xci[ind_c])), human_xci=c(rep("Escape",3),rep("Inactive",3),rep("Variable",3)),group_xci=rep(c("Escape","Inactive","Variable"),3))
  df_g3$tot=length(c)
  df_g3$group="Group3"
  
  df_g=rbind(df_g1,df_g2,df_g3)
  df_g$perc=100*df_g$nums/df_g$tot
  df_g$lower=100*apply(as.matrix(df_g), 1, function(x) {prop.test(as.numeric(x[1]),as.numeric(x[4]))$conf.int[1]})
  df_g$upper=100*apply(as.matrix(df_g), 1, function(x) {prop.test(as.numeric(x[1]),as.numeric(x[4]))$conf.int[2]})
  df_g$comb_group=paste(df_g$group,"_",df_g$group_xci,sep="")
  
  df_g$all_lower=0
  df_g$all_upper=0
  
  df_g$all_lower[which(df_g$comb_group=="Group1_Escape")]=sum(df_g$lower[which(df_g$comb_group=="Group1_Escape")])
  df_g$all_lower[which(df_g$comb_group=="Group2_Escape")]=sum(df_g$lower[which(df_g$comb_group=="Group2_Escape")])
  df_g$all_lower[which(df_g$comb_group=="Group3_Escape")]=sum(df_g$lower[which(df_g$comb_group=="Group3_Escape")])
  df_g$all_lower[which(df_g$comb_group=="Group1_Inactive")]=sum(df_g$lower[which(df_g$comb_group=="Group1_Inactive")])
  df_g$all_lower[which(df_g$comb_group=="Group2_Inactive")]=sum(df_g$lower[which(df_g$comb_group=="Group2_Inactive")])
  df_g$all_lower[which(df_g$comb_group=="Group3_Inactive")]=sum(df_g$lower[which(df_g$comb_group=="Group3_Inactive")])
  df_g$all_lower[which(df_g$comb_group=="Group1_Variable")]=sum(df_g$lower[which(df_g$comb_group=="Group1_Variable")])
  df_g$all_lower[which(df_g$comb_group=="Group2_Variable")]=sum(df_g$lower[which(df_g$comb_group=="Group2_Variable")])
  df_g$all_lower[which(df_g$comb_group=="Group3_Variable")]=sum(df_g$lower[which(df_g$comb_group=="Group3_Variable")])
  
  df_g$all_upper[which(df_g$comb_group=="Group1_Escape")]=sum(df_g$upper[which(df_g$comb_group=="Group1_Escape")])
  df_g$all_upper[which(df_g$comb_group=="Group2_Escape")]=sum(df_g$upper[which(df_g$comb_group=="Group2_Escape")])
  df_g$all_upper[which(df_g$comb_group=="Group3_Escape")]=sum(df_g$upper[which(df_g$comb_group=="Group3_Escape")])
  df_g$all_upper[which(df_g$comb_group=="Group1_Inactive")]=sum(df_g$upper[which(df_g$comb_group=="Group1_Inactive")])
  df_g$all_upper[which(df_g$comb_group=="Group2_Inactive")]=sum(df_g$upper[which(df_g$comb_group=="Group2_Inactive")])
  df_g$all_upper[which(df_g$comb_group=="Group3_Inactive")]=sum(df_g$upper[which(df_g$comb_group=="Group3_Inactive")])
  df_g$all_upper[which(df_g$comb_group=="Group1_Variable")]=sum(df_g$upper[which(df_g$comb_group=="Group1_Variable")])
  df_g$all_upper[which(df_g$comb_group=="Group2_Variable")]=sum(df_g$upper[which(df_g$comb_group=="Group2_Variable")])
  df_g$all_upper[which(df_g$comb_group=="Group3_Variable")]=sum(df_g$upper[which(df_g$comb_group=="Group3_Variable")])
  
  
  df_human=data.frame(nums=c(sum(human_xci=="Escape"),sum(human_xci=="Inactive"),sum(human_xci=="Variable")),human_xci=c("Escape","Inactive","Variable"),group_xci=c("Escape","Inactive","Variable"),
                      tot=c(rep(length(human_xci),3)),group=rep("Human",3))
  df_human$perc=100*df_human$nums/df_human$tot
  df_human$lower=100*apply(as.matrix(df_human), 1, function(x) {prop.test(as.numeric(x[1]),as.numeric(x[4]))$conf.int[1]})
  df_human$upper=100*apply(as.matrix(df_human), 1, function(x) {prop.test(as.numeric(x[1]),as.numeric(x[4]))$conf.int[2]})
  df_human$comb_group=paste(df_human$group,"_",df_human$group_xci,sep="")
  df_human$all_lower=df_human$lower
  df_human$all_upper=df_human$upper
  
  df_gg=rbind(df_human,df_g)
  df_gg$comb_group=factor(df_gg$comb_group,levels=c("Human_Escape","Human_Variable", "Human_Inactive",  "Group1_Escape",
                                                    "Group1_Variable", "Group1_Inactive", "Group2_Escape", "Group2_Variable",
                                                    "Group2_Inactive", "Group3_Escape", "Group3_Variable", "Group3_Inactive"))
  df_gg$group=factor(df_gg$group,levels=unique(df_gg$group))
  df_gg$human_xci=factor(df_gg$human_xci,levels=c("Escape","Variable","Inactive"))
  
  cols=c('red',"darkgoldenrod1",'blue')
  
  p.xci=ggplot(data=df_gg, aes(x=comb_group, y=perc,fill=human_xci,color=group)) +
    geom_bar(position='stack',stat="identity",aes(fill=human_xci,color=group),linewidth=0.8) +
    geom_errorbar(data = df_gg,
                  aes(x = comb_group, ymax = all_upper, ymin = all_lower), 
                  width = 0.2) +
    ylim(0,100) +
    theme_minimal() +
    scale_fill_manual(values=cols,breaks=c("Escape","Variable","Inactive"),labels=c("Escape","Variable","Inactive")) +
    scale_color_manual(name="Female\ngroups",values=c("black","#CC0033", "mediumorchid","#FF99CC"),labels=c("Human","Group 1","Group 2", "Group 3"), breaks=c("Human","Group1","Group2", "Group3")) +
    xlab("") +
    ylab("Percentage of chr-X genes") +
    labs(fill="Human\nXCI status") +
    theme(text = element_text(size=20)) +
    guides(color="none") +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, color=c(rep("black",3),rep("#CC0033",3),rep("mediumorchid",3),rep("#FF99CC",3)))) 
  
  ggsave(paste("/Users/topah/Desktop/hipsci_codes/figures/xci_",sub("\\.", "", as.character(th_esc)),".pdf",sep=""), p.xci, width=18,height=18,units="cm",limitsize = FALSE)
  
  print(paste("Percentage of human escape: ",df_gg$perc[which(df_gg$comb_group=="Human_Escape")],sep=""))
  print(paste("Percentage of human variable: ",df_gg$perc[which(df_gg$comb_group=="Human_Variable")],sep=""))
  print(paste("Percentage of human inactive: ",df_gg$perc[which(df_gg$comb_group=="Human_Inactive")],sep=""))
  
}