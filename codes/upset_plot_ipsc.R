# -----------------------------------------------------------
# Script Name: upset_plot_ipsc.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

get_counts <- function(G,cond,A=rep(c(1,-1,0),3)) {
  if ((sum(cond[1:3])!=1 | sum(cond[4:6])!=1 | sum(cond[7:9])!=1)) {
    r=0
  } else {
    cond=A[which(cond==1)]
    ind=vector("list", length(cond))
    for (i in 1:length(cond)) {
      if (cond[i]==(-1)) {
        ind[[i]]=which(G[[i]]$adj.P.Val<0.05 & G[[i]]$logFC<0)
      }
      if (cond[i]==(0)) {
        ind[[i]]=which(G[[i]]$adj.P.Val>=0.05)
      }
      if (cond[i]==(1)) {
        ind[[i]]=which(G[[i]]$adj.P.Val<0.05 & G[[i]]$logFC>0)
      }
    }
    r=length(Reduce(intersect, ind))
  }
  return(r)
}

upset_plot_ipsc <- function(path,chr) {
  # set chr="X" or chr="aut"
  
  #library(ComplexUpset)
  #library(ggplot2)
  #library(tidyr)
  hipsci_datadir = file.path(path,"data")
  hipsci_resultsdir = file.path(path,"results")
  hipsci_figuresdir = file.path(path,"figures")
  plotname=file.path(hipsci_figuresdir,paste("ipsc_upset_",chr,".pdf",sep=""))
  
  l <- list(G1up = c(TRUE,FALSE), G1down = c(TRUE,FALSE), G1no = c(TRUE,FALSE),
          G2up = c(TRUE,FALSE), G2down = c(TRUE,FALSE), G2no = c(TRUE,FALSE),
          G3up = c(TRUE,FALSE), G3down = c(TRUE,FALSE), G3no = c(TRUE,FALSE))
  l=expand.grid(l)
  SS=c("G1:up","G1:down","G1:notsignificant","G2:up","G2:down","G2:notsignificant","G3:up","G3:down","G3:notsignificant")
  names(SS) <- SS
  
  f1=read.table(file.path(hipsci_resultsdir,"top.table_final_sva_dream_ipsc_MG1_FALSE_TRUE.txt"))
  f2=read.table(file.path(hipsci_resultsdir,"top.table_final_sva_dream_ipsc_MG2_FALSE_TRUE.txt"))
  f3=read.table(file.path(hipsci_resultsdir,"top.table_final_sva_dream_ipsc_MG3_FALSE_TRUE.txt"))
  f1$logFC=(-f1$logFC)
  f2$logFC=(-f2$logFC)
  f3$logFC=(-f3$logFC)
  if (chr=="X") {
    ind_chr=which(f1$chrs=="X")
  } else if (chr=="aut") {
    ind_chr=which(f1$chrs!="X" & f1$chrs!="Y" & f1$chrs!="MT")
  } else {
    ind_chr=(1:dim(f1)[1])
  }
  f1=f1[ind_chr,]
  f2=f2[ind_chr,]
  f3=f3[ind_chr,]
  
  G=list(f1,f2,f3)
  count=matrix(0,dim(l)[1],1)
  for (i in 1:dim(l)[1]) {
    cond=l[i,]
    count[i]=get_counts(G,cond)
  }
  CC=cbind(l,count)
  DD=as_tibble(CC)
  #DD %>% print(n = nrow(DD))
  indvs = DD %>% uncount(count) 
  rm.ind=which(indvs$G1no==TRUE & indvs$G2no==TRUE & indvs$G3no==TRUE)
  if (length(rm.ind)!=0) {
    indvs=indvs[-rm.ind,]
  }
  colnames(indvs)=SS
  
  if (chr=="X") {
    indvs$par=FALSE
    indvs$par[which(indvs[,1]==FALSE & indvs[,2]==FALSE & indvs[,4]==TRUE & indvs[,5]==FALSE  & indvs[,7]==TRUE & indvs[,8]==FALSE)][1]=TRUE
    indvs$par[which(indvs[,1]==FALSE & indvs[,2]==FALSE &  indvs[,4]==FALSE & indvs[,5]==FALSE  & indvs[,7]==TRUE & indvs[,8]==FALSE )][1]=TRUE
    indvs$par[which(indvs[,1]==FALSE & indvs[,2]==TRUE &  indvs[,4]==FALSE & indvs[,5]==FALSE  & indvs[,7]==FALSE & indvs[,8]==FALSE )][1:8]=TRUE
    indvs$par[which(indvs[,1]==FALSE & indvs[,2]==FALSE &  indvs[,4]==FALSE & indvs[,5]==FALSE  & indvs[,7]==FALSE & indvs[,8]==TRUE )][1]=TRUE
    indvs$par[which(indvs[,1]==FALSE & indvs[,2]==TRUE &  indvs[,4]==FALSE & indvs[,5]==TRUE  & indvs[,7]==FALSE & indvs[,8]==TRUE )][1:2]=TRUE
    indvs$par[which(indvs[,1]==FALSE & indvs[,2]==TRUE & indvs[,4]==FALSE & indvs[,5]==FALSE  & indvs[,7]==FALSE & indvs[,8]==TRUE )][1]=TRUE
    
    PAR=as.factor(indvs$par)
    p=upset(data = indvs[,c(1,2,4,5,7,8)], intersect = SS[c(1,2,4,5,7,8)],  # Remove notsignificant rows
          base_annotations=list(
            'Intersection size'=intersection_size(counts=TRUE,mapping=aes(fill=PAR)) +
              scale_fill_manual(values=c('TRUE'='gray60','FALSE'='gray26'))),
          name="", 
          min_size = 0,
          width_ratio = 0.125,
          themes=upset_modify_themes(list('intersections_matrix'=theme(axis.text.y=element_text(size=12))))) +
    labs(title = "",caption = "")
  } else {
    p=upset(data = indvs[,c(1,2,4,5,7,8)], intersect = SS[c(1,2,4,5,7,8)],  # Remove notsignificant rows
          name="", 
          min_size = 0,
          width_ratio = 0.125,
          themes=upset_modify_themes(list('intersections_matrix'=theme(axis.text.y=element_text(size=12))))) +
      labs(title = "",caption = "")
  }
  if (chr=="X") {
    ggsave(plotname, p, width=25,height=13,units="cm",limitsize = FALSE)
  } else if (chr=="aut") {
    ggsave(plotname, p, width=35,height=18,units="cm",limitsize = FALSE)
  }
}

path="~/Documents/Git/github-work/hipsci-xci-final2"
upset_plot_ipsc(path,chr="X")
upset_plot_ipsc(path,chr="aut")