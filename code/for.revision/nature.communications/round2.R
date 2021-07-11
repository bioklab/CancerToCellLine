require(plyr)
require(dplyr)
require(pheatmap)
require(ggplot2)
library(ComplexHeatmap)
require(foreach)
source('client-side/code/for.figure/ggplot.style.R')
load('~/Project/Cancer2CellLine/client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.mutation.profile.RData')
load('~/Project/Cancer2CellLine/client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
load('~/Project/Cancer2CellLine/client-side/output/for.revision.nature.communications.round1.R.output//for.revision.nature.communications.round1.RData')

######################################################################################################################
#   
#  Fig 2, revised. Not used in the final submission
#  
######################################################################################################################


############################################################################################################
###### Re-generate Fig 2a, Use R Oncoprint package to show the mutation porifles of hotspot mutated and DE mutated genes  
############################################################################################################

# #### Compute differentially-mutated gene between MET500 and TCGA breast cancer samples
# 
# MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)
# common.gene.list     <- intersect(names(GDC.TCGA.brca.gene.mutation.freq),names(MET500.breast.cancer.gene.mutation.freq)) 
# dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
#   q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
#   prob    <- GDC.TCGA.brca.gene.mutation.freq[g]    
#   if(prob ==0) {
#     prob <- min(GDC.TCGA.brca.gene.mutation.freq[GDC.TCGA.brca.gene.mutation.freq > 0])  
#     #actually, the prob value equals to  1/ncol(GDC.TCGA.brca.gene.mutation.profile)  
#     
#   }
#   p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = FALSE)#right side p-value,metastatic cancer has more mutation burdens
#   data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
# }
# rownames(dm.df)  <- common.gene.list
# dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
# dm.df            <- dm.df[order(dm.df$fdr),]
# dm.gene          <- rownames(dm.df)[dm.df$fdr  < 0.001 ]
# MET500.vs.GDC.dm.df <- dm.df











#### Check the overlap between highly mutated and differntially mutated genes, and pool them together
library(gplots)
venn(list(highly.mutated.gene=MET500.breast.cancer.top.mutated.gene,differentialy.mutated.gene=dm.gene))
pooled.gene              <- c(dm.gene,setdiff(MET500.breast.cancer.top.mutated.gene,dm.gene)) %>% unique

TCGA.mutation.cnt        <- apply(TCGA.breast.cancer.gene.mutation.profile[pooled.gene,],1,sum) 
TCGA.mutation.frequency  <- apply(TCGA.breast.cancer.gene.mutation.profile[pooled.gene,],1,sum) / ncol(TCGA.breast.cancer.gene.mutation.profile)

MET500.mutation.cnt        <- apply(MET500.breast.cancer.gene.mutation.profile[pooled.gene,],1,sum) 
MET500.mutation.frequency  <- apply(MET500.breast.cancer.gene.mutation.profile[pooled.gene,],1,sum) / ncol(MET500.breast.cancer.gene.mutation.profile)


CCLE.mutation.cnt        <- apply(CCLE.breast.cancer.gene.mutation.profile[pooled.gene,],1,sum) 
CCLE.mutation.frequency  <- apply(CCLE.breast.cancer.gene.mutation.profile[pooled.gene,],1,sum) / ncol(CCLE.breast.cancer.gene.mutation.profile)


print(names(CCLE.mutation.frequency)[CCLE.mutation.frequency == 0]) # genes which ARE NOT mutated in any cell lines


### Rearrange MET500 and TCGA samples according to clustering results of mutation profile,rearrange CCLE samples according to sites (metastaic or primary) 
MET500.dist              <- dist(MET500.breast.cancer.gene.mutation.profile[pooled.gene,] %>% t,method='binary')
MET500.hclust.rs         <- hclust(MET500.dist)
MET500.rearranged.sample <- MET500.hclust.rs$labels[MET500.hclust.rs$order]

TCGA.dist                <- dist(TCGA.breast.cancer.gene.mutation.profile[pooled.gene,] %>% t,method='binary')
TCGA.hclust.rs           <- hclust(TCGA.dist)
TCGA.rearranged.sample   <- TCGA.hclust.rs$labels[TCGA.hclust.rs$order]

CCLE.rearranged.sample   <- c(CCLE.breast.cancer.metastatic.cell.line,CCLE.breast.cancer.non.metastatic.cell.line)


### Rearrange pooled gene according to clustering results
pooled.gene.dist         <- dist(cbind(MET500.breast.cancer.gene.mutation.freq[pooled.gene],CCLE.breast.cancer.gene.mutation.freq[pooled.gene]),method='euclidean')
pooled.gene.rs           <- hclust(pooled.gene.dist)
pooled.gene              <- pooled.gene.rs$labels[pooled.gene.rs$order]


#### Assemble the mutation profile matrix for the pooled genes
pooled.mutation.profile      <- cbind(MET500.breast.cancer.gene.mutation.profile[pooled.gene,MET500.rearranged.sample],
                                      CCLE.breast.cancer.gene.mutation.profile[pooled.gene,  CCLE.rearranged.sample],
                                      TCGA.breast.cancer.gene.mutation.profile[pooled.gene,  TCGA.rearranged.sample]
)

#### Show the oncoprint and save it
source.vec <- c( rep(times= MET500.rearranged.sample %>% length, x='MET500'),
                 rep(times= CCLE.rearranged.sample   %>% length, x='CCLE'),
                 rep(times= TCGA.rearranged.sample   %>% length, x='TCGA')
                 
)
names(source.vec) <- c(MET500.rearranged.sample,CCLE.rearranged.sample,TCGA.rearranged.sample)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h, gp = gpar(fill = "#CCCCCC", col = NA))
  },
  
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.8, gp = gpar(fill = "black", col = NA))
  }
)
col <- c("MUT" = "black",'background' = '#CCCCCC')

library(circlize)
pdf("/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/oncoprint.MET500.CCLE.revised.pdf",width = 30,height=15 )
col.annotation  <-  HeatmapAnnotation(type=source.vec[c(MET500.rearranged.sample,CCLE.rearranged.sample)],which='column',col = list(type = c("MET500" =  "#800000", "TCGA" = "#008080","CCLE" = "#e6beff")) ,show_legend = FALSE)
oncoPrint(pooled.mutation.profile[,c(MET500.rearranged.sample,CCLE.rearranged.sample)], get_type = function(x) {ifelse(x == 1,'MUT','background')},row_names_side = 'left',
          show_heatmap_legend = FALSE,
          alter_fun = alter_fun, col = col, row_order = NULL, column_order=NULL, remove_empty_columns = T,
          show_column_names = FALSE,show_row_barplot = F,top_annotation = col.annotation,show_pct = FALSE,row_names_gp = gpar(fontsize = 13,fontface=2)
) + 
  Heatmap(cbind(MET500.breast.cancer.gene.mutation.freq[pooled.gene],CCLE.breast.cancer.gene.mutation.freq[pooled.gene]),show_column_dend = FALSE, 
          col=colorRamp2(c(0,1),c('blue','red')),
          row_names_gp = gpar(fontsize = 13,fontface=2),
          width = unit(4, "cm"),show_heatmap_legend = FALSE,
          top_annotation=HeatmapAnnotation( show_legend = FALSE,type=c('MET500','CCLE'),which='column',col = list(type = c("MET500" =  "#800000", "TCGA" = "#008080","CCLE" = "#e6beff")) )
  )
dev.off()



######################################################################################################################
#   
#  Fig S1, revised. NOT used in the final revision
#  
######################################################################################################################


###########################################################
###### Fig S1b revised, draw the valcano plot of differential mutation analysis between MET500 and TCGA
###########################################################
dm.df$is.de.gene <- ifelse(rownames(dm.df) %in% dm.gene,'Y','N')
draw.df          <- dm.df
ggplot(draw.df,aes(x=log2(ratio+1),y=-1 * log10(fdr),col=is.de.gene)) + geom_point(size=6.0,show.legend=F) + 
  xlab('log2((MET500/TCGA)+1)') + ylab('-log10(fdr)') + xlim(0,6.5) + geom_hline(yintercept = 3,linetype="dashed") + 
  ggplot.style +scale_color_manual(values=c('Y' = 'red','N' = 'grey')) 
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/MET500.vs.TCGA.differential.mutation.frequency.volcano.plot.revised.pdf',width=20,height=10)


###########################################################
###### Fig S1c revised, draw the mutation frequency of cell line hypermutated genes
###########################################################

CCLE.highly.mutated.gene           <-  names(CCLE.breast.cancer.gene.mutation.freq)[CCLE.breast.cancer.gene.mutation.freq >= 0.5]
CCLE.specific.highly.mutated.gene  <-  setdiff(CCLE.highly.mutated.gene,pooled.gene)
col.annotation                     <-  HeatmapAnnotation(type=c('MET500','TCGA','CCLE'),which='column',col = list(type = c("MET500" =  "#800000", "TCGA" = "#008080","CCLE" = "#e6beff")) ,show_legend = FALSE)
pooled.gene                        <- c(dm.gene,setdiff(MET500.breast.cancer.top.mutated.gene,dm.gene)) %>% unique

require(circlize)
pdf("/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/CCLE.specific.highly.mutated.gene.mutation.frequency.heatmap.revised.pdf",width = 30,height=15 )
Heatmap(log10(0.01+cbind(MET500.breast.cancer.gene.mutation.freq[CCLE.specific.highly.mutated.gene],GDC.TCGA.brca.gene.mutation.freq[CCLE.specific.highly.mutated.gene],CCLE.breast.cancer.gene.mutation.freq[CCLE.specific.highly.mutated.gene]) ),show_heatmap_legend = FALSE, 
        top_annotation = col.annotation,cluster_rows = FALSE,cluster_columns = FALSE,
        col=colorRamp2(c(-2,-0.05),c('blue','red')),row_names_gp = gpar(fontsize = 35,fontface=2),
        width = unit(20, "in")
)
dev.off()

###########################################################
######Fig S11 c and d NOT USED in the final revision
###########################################################

## re-download the most recent cBioPortal data!!!!! Well, not used in the final revision, just to check consitency with 
# data downloaded in April 2018

# study.id   <- 'brca_tcga'
# mycgds     <-  CGDS("http://www.cbioportal.org/")
# 
# case.list  <- 'brca_tcga_sequenced'
# profile.id <- 'brca_tcga_mutations'
# clinical.data                           <- getClinicalData(mycgds,case.list)
# Breast.Invasive.Ductal.Carcinoma.sample <- rownames(clinical.data)[clinical.data$CANCER_TYPE_DETAILED == 'Breast Invasive Ductal Carcinoma']
# 
# 
# gene.list <- rownames(GDC.TCGA.brca.gene.mutation.profile)
# gene.list[gene.list == 'NOV'] <- 'CCN3' # Gene symbol is invalid in cBioPortal, replace it with its alias CCN3
# tmp3 <- foreach(i = seq(from = 1,to=length(gene.list),by=10)) %do% {
#   low <- i
#   up  <- i + 9
#   if(up > length(gene.list) ) {
#     up <- length(gene.list)
#   }
#   pool.gene <- gene.list[low:up]
#   tmp      <- getProfileData(mycgds,pool.gene,profile.id,case.list)
#   tmp      <- tmp[rownames(tmp) %in% Breast.Invasive.Ductal.Carcinoma.sample,]
#   TCGA.breast.cancer.gene.mutation.profile <- as.matrix(tmp) 
#   TCGA.breast.cancer.gene.mutation.profile[is.na(TCGA.breast.cancer.gene.mutation.profile)]    <- '0'
#   TCGA.breast.cancer.gene.mutation.profile[TCGA.breast.cancer.gene.mutation.profile == 'NaN' ] <- '0'
#   TCGA.breast.cancer.gene.mutation.profile[TCGA.breast.cancer.gene.mutation.profile != '0']    <- '1'
#   tmp2 <- matrix(data=as.integer(TCGA.breast.cancer.gene.mutation.profile), 
#                  ncol=ncol(TCGA.breast.cancer.gene.mutation.profile), 
#                  nrow=nrow(TCGA.breast.cancer.gene.mutation.profile)
#   )
#   rownames(tmp2) <- rownames(TCGA.breast.cancer.gene.mutation.profile)
#   colnames(tmp2) <- colnames(TCGA.breast.cancer.gene.mutation.profile)
#   tmp2
# }
# 
# row.name <- rownames(tmp3[[1]])
# tmp4 <- foreach(l = tmp3,.combine='cbind') %do% {
#   l[row.name,]  
# }
# 
# 
# TCGA.breast.cancer.gene.mutation.profile <- tmp4 %>% t
# flag                                     <- rownames(TCGA.breast.cancer.gene.mutation.profile) == 'NKX2.1'
# rownames(TCGA.breast.cancer.gene.mutation.profile)[flag] <- 'NKX2-1' #hmm, cBioPortal returns NKX2-1 for gene NKX2.1, replace back
# flag                                     <- rownames(TCGA.breast.cancer.gene.mutation.profile) == 'MRTFA'
# rownames(TCGA.breast.cancer.gene.mutation.profile)[flag] <- 'MKL1' #hmm, cBioPortal returns MRTFA for gene MKL1, replace back
# flag                                     <- rownames(TCGA.breast.cancer.gene.mutation.profile) == 'LTO1'
# rownames(TCGA.breast.cancer.gene.mutation.profile)[flag] <- 'ORAOV1' #hmm, cBioPortal returns LTO1 for gene ORAOV1, replace back
# flag                                     <- rownames(TCGA.breast.cancer.gene.mutation.profile) == 'CCN3'
# rownames(TCGA.breast.cancer.gene.mutation.profile)[flag] <- 'NOV' #hmm, cBioPortal returns CCN3 for gene NOV, replace back
# 
# 
# TCGA.breast.cancer.gene.mutation.freq    <- apply(TCGA.breast.cancer.gene.mutation.profile,1,function(x) sum(x)/ncol(TCGA.breast.cancer.gene.mutation.profile)) %>% sort(decreasing = TRUE)
# TCGA.breast.cancer.top.mutated.gene      <- names(TCGA.breast.cancer.gene.mutation.freq)[TCGA.breast.cancer.gene.mutation.freq >= 0.05]
# 
# cBioPortal.TCGA.brca.gene.mutation.freq    <- TCGA.breast.cancer.gene.mutation.freq
# cBioPortal.TCGA.brca.gene.mutation.profile <- TCGA.breast.cancer.gene.mutation.profile
# 
# 
# MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)
# common.gene.list     <- intersect(names(cBioPortal.TCGA.brca.gene.mutation.freq),names(MET500.breast.cancer.gene.mutation.freq)) # Well, this operation only removes gene ZNF668, whose mutation data is not returned by cBioportal
# dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
#   q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
#   prob    <- cBioPortal.TCGA.brca.gene.mutation.freq[g]    
#   if(prob ==0) {
#     prob <- min(cBioPortal.TCGA.brca.gene.mutation.freq[cBioPortal.TCGA.brca.gene.mutation.freq > 0])  
#   }
#   p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F)#right side p-value,metastatic cancer has more mutation burdens
#   data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
# }
# rownames(dm.df)  <- common.gene.list
# dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
# dm.df            <- dm.df[order(dm.df$fdr),]
# dm.gene          <- rownames(dm.df)[dm.df$fdr  < 0.001 ]
# MET500.vs.cBioPortal.dm.df <- dm.df
# 
# save(file='client-side/output/for.revision.nature.communications.round2.R.output/for.revision.nature.communications.round2.RData',list=c('cBioPortal.TCGA.brca.gene.mutation.freq','cBioPortal.TCGA.brca.gene.mutation.profile'))




###############################################################
######### show the mutation frequency are consistent, Fig S11, panel c. NOT used in final submission
###############################################################
plot(x=cBioPortal.TCGA.brca.gene.mutation.freq,GDC.TCGA.brca.gene.mutation.freq[names(cBioPortal.TCGA.brca.gene.mutation.freq)],xlab='cBioPortal',ylab='GDC',cex.lab=1.5,cex.axis=1.5)
draw.df <- data.frame(cBioPortal = cBioPortal.TCGA.brca.gene.mutation.freq,GDC = GDC.TCGA.brca.gene.mutation.freq[names(cBioPortal.TCGA.brca.gene.mutation.freq)])
ggplot(draw.df,aes(x=cBioPortal,y=GDC)) + geom_point(size=5) + ggplot.style + xlim(c(0,0.5)) + ylim(c(0,0.5)) + geom_abline(slope = 1,intercept = 0)
cor(draw.df$cBioPortal,draw.df$GDC,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.2/cBioportal.vs.GDC.mutation.frequency.pdf',width=20,height=20)


###############################################################
######### show the p-values are consistent, Fig S11, panel d. NOT used in final submission
###############################################################
draw.df <- data.frame(cBioPortal=-1 * (MET500.vs.cBioPortal.dm.df$p.value %>% log10),GDC= -1 * (MET500.vs.GDC.dm.df[rownames(MET500.vs.cBioPortal.dm.df),'p.value'] %>% log10))
ggplot(draw.df,aes(x=cBioPortal,y=GDC)) + geom_point(size=5) + ggplot.style + xlim(c(0,13)) + ylim(c(0,13)) + geom_abline(slope = 1,intercept = 0) + xlab('MET500.vs.cBioPortal') + ylab('MET500.vs.GDC')
cor(draw.df$cBioPortal,draw.df$GDC,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.2/p.value.consistency.pdf',width=20,height=20)





##############################################################################################################################
####### Below, I addressed the issues proposed by reviewer 3 in the second round of review
##############################################################################################################################

######################################################################
#
# Question 1: No computation is needed here
#
######################################################################




######################################################################
#
# Question 2:  MET500.vs.Stage.IIA  MET500.vs.Stage.IIB MET500.vs.StageIIIA, results highly correlated, demonstrating that stage is not a severe confounding factor
# Subtype-specific analysis, hard to say for Basal-like due to small number of samples.
#
######################################################################
# re-compute differential mutation with the TCGA mutation data downloaded from April
MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)
common.gene.list     <- intersect(names(TCGA.breast.cancer.gene.mutation.freq),names(MET500.breast.cancer.gene.mutation.freq)) # Well, this operation only removes gene ZNF668, whose mutation data is not returned by cBioportal
dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
  q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
  prob    <- TCGA.breast.cancer.gene.mutation.freq[g]    
  if(prob ==0) {
    prob <- min(TCGA.breast.cancer.gene.mutation.freq[TCGA.breast.cancer.gene.mutation.freq > 0])  
  }
  p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F)#right side p-value,metastatic cancer has more mutation burdens
  data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
dm.gene          <- rownames(dm.df)[dm.df$fdr  < 0.001 ]
MET500.vs.TCGA.raw.dm.df <- dm.df


load('~/Project/Cancer2CellLine/client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.mutation.profile.RData')
require(cgdsr)
study.id   <- 'brca_tcga'
mycgds     <-  CGDS("http://www.cbioportal.org/")

case.list  <- 'brca_tcga_sequenced'
profile.id <- 'brca_tcga_mutations'
clinical.data                           <- getClinicalData(mycgds,case.list)
Breast.Invasive.Ductal.Carcinoma.sample <- rownames(clinical.data)[clinical.data$CANCER_TYPE_DETAILED == 'Breast Invasive Ductal Carcinoma']
stage.vec                               <- clinical.data[Breast.Invasive.Ductal.Carcinoma.sample,'AJCC_PATHOLOGIC_TUMOR_STAGE']
names(stage.vec)                        <- Breast.Invasive.Ductal.Carcinoma.sample

stage.mutation.freq.list <- foreach(stage = unique(stage.vec)) %do% {
  sample.id <- names(stage.vec)[stage.vec == stage]
  #gene.mutation.freq <- apply(GDC.TCGA.brca.gene.mutation.profile[,sample.id],1,function(x) sum(x)/length(x))
  gene.mutation.freq <- apply(TCGA.breast.cancer.gene.mutation.profile[,sample.id],1,function(x) sum(x)/length(x))
  
  gene.mutation.freq
}
names(stage.mutation.freq.list) <- unique(stage.vec)

######## MET500.vs.Stage.II.A TCGA samples
common.gene.list     <- intersect(names(TCGA.breast.cancer.gene.mutation.freq),names(MET500.breast.cancer.gene.mutation.freq)) 

Stage.IIA.mutation.freq <-  stage.mutation.freq.list[['Stage IIA']] 
MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)

dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {  
  q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
  prob    <- Stage.IIA.mutation.freq[g]    
  if(prob ==0) {
    prob <- min(Stage.IIA.mutation.freq[Stage.IIA.mutation.freq > 0])  
  }
  p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F)#right side p-value,metastatic cancer has more mutation burdens
  data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
dm.gene          <- rownames(dm.df)[dm.df$fdr  < 0.001 ]
MET500.vs.Stage.IIA.dm.df <- dm.df

#plot(x=-1 * MET500.vs.Stage.IIA.dm.df$p.value %>% log10,y=-1 * MET500.vs.GDC.dm.df[rownames(MET500.vs.Stage.IIA.dm.df),'p.value'] %>% log10,xlim=c(0,11),ylim=c(0,11))
plot(x=-1 * MET500.vs.Stage.IIA.dm.df$p.value %>% log10,y=-1 * MET500.vs.TCGA.raw.dm.df[rownames(MET500.vs.Stage.IIA.dm.df),'p.value'] %>% log10,xlim=c(0,11),ylim=c(0,11))
lines(c(0,11),c(0,11))
draw.df <- data.frame( x =-1 * (MET500.vs.Stage.IIA.dm.df$p.value %>% log10),y= -1 * (MET500.vs.TCGA.raw.dm.df[rownames(MET500.vs.Stage.IIA.dm.df),'p.value'] %>% log10))
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=5) + ggplot.style + xlim(c(0,13)) + ylim(c(0,13)) + geom_abline(slope = 1,intercept = 0) + xlab('MET500.vs.TCGA (Stage IIA)') + ylab('MET500.vs.TCGA')
cor(draw.df$x,draw.df$y,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.2/p.value.StageIIA.pdf',width=20,height=20)







Stage.IIB.mutation.freq <-  stage.mutation.freq.list[['Stage IIB']] 
MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)
#common.gene.list     <- intersect(names(cBioPortal.TCGA.brca.gene.mutation.freq),names(MET500.breast.cancer.gene.mutation.freq)) 
dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
  q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
  prob    <- Stage.IIB.mutation.freq[g]    
  if(prob ==0) {
    prob <- min(Stage.IIB.mutation.freq[Stage.IIB.mutation.freq > 0])  
  }
  p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F) #right side p-value,metastatic cancer has more mutation burdens
  data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
dm.gene          <- rownames(dm.df)[dm.df$fdr  < 0.001 ]
MET500.vs.Stage.IIB.dm.df <- dm.df
plot(x=-1 * MET500.vs.Stage.IIB.dm.df$p.value %>% log10,y=-1 * MET500.vs.TCGA.raw.dm.df[rownames(MET500.vs.Stage.IIB.dm.df),'p.value'] %>% log10,xlim=c(0,11),ylim=c(0,11))
lines(c(0,11),c(0,11))
draw.df <- data.frame( x =-1 * (MET500.vs.Stage.IIB.dm.df$p.value %>% log10),y= -1 * (MET500.vs.TCGA.raw.dm.df[rownames(MET500.vs.Stage.IIB.dm.df),'p.value'] %>% log10))
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=5) + ggplot.style + xlim(c(0,13)) + ylim(c(0,13)) + geom_abline(slope = 1,intercept = 0) + xlab('MET500.vs.TCGA (Stage IIB)') + ylab('MET500.vs.TCGA')
cor(draw.df$x,draw.df$y,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.2/p.value.StageIIB.pdf',width=20,height=20)




Stage.IIIA.mutation.freq <-  stage.mutation.freq.list[['Stage IIIA']] 
MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)
#common.gene.list     <- intersect(names(cBioPortal.TCGA.brca.gene.mutation.freq),names(MET500.breast.cancer.gene.mutation.freq)) 
dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
  q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
  prob    <- Stage.IIIA.mutation.freq[g]    
  if(prob ==0) {
    prob <- min(Stage.IIIA.mutation.freq[Stage.IIIA.mutation.freq > 0])  
  }
  p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F)#right side p-value,metastatic cancer has more mutation burdens
  data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
dm.gene          <- rownames(dm.df)[dm.df$fdr  < 0.001 ]
MET500.vs.Stage.IIIA.dm.df <- dm.df
plot(x=-1 * MET500.vs.Stage.IIIA.dm.df$p.value %>% log10,y=-1 * MET500.vs.TCGA.raw.dm.df[rownames(MET500.vs.Stage.IIIA.dm.df),'p.value'] %>% log10,xlim=c(0,11),ylim=c(0,11))
lines(c(0,11),c(0,11))
draw.df <- data.frame( x =-1 * (MET500.vs.Stage.IIIA.dm.df$p.value %>% log10),y= -1 * (MET500.vs.TCGA.raw.dm.df[rownames(MET500.vs.Stage.IIIA.dm.df),'p.value'] %>% log10))
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=5) + ggplot.style + xlim(c(0,13)) + ylim(c(0,13)) + geom_abline(slope = 1,intercept = 0) + xlab('MET500.vs.TCGA (Stage IIIA)') + ylab('MET500.vs.TCGA')
cor(draw.df$x,draw.df$y,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.2/p.value.StageIIIA.pdf',width=20,height=20)





# Stage.IIA.mutation.freq  <-  stage.mutation.freq.list[['Stage IIA']] 
# Stage.IA.mutation.freq   <-  stage.mutation.freq.list[['Stage IA']] 
# dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
#   q       <- GDC.TCGA.brca.gene.mutation.profile[g,names(stage.vec)[stage.vec == 'Stage IA']] %>% sum
#   prob    <- Stage.IIA.mutation.freq[g]    
#   if(prob ==0) {
#     prob <- min(Stage.IIA.mutation.freq[Stage.IIA.mutation.freq > 0])  
#   }
#   p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = T)# left side p-value,lower stage has lower mutation burdens
#   data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
# }
# rownames(dm.df)  <- common.gene.list
# dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
# dm.df            <- dm.df[order(dm.df$fdr),]
# dm.gene          <- rownames(dm.df)[dm.df$fdr  < 0.001 ]
# Stage.IA.vs.IIA.dm.df <- dm.df
# plot(x=-1 * dm.df$p.value %>% log10,y=-1 * MET500.vs.GDC.dm.df[rownames(dm.df),'p.value'] %>% log10)


load('server-side/RData/TCGA.breast.cancer.RData')
load('server-side/RData/MET500.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')

require(genefu)
pam50.gene.df  <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#')
pam50.gene     <- pam50.gene.df$ensemble.gene.id %>% as.character
combined.data               <- cbind(MET500.log2.fpkm.matrix[pam50.gene,MET500.breast.cancer.polyA.sample],TCGA.breast.cancer.log2.fpkm.matrix[pam50.gene,])
pam50.gene.df               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#')
colnames(pam50.gene.df)[1]  <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]  <- 'probe' #damn it, genefu package doesnot mentioning this!
pam50.gene.df$EntrezGene.ID <- as.character(pam50.gene.df$EntrezGene.ID)
pam50.gene.expr             <- combined.data %>% t 
annot.matrix                <- pam50.gene.df[,1:2] %>% as.matrix
rownames(annot.matrix)      <- annot.matrix[,'probe']
pam50.subtype.rs            <- intrinsic.cluster.predict(sbt.model = pam50.robust,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix ,do.mapping = TRUE )


pam50.subtype.vec        <- pam50.subtype.rs[['subtype']]
names(pam50.subtype.vec) <- gsub(pattern = '\\-',replacement = '\\.',x = names(pam50.subtype.vec))

flag             <- pam50.subtype.vec != 'Basal' & (names(pam50.subtype.vec) %in% Breast.Invasive.Ductal.Carcinoma.sample)
TCGA.non.basal.sample <- names(pam50.subtype.vec)[flag]

flag             <- pam50.subtype.vec == 'Basal' & (names(pam50.subtype.vec) %in% Breast.Invasive.Ductal.Carcinoma.sample)
TCGA.basal.sample <- names(pam50.subtype.vec)[flag]

GDC.TCGA.brca.gene.mutation.freq.non.basal.sample <- apply(GDC.TCGA.brca.gene.mutation.profile[,TCGA.non.basal.sample],1,function(x) sum(x)/length(x))
GDC.TCGA.brca.gene.mutation.freq.basal.sample     <- apply(GDC.TCGA.brca.gene.mutation.profile[,TCGA.basal.sample],    1,function(x) sum(x)/length(x))




MET500.run.meta <- read.csv('client-side/meta.data/MET500.RNASeq.meta.csv',header=TRUE,stringsAsFactors=FALSE)

flag                    <- pam50.subtype.vec != 'Basal' & (grepl(x=names(pam50.subtype.vec),pattern = 'SRR'))
MET500.non.basal.sample <- names(pam50.subtype.vec)[flag]
idx                     <- match(MET500.non.basal.sample,table = MET500.run.meta$Run)
MET500.non.basal.cohort <- MET500.run.meta$MET500.id[idx] %>% unique
MET500.non.basal.cohort <- intersect(MET500.non.basal.cohort,colnames(MET500.breast.cancer.gene.mutation.profile))


flag                    <- pam50.subtype.vec == 'Basal' & (grepl(x=names(pam50.subtype.vec),pattern = 'SRR'))
MET500.basal.sample     <- names(pam50.subtype.vec)[flag]
idx                     <- match(MET500.basal.sample,table = MET500.run.meta$Run)
MET500.basal.cohort     <- MET500.run.meta$MET500.id[idx] %>% unique
MET500.basal.cohort     <- intersect(MET500.basal.cohort,colnames(MET500.breast.cancer.gene.mutation.profile))


MET500.sample.number <- length(MET500.non.basal.cohort)
dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
  q       <- MET500.breast.cancer.gene.mutation.profile[g,MET500.non.basal.cohort] %>% sum
  prob    <- GDC.TCGA.brca.gene.mutation.freq.non.basal.sample[g]    
  if(prob ==0) {
    prob <- min(GDC.TCGA.brca.gene.mutation.freq.non.basal.sample[GDC.TCGA.brca.gene.mutation.freq.non.basal.sample > 0])  
  }
  p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F)# right side p-value, metastatic samples, higher mutation burden
  data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
MET500.vs.GDC.non.basal.dm.df <- dm.df


MET500.sample.number <- length(MET500.basal.cohort)
dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
  q       <- MET500.breast.cancer.gene.mutation.profile[g,MET500.basal.cohort] %>% sum
  prob    <- GDC.TCGA.brca.gene.mutation.freq.basal.sample[g]    
  if(prob ==0) {
    prob <- min(GDC.TCGA.brca.gene.mutation.freq.basal.sample[GDC.TCGA.brca.gene.mutation.freq.basal.sample > 0])  
  }
  p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F)# right side p-value, metastatic samples, higher mutation burden
  data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
MET500.vs.GDC.basal.dm.df <- dm.df


draw.df <- data.frame( x =-1 * (MET500.vs.GDC.non.basal.dm.df$p.value %>% log10),y= -1 * (MET500.vs.GDC.basal.dm.df[rownames(MET500.vs.GDC.non.basal.dm.df),'p.value'] %>% log10))
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=5) + ggplot.style + xlim(c(0,13)) + ylim(c(0,13)) + geom_abline(slope = 1,intercept = 0) + xlab('MET500.vs.TCGA (non-Basal)') + ylab('MET500.vs.TCGA (Basal)')
cor(draw.df$x,draw.df$y,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.2/p.value.basal.vs.non.basal.pdf',width=20,height=20)


source('client-side/code/for.revision/nature.communicaionts/round2.DE.R') # perform DE analysis with RUVg inferred factor values
