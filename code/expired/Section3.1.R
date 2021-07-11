load('server-side/RData/CCLE.RData')
load('server-side/RData/MET500.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
source('client-side/code/util.R')
require(plyr)
require(dplyr)
require(foreach)



# Organize organoid RNASeq data
BreastCancerBioBank_RNAseq_readCounts_RPKM           <- read.delim("~/Project/Cancer2CellLine/client-side/Data/organoid.expr/BCBB_RNAseq_SummaryLevelData/BreastCancerBioBank_RNAseq_readCounts_RPKM.txt", stringsAsFactors=FALSE)
rownames(BreastCancerBioBank_RNAseq_readCounts_RPKM) <- BreastCancerBioBank_RNAseq_readCounts_RPKM$X
BreastCancerBioBank_RNAseq_readCounts_RPKM$X         <- NULL
organoid.rpkm.matrix                                 <- as.matrix(BreastCancerBioBank_RNAseq_readCounts_RPKM)

remove.version <- function(x) {
    tmp <- strsplit(x=x,split = "\\.") %>% unlist 
    tmp[1]
}
organoid.line                   <- sapply(colnames(organoid.rpkm.matrix),remove.version)
organoid.line[25]               <- 'W894.10B01A0' #hmm, two replicates for W894?
organoid.line                   <- paste(organoid.line,'_BREAST_ORGANOID',sep='')
colnames(organoid.rpkm.matrix)  <- organoid.line
organoid.log2.rpkm.matrix       <- log2(organoid.rpkm.matrix+1)



# Pick out most varied genes across CCLE cell lines
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.rna.seq.marker.gene.1000                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]


# Correlating cell lines and organoids to MET500 samples (different subtypes) 
s                   <- intersect(rownames(CCLE.log2.rpkm.matrix),rownames(organoid.rpkm.matrix))
merged.model.matrix <- cbind(CCLE.log2.rpkm.matrix[s,],organoid.log2.rpkm.matrix[s,])

MET500.vs.CCLE.and.organoid.polyA.LumB.expression.correlation.result  <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.LumB.sample],  expr.of.cell.lines = merged.model.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.Basal.sample], expr.of.cell.lines = merged.model.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.and.organoid.polyA.Her2.expression.correlation.result  <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.Her2.sample],  expr.of.cell.lines = merged.model.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.and.organoid.polyA.LumA.expression.correlation.result  <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.LumA.sample],  expr.of.cell.lines = merged.model.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.and.organoid.polyA.non.Basal.expression.correlation.result  <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,c(MET500.breast.cancer.polyA.LumA.sample,MET500.breast.cancer.polyA.LumB.sample,MET500.breast.cancer.polyA.Her2.sample)],  expr.of.cell.lines = merged.model.matrix,CCLE.rna.seq.marker.gene.1000)


# Determine subtypes of organoids
require(genefu)
pam50.gene.df               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#')
colnames(pam50.gene.df)[1]  <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]  <- 'probe' #damn it, genefu package doesnot mentioning this!
pam50.gene.df$EntrezGene.ID <- as.character(pam50.gene.df$EntrezGene.ID)
pam50.gene.expr             <- organoid.rpkm.matrix[pam50.gene.df$probe %>% as.character,] %>% t 
annot.matrix                <- pam50.gene.df[,1:2] %>% as.matrix
rownames(annot.matrix)      <- annot.matrix[,'probe']
pam50.subtype.rs            <- intrinsic.cluster.predict(sbt.model = pam50.robust,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix ,do.mapping = TRUE )

LumA.organoid               <- colnames(organoid.rpkm.matrix)[pam50.subtype.rs$subtype == 'LumA']
LumB.organoid               <- colnames(organoid.rpkm.matrix)[pam50.subtype.rs$subtype == 'LumB']
Her2.organoid               <- colnames(organoid.rpkm.matrix)[pam50.subtype.rs$subtype == 'Her2']
Basal.organoid              <- colnames(organoid.rpkm.matrix)[pam50.subtype.rs$subtype == 'Basal']



### Panel (a), visualize of organoid ranking results between different subtypes 
require('PerformanceAnalytics')
source('client-side/code/for.figure/chart.Correlation.R') # I steal the code from PerformanceAnalytics package and modify it,remove the weired stars

mat <- cbind( MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result$cell.line.median.cor[organoid.line], 
              MET500.vs.CCLE.and.organoid.polyA.LumB.expression.correlation.result$cell.line.median.cor[organoid.line],
              MET500.vs.CCLE.and.organoid.polyA.LumA.expression.correlation.result$cell.line.median.cor[organoid.line],
              MET500.vs.CCLE.and.organoid.polyA.Her2.expression.correlation.result$cell.line.median.cor[organoid.line]
)
colnames(mat) <- c('Basal-like','LuminalB','LuminalA','Her2\nenriched')

pdf(file  = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3.1/MET500.result.correlation.between.subtype.pdf',width=20,height=20)
chart.Correlation(mat, histogram=FALSE, pch=19,method = 'spearman',cex.labels=10)
dev.off()










df1 <- data.frame(median.cor = MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result$cell.line.median.cor[Basal.organoid],
                  type = rep(x='organoid',times = Basal.organoid %>% length)
                  )

df2 <- data.frame(median.cor = MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.Basal.cell.line],
                  type = rep(x='cell.line',times = CCLE.breast.cancer.Basal.cell.line %>% length)
)
require(ggplot2)
ggplot(rbind(df1,df2)) + geom_boxplot(aes(x=type,y=median.cor))




non.Basal.organoid <- c(LumA.organoid,LumB.organoid,Her2.organoid)
CCLE.breast.cancer.non.Basal.cell.line <- c(CCLE.breast.cancer.LumA.cell.line,CCLE.breast.cancer.LumB.cell.line,CCLE.breast.cancer.Her2.cell.line)

df1 <- data.frame(median.cor = MET500.vs.CCLE.and.organoid.polyA.non.Basal.expression.correlation.result$cell.line.median.cor[non.Basal.organoid],
                  type = rep(x='organoid',times = non.Basal.organoid %>% length)
)

df2 <- data.frame(median.cor = MET500.vs.CCLE.and.organoid.polyA.non.Basal.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.non.Basal.cell.line],
                  type = rep(x='cell.line',times = CCLE.breast.cancer.non.Basal.cell.line %>% length)
)
require(ggplot2)
ggplot(rbind(df1,df2)) + geom_boxplot(aes(x=type,y=median.cor))




# For each of the subtype, show that most correlated organoid is better than most correlated cell line
#LumA
cell.line <- 'MDAMB415_BREAST'
organoid  <- 'MMC01031_BREAST_ORGANOID'
rs        <- MET500.vs.CCLE.and.organoid.polyA.LumA.expression.correlation.result

df1 <- data.frame(correlation = rs$correlation.matrix[,organoid],
                  type = rep(x='MMC01031',times = nrow(rs$correlation.matrix))
)

df2 <- data.frame(correlation = rs$correlation.matrix[,cell.line],
                  type = rep(x='MDAMB415',times = nrow(rs$correlation.matrix))
)
require(ggplot2)
ggplot(rbind(df1,df2)) + geom_boxplot(aes(x=type,y=correlation))
wilcox.test(df1$correlation,df2$correlation,paired = TRUE)





#LumB
cell.line <- 'BT483_BREAST'
organoid  <- 'MMC01031_BREAST_ORGANOID'
rs        <- MET500.vs.CCLE.and.organoid.polyA.LumB.expression.correlation.result

df1 <- data.frame(correlation = rs$correlation.matrix[,organoid],
                  type = rep(x='MMC01031',times = nrow(rs$correlation.matrix))
)

df2 <- data.frame(correlation = rs$correlation.matrix[,cell.line],
                  type = rep(x='BT483',times = nrow(rs$correlation.matrix))
)
require(ggplot2)
ggplot(rbind(df1,df2)) + geom_boxplot(aes(x=type,y=correlation))
wilcox.test(df1$correlation,df2$correlation,paired = TRUE)


#Her2
cell.line <- 'EFM192A_BREAST'
organoid  <- 'MMC01031_BREAST_ORGANOID'
rs        <- MET500.vs.CCLE.and.organoid.polyA.Her2.expression.correlation.result

df1 <- data.frame(correlation = rs$correlation.matrix[,organoid],
                  type = rep(x='MMC01031',times = nrow(rs$correlation.matrix))
)

df2 <- data.frame(correlation = rs$correlation.matrix[,cell.line],
                  type = rep(x='EFM192A',times = nrow(rs$correlation.matrix))
)
require(ggplot2)
ggplot(rbind(df1,df2)) + geom_boxplot(aes(x=type,y=correlation))
wilcox.test(df1$correlation,df2$correlation,paired = TRUE)



#Basal
cell.line <- 'HCC70_BREAST'
organoid  <- 'W1009_BREAST_ORGANOID'
rs        <- MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result

df1 <- data.frame(correlation = rs$correlation.matrix[,organoid],
                  type = rep(x='W1009',times = nrow(rs$correlation.matrix))
)

df2 <- data.frame(correlation = rs$correlation.matrix[,cell.line],
                  type = rep(x='HCC70',times = nrow(rs$correlation.matrix))
)
require(ggplot2)
ggplot(rbind(df1,df2)) + geom_boxplot(aes(x=type,y=correlation))
wilcox.test(df1$correlation,df2$correlation,paired = TRUE)



## Pick out good organoid lins
all.organoid.line <- colnames(organoid.log2.rpkm.matrix)
CCLE.non.breast.cancer.cell.line <- setdiff(names(MET500.vs.CCLE.and.organoid.polyA.LumB.expression.correlation.result$cell.line.median.cor),c(CCLE.breast.cancer.cell.line,all.organoid.line))


m                                <- MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% median
s                                <- MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% mad                      # fit the normal distribtuion 
p.value                          <- pnorm(q=MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result$cell.line.median.cor[all.organoid.line],mean = m,sd = s,lower.tail = FALSE) # compute p-value
fdr.vec                          <- p.adjust(p.value,method='fdr')
MET500.Basal.good.organoid       <- names(fdr.vec)[fdr.vec <= 0.01]

m                                <- MET500.vs.CCLE.and.organoid.polyA.LumA.expression.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% median
s                                <- MET500.vs.CCLE.and.organoid.polyA.LumA.expression.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% mad                      # fit the normal distribtuion 
p.value                          <- pnorm(q=MET500.vs.CCLE.and.organoid.polyA.LumA.expression.correlation.result$cell.line.median.cor[all.organoid.line],mean = m,sd = s,lower.tail = FALSE) # compute p-value
fdr.vec                          <- p.adjust(p.value,method='fdr')
MET500.Basal.good.organoid       <- names(fdr.vec)[fdr.vec <= 0.01]






