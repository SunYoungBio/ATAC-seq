
BiocManager::install("DiffBind")

library(DiffBind)
samples <- read.csv(file.path(system.file("extra", package="DiffBind"),
                              "tamoxifen.csv"))

mysamples <- read.table("sampleSheet.txt",header = T)
getwd()
setwd("~/project/wangshuai/ATAC/190429_A00679_0089_AHJWL3DSXX/re1.re2/macs2/diffBind/")
#dbObj <- dba(sampleSheet="sampleSheet.txt")
dbObj <- dba(sampleSheet=mysamples)
mode(dbObj)
dbObj
plot(dbObj)

dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj)
dbObj

#差异分析 Establishing a contrast 
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
dbObj
#  summary of results
dba.show(dbObj, bContrasts=T)
#  overlapping peaks identified by the two different tools (DESeq2 and edgeR)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)

# EdgeR
out <- as.data.frame(comp1.edgeR)
write.table(out, file="all_peaks_edgeR.txt", 
            sep="\t", quote=F, col.names = F, row.names = F)
# DESeq2
out <- as.data.frame(comp1.deseq)
write.table(out, file="all_peaks_deseq2.txt", 
            sep="\t", quote=F, col.names = F, row.names = F)

# Create bed files for each keeping only significant peaks (p < 0.05)
# EdgeR
out <- as.data.frame(comp1.edgeR)
head(out)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "Fold", "FDR")]
write.table(edge.bed, file="Osbon26VSnip26_edgeR_sig.bed", 
            sep="\t", quote=F, row.names=F, col.names=F)

# DESeq2
out <- as.data.frame(comp1.deseq)
deseq.bed <- out[ which(out$FDR < 0.05), 
                  c("seqnames", "start", "end", "Fold", "FDR")]
write.table(deseq.bed, file="Osbon26VSnip26_deseq2_sig.bed", 
            sep="\t", quote=F, row.names=F, col.names=F)

