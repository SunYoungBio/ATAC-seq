#source ("https://bioconductor.org/biocLite.R")
#BiocManager::install("ChIPseeker")
#载入各种包
library("ChIPseeker")
library(clusterProfiler)
#library("org.Hs.eg.db")
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#library("TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts")
#lincRNA_txdb=TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts



library(GenomicFeatures)
#txdb=makeTxDbFromGFF("genomes/tigr7/gff3_files/tigr7.gff3")
#saveDb(txdb,"TxDb.tigr7")
txdb=loadDb("~/R/TxDb.tigr7")

seqlevels(txdb)
cols=columns(txdb)
keytypes(txdb)
keys=tail(keys(txdb))
select(txdb, keys = keys, columns="TXNAME", keytype="GENEID")
aa=select(txdb, keys=keys, columns=cols, keytype="GENEID")

nip26.uniq=readPeakFile("project/wangshuai/ATAC/190429_A00679_0089_AHJWL3DSXX/re1.re2/macs2/nip26.uniq.bed")
bon22.uniq=readPeakFile("project/wangshuai/ATAC/190429_A00679_0089_AHJWL3DSXX/re1.re2/macs2/bon22.uniq.bed")
#covplot(nip26)
setwd("project/wangshuai/ATAC/osbon.nip.35/")
setwd("project/wangshuai/ATAC/190429_A00679_0089_AHJWL3DSXX/re1.re2/macs2/diffBind/")
peak1=readPeakFile("Osbon26downPeak.bed")
peak2=readPeakFile("Osbon26upPeak.bed")

peaksList=list(Nip26=peak1,Bon26=peak2)
peaksList=list(Osbon26downPeak=peak1,Osbon26upPeak=peak2)

promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=0)
#tagMatrix.nip26 <- getTagMatrix(nip26, windows=promoter)

#tagMatrixList <- lapply(peaks.nip.bon, getTagMatrix, windows=promoter)
#annotatePeak传入annoDb参数,可进行基因ID转换（Entrez，ENSEMBL，SYMBOL，GENENAME）
options(ChIPseeker.downstreamDistance = 1000)#设置下游区间，最好与promoter区间相同。
peakAnno <- annotatePeak(peak1, 
                         tssRegion=c(-1000, 0),
                         overlap = "all",#注意：all输出离geneBody最近的基因，
                         #否则有可能默认离TSS最近，导致一个基因与远端peak overlap，
                         #注释到是下一个基因。
                         TxDb=txdb)
peakAnno
plotAnnoPie(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
#peakAnno@detailGenomicAnnotation

upsetplot(peakAnno, vennpie=TRUE)
plotAnnoBar(peakAnno)
mode(peakAnno)
str(peakAnno)
bon=as.data.frame(peakAnno)
write.table(bon,"Osbon26VSnip26_edgeR_sig_peakAnno.txt",sep='\t',quote = F)

plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

peakAnnoList <- lapply(files, annotatePeak, 
                       TxDb=txdb,tssRegion=c(-3000, 3000))
plotAnnoBar(peakAnnoList)

plotDistToTSS(peakAnnoList)

genes <- lapply(peakAnnoList, function(i) 
  as.data.frame(i)$geneId)
vennplot(genes[2:4], by='Vennerable')


getwd()

peakAnnoList <- lapply(peaksList, annotatePeak, TxDb=txdb,
                       tssRegion=c(-1000, 0), verbose=FALSE,
                       addFlankGeneInfo=TRUE, 
                       overlap = "all",
                       flankDistance=3000#,
                       #annoDb="org.Hs.eg.db"
                       )

peakAnnoList
plotAnnoBar(peakAnnoList)
#plotAnnoPie(peakAnnoList)
plotDistToTSS(peakAnnoList,title="Distribution of transcription factor-binding loci \n relative to TSS")
genes <- lapply(peakAnnoList, function(i) 
  as.data.frame(i)$geneId)
vennplot(genes[2:4], by='Vennerable')

# Create a list with genes from each sample
gene = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
# Run GO enrichment analysis 
ego <- enrichGO(gene = entrez, 
                keytype = "ENTREZID", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
# Dotplot visualization
dotplot(ego, showCategory=50)
# Multiple samples KEGG analysis
compKEGG <- compareCluster(geneCluster = gene, 
                           fun = "enrichKEGG",
                           organism = "human",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")

# Output peakAnnolist file
save(peakAnnoList,file="peakAnnolist.rda")
write.table(as.data.frame(peakAnnoList$Osbon26downPeak),file="Osbon26down.PeakAnno",
            sep='\t',
            quote = F, row.names = F)
write.table(as.data.frame(peakAnnoList$Osbon26upPeak),file="Osbon26up.PeakAnno",
            sep='\t',quote = F,row.names = F)

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "results/clusterProfiler_Nanog.csv")

