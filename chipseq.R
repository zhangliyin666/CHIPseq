?BiocManager
BiocManager::install("ChIPpeakAnno")
library(ChIPpeakAnno)
?ChIPpeakAnno
#1.找overlap
#读取文件所的所在路径
bed <- system.file("extdata","MACS_output.bed", package="ChIPpeakAnno")
#将bed文件转换为 Ranges格式
gr1 <- toGRanges(bed,format="BED", header=FALSE)
#读取文件所的所在路径
gff<-system.file("extdata","GFF_peaks.gff",package="ChIPpeakAnno")
#将gff文件转换为GRanges
gr2<-toGRanges(gff,format="GFF",header=FALSE,skip=3)
gr1
#library(GenomicFeatures)
#txdb <- makeTxDbFromGFF('H:/R-4.0.3/library/ChIPpeakAnno/extdata/GFF_peaks.gff')
#anno <- toGRanges(txdb, format='gene')
#anno[1:2]
#查看格式，找重合需要格式为numeric
class(gr1$score)
class(gr2$score)
ol <- findOverlapsOfPeaks(gr1, gr2)
#找gr1和gr2都有的peak
length(ol$peaklist$`gr1///gr2`)
#仅gr1有的peak
length(ol$peaklist$gr1)
#仅gr2有的peak
length(ol$peaklist$gr2)
#查看重叠peak的信息（前两行）
ol$peaklist[["gr1///gr2"]][1:2]
# 或 ol$peaklist$`gr1///gr2`
#韦恩图数据
ol$venn_cnt
#作韦恩图
makeVennDiagram(ol, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2")) # label color, keep same as circle border colors
#2.准备注释数据（跳转5）
BiocManager::install("EnsDb.Hsapiens.v75")
BiocManager::install("EnsDb.Mmusculus.v79")
library(EnsDb.Hsapiens.v75) ##(hg19)
annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")
#3.结合位点可视化
#对基因组feature（例如转录起始位点TSS）附近的重叠peaks进行可视化
overlaps <- ol$peaklist[["gr1///gr2"]]
binOverFeature(overlaps, annotationData=annoData,
               radius=5000, nbins=20, FUN=length, errFun=0,
               ylab="count", 
               main="Distribution of aggregated peak numbers around TSS")
#对peaks的分布区域进行概括
#可能一个peak会跨越多种类型的基因feature信息，因此可能看到输出结果中注释到的feature数量比输入的peaks数量更多
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
aCR<-assignChromosomeRegion(overlaps, nucleotideLevel=FALSE, 
                            precedence=c("Promoters", "immediateDownstream", 
                                         "fiveUTRs", "threeUTRs", 
                                         "Exons", "Introns"), 
                            TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
barplot(aCR$percentage, las=3)
#4.提取peak周围的序列
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
seq <- getAllPeakSequence(gr1, upstream=20, downstream=20, genome=Hsapiens)
#5.进行motif分析
#用1号染色体的碱基分布当做背景
freqs <- oligoFrequency(Hsapiens$chr1, MarkovOrder=3)
#oligoLength规定了motif的长度
os <- oligoSummary(seq, oligoLength=6, MarkovOrder=3,
                   quickMotif=TRUE, freqs=freqs)
zscore <- sort(os$zscore)
#绘制所有6个碱基组合的频率分布图
h <- hist(zscore, breaks=100, xlim=c(-50, 50), main="Histogram of Z-score")
#频率最大的碱基组合即为motif的结果
text(zscore[length(zscore)], max(h$counts)/10,
     labels=names(zscore[length(zscore)]), adj=1)
#绘制sequence logo
BiocManager::install("motifStack")
library(motifStack)
pfms <- mapply(function(.ele, id)
  new("pfm", mat=.ele, name=paste("SAMPLE motif", id)),
  os$motifs, 1:length(os$motifs))
motifStack(pfms[[1]])
#5.进行基因注释
overlaps.anno <- annotatePeakInBatch(gr1,
                                     AnnotationData=annoData,
                                     output="nearestLocation"
)
library(org.Hs.eg.db)
overlaps.anno <- addGeneIDs(overlaps.anno,
                            "org.Hs.eg.db",
                            IDs2Add = "entrez_id")
pie1(table(overlaps.anno$insideFeature),radius = 3,cex = 0.8)

over <- getEnrichedGO(overlaps.anno, orgAnn="org.Hs.eg.db",
                      maxP=.05, minGOterm=10,
                      multiAdjMethod="BH", condense=TRUE)
