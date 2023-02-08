##2022-12-12
## This script is for obtaining the tet2 regulating regions
## siTET2 vs siNT, 5hmC data, criteria: ratio<0.8 (DOWN) || ratio>1/0.8(UP)
## overlap tet2 ChIP-seq regions with 5hmC altered regions
## and save the UP and DOWN regions seperately

## then send the tet regulating regions to WSL

## also explore the overlap between tet2 regulating regions and different histone marks

library(openxlsx)
library(dplyr)
library(GenomicRanges)
## tet2 regulating regions (tet2 binding region + 5hmC altered region)
## current criteria: siTET2/siNT ratio < 0.8
tet2_binding_sites_GR <-
  readRDS("E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\tet2_binding_sites_GR.rds")

siTET2_hmC_down_T1 <- 
  read.table("E:\\20220921-WSL-5hmc\\data\\GSE153252_siTET2_5hmC\\hmc_result_table_13_14.xls",
             sep = "\t",header = T,stringsAsFactors = F) %>%
  #filter(read_count.ratio<0.8 & X2.read_count>5) %>%
  filter(read_count.ratio<0.8) %>%
  mutate(start=position,end=position+1) %>%
  dplyr::select(chrom,start,end,strand) %>%
  makeGRangesFromDataFrame()
siTET2_hmC_down_T1 <-
  siTET2_hmC_down_T1[seqnames(siTET2_hmC_down_T1) %in% paste0("chr",c(seq(1,22,by=1),"M"))]
seqlevels(siTET2_hmC_down_T1) <- paste0("chr",c(seq(1,22,by=1),"M"))

siTET2_hmC_down_T2 <- 
  read.table("E:\\20220921-WSL-5hmc\\data\\GSE153252_siTET2_5hmC\\hmc_result_table_15_16.xls",
             sep = "\t",header = T,stringsAsFactors = F) %>%
  #filter(read_count.ratio<0.8 & X4.read_count>5) %>%
  filter(read_count.ratio<0.8 ) %>%
  mutate(start=position,end=position+1) %>%
  dplyr::select(chrom,start,end,strand) %>%
  makeGRangesFromDataFrame()
siTET2_hmC_down_T2 <-
  siTET2_hmC_down_T2[seqnames(siTET2_hmC_down_T2) %in% paste0("chr",c(seq(1,22,by=1),"M"))]
seqlevels(siTET2_hmC_down_T2) <- paste0("chr",c(seq(1,22,by=1),"M"))

length(siTET2_hmC_down_T1) ##134478
length(siTET2_hmC_down_T2) ##267327

## overlap between two experiments
siTET2_5hmC_down <- intersect(siTET2_hmC_down_T1,siTET2_hmC_down_T2)
length(siTET2_5hmC_down) ## 22148

length(findOverlaps(siTET2_5hmC_down,tet2_binding_sites_GR)) ##383

 tet2_5hmC_down_region <- 
  tet2_binding_sites_GR[unique(as.data.frame(findOverlaps(siTET2_5hmC_down,tet2_binding_sites_GR))$subjectHits),] %>%
   makeGRangesFromDataFrame()
 
length(tet2_regulating_region) ##351

## annotate the peaks and save the file
peakAnno.edb <- 
  annotatePeak(tet2_5hmC_down_region, tssRegion = c(-2000, 2000),
               genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
               #ignoreUpstream = TRUE,ignoreDownstream = TRUE,
               TxDb = txdb, annoDb = "org.Hs.eg.db")
write.table(as.data.frame(peakAnno.edb),
            file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\tet2_ChIPSeq_5hmC_DOWN_ol.txt",
            quote = F,
            row.names = F,
            sep = "\t")

## tet2 binding region overlap with histone marks
library(ChIPseeker)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(org.Hs.eg.db)
library(reshape2)
library(dplyr)
library(GenomicRanges)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakFile_ls <- list()
for(chip_file in list.files("E:\\20220921-WSL-5hmc\\data\\ENCODE_MCF7_ChIPseq",full.names = T)){
  BaseName <- str_split(basename(chip_file),pattern = "\\.",simplify = T)[1]
  peakFile_ls[BaseName] <- readPeakFile(chip_file)
  print(paste0("Successfully read file ",BaseName))
}
saveRDS(peakFile_ls,
        file = "E:\\20220921-WSL-5hmc\\analysis\\tet2_enhancer_analysis\\MCF7_histone_ChIPseq_peak.rds")

PeakOverlap <-
  enrichPeakOverlap(queryPeak   = tet2_binding_sites_GR,
                    targetPeak    = peakFile_ls,
                    TxDb          = txdb,
                    pAdjustMethod = "BH",
                    nShuffle      = 50,
                    chainFile     = NULL,
                    verbose       = FALSE)
PeakOverlap
write.table(PeakOverlap,
            file = "E:\\20220921-WSL-5hmc\\analysis\\tet2_enhancer_analysis\\tet2Binding_histone_ol_test.txt",
            quote = F,sep = "\t",row.names = F)


PeakOverlap_2 <-
  enrichPeakOverlap(queryPeak   = tet2_5hmC_down_region,
                    targetPeak    = peakFile_ls,
                    TxDb          = txdb,
                    pAdjustMethod = "BH",
                    nShuffle      = 50,
                    chainFile     = NULL,
                    verbose       = FALSE)
PeakOverlap_2
write.table(PeakOverlap_2,
            file = "E:\\20220921-WSL-5hmc\\analysis\\tet2_enhancer_analysis\\tet2_5hmCDOWN_histone_ol_test.txt",
            quote = F,sep = "\t",row.names = F)

## 5hmC up-regulated regions upon tet2 knock-down
siTET2_hmC_up_T1 <- 
  read.table("E:\\20220921-WSL-5hmc\\data\\GSE153252_siTET2_5hmC\\hmc_result_table_13_14.xls",
             sep = "\t",header = T,stringsAsFactors = F) %>%
  #filter(read_count.ratio<0.8 & X2.read_count>5) %>%
  filter(read_count.ratio > (1/0.8)) %>%
  mutate(start=position,end=position+1) %>%
  dplyr::select(chrom,start,end,strand) %>%
  makeGRangesFromDataFrame()
siTET2_hmC_up_T1 <-
  siTET2_hmC_up_T1[seqnames(siTET2_hmC_up_T1) %in% paste0("chr",c(seq(1,22,by=1),"M"))]
seqlevels(siTET2_hmC_up_T1) <- paste0("chr",c(seq(1,22,by=1),"M"))

siTET2_hmC_up_T2 <- 
  read.table("E:\\20220921-WSL-5hmc\\data\\GSE153252_siTET2_5hmC\\hmc_result_table_15_16.xls",
             sep = "\t",header = T,stringsAsFactors = F) %>%
  #filter(read_count.ratio<0.8 & X4.read_count>5) %>%
  filter(read_count.ratio > (1/0.8)) %>%
  mutate(start=position,end=position+1) %>%
  dplyr::select(chrom,start,end,strand) %>%
  makeGRangesFromDataFrame()
siTET2_hmC_up_T2 <-
  siTET2_hmC_up_T2[seqnames(siTET2_hmC_up_T2) %in% paste0("chr",c(seq(1,22,by=1),"M"))]
seqlevels(siTET2_hmC_up_T2) <- paste0("chr",c(seq(1,22,by=1),"M"))

length(siTET2_hmC_up_T1) ##1819301
length(siTET2_hmC_up_T2) ##1606487

## overlap between two experiments
siTET2_5hmC_up <- intersect(siTET2_hmC_up_T1,siTET2_hmC_up_T2)
length(siTET2_5hmC_up) ## 1313113

length(findOverlaps(siTET2_5hmC_up,tet2_binding_sites_GR)) ##10945

tet2_5hmC_up_region <- 
  tet2_binding_sites_GR[unique(as.data.frame(findOverlaps(siTET2_5hmC_up,tet2_binding_sites_GR))$subjectHits),] %>%
  makeGRangesFromDataFrame()

length(tet2_5hmC_up_region) ##4290
head(tet2_5hmC_up_region)
## annotate the tet2 regulation region and save
peakAnno.edb <- 
  annotatePeak(tet2_5hmC_up_region, tssRegion = c(-2000, 2000),
               genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
               #ignoreUpstream = TRUE,ignoreDownstream = TRUE,
               TxDb = txdb, annoDb = "org.Hs.eg.db")
write.table(as.data.frame(peakAnno.edb),
            file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\tet2_ChIPSeq_5hmC_UP_ol.txt",
            quote = F,
            row.names = F,
            sep = "\t")

## tet2 binding region overlap with histone marks
PeakOverlap_3 <-
  enrichPeakOverlap(queryPeak   = tet2_5hmC_up_region,
                    targetPeak    = peakFile_ls,
                    TxDb          = txdb,
                    pAdjustMethod = "BH",
                    nShuffle      = 50,
                    chainFile     = NULL,
                    verbose       = FALSE)
PeakOverlap_3
write.table(PeakOverlap_3,
            file = "E:\\20220921-WSL-5hmc\\analysis\\tet2_enhancer_analysis\\tet2_5hmCUP_histone_ol_test.txt",
            quote = F,sep = "\t",row.names = F)


## 2022-12-13
## get the enhancer regions, overlap between H3K4me1 and H3K27ac

length(peakFile_ls$H3K27ac)
length(peakFile_ls$H3K4me1)

length(intersect(peakFile_ls$H3K27ac,peakFile_ls$H3K4me1))

# Prepare the promotor regions
promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
# Calculate the tag matrix
tagMatrixList <- getTagMatrix(intersect(peakFile_ls$H3K27ac,peakFile_ls$H3K4me1), 
                              windows=promoter)

## Profile plots
plotAvgProf(tagMatrixList, xlim=c(-2000, 2000), conf=0.95,resample=500)

