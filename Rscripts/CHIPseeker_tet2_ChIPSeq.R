## Exploration on tet2 ChIP-seq data
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(org.Hs.eg.db)
library(reshape2)
library(dplyr)
library(GenomicRanges)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

##Load the data
## tet2 ChIP-seq data
tet2_binding_sites <-
  read.table("E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\hglft_genome_32e07_39b920.bed",
           sep = "\t", header = F)
colnames(tet2_binding_sites) <- c("seqnames","start","end","name","addin")
tet2_binding_sites_GR <- makeGRangesFromDataFrame(tet2_binding_sites)
saveRDS(tet2_binding_sites_GR,
        file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\tet2_binding_sites_GR.rds")
tet2_binding_window<- 
  makeBioRegionFromGranges(gr=tet2_binding_sites_GR,by="peak",type="body",
                           upstream = rel(0.2),
                           downstream = rel(0.2))

## annotate the peaks with ChIPseeker
tmp <- readPeakFile("E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\hglft_genome_32e07_39b920.bed",
                    header=F)
peakAnno.edb <- 
  annotatePeak(tmp, tssRegion = c(-2000, 2000),
               genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
               #ignoreUpstream = TRUE,ignoreDownstream = TRUE,
               TxDb = txdb, annoDb = "org.Hs.eg.db")
write.table(as.data.frame(peakAnno.edb),
            file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\tet2_ChIPSeq_annotated.txt",
            quote = F,
            row.names = F,
            sep = "\t")

## plot the average plot for differential peaks
out_path <- "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\diffPeak_tet2_avg_plots"
dir.create(out_path)
for(diff_peak in list.files("E:\\20220921-WSL-5hmc\\data\\diff_peak_new",pattern = "*deseq2.txt$",full.names = TRUE)){
  BaseName <- paste(str_split(basename(diff_peak),"_",simplify = T)[1:3],collapse = "-")
  
  tmp_up_df <- read.table(diff_peak,header = T,stringsAsFactors = F) %>%
    filter(Fold>1 & p.value<0.05) %>%
    makeGRangesFromDataFrame()
  tmp_down_df <- read.table(diff_peak,header = T,stringsAsFactors = F) %>%
    filter(Fold<(-1) & p.value<0.05) %>%
    makeGRangesFromDataFrame()
  tagMatrix_up <- getTagMatrix(tmp_up_df,
                               windows = tet2_binding_window,
                               upstream = rel(0.2),
                               downstream = rel(0.2),
                               nbin = 100)
  tagMatrix_down <- getTagMatrix(tmp_down_df, 
                                 windows = tet2_binding_window,
                                 nbin = 100)
  
  if(!is.null(ncol(tagMatrix_up)) & dim(tagMatrix_up)[1]>1){
    out_file_name <- paste(c(BaseName,"_UP_","tagAvgProf"),collapse = "")
    pdf(paste0(out_path,"\\",out_file_name,".pdf"),width = 10,height = 8)
#    plot(plotAvgProf(tagMatrix_up, xlim=c(-100, 100),
#                     xlab="tet2 binding region (5'->3')",
#                     ylab = "Read Count Frequency"))
   
    print(plotPeakProf(tagMatrix_up,conf = 0.95))
    dev.off()
  }
  if(!is.null(ncol(tagMatrix_down)) & dim(tagMatrix_down)[1]>1){
  ## plot 
  out_file_name <- paste(c(BaseName,"_DOWN_","tagAvgProf"),collapse = "")
  pdf(paste0(out_path,"\\",out_file_name,".pdf"),width = 10,height = 8)
#  plot(plotAvgProf(tagMatrix_down, xlim=c(-100, 99),
#                  xlab="tet2 binding region (5'->3')", ylab = "Read Count Frequency"))
  print(plotPeakProf(tagMatrix_down,conf = 0.95))
  dev.off()
  }
}

## overlap of tet2 ChIP-seq with DMRs and DHMRs
target_ls <- list()
for(diff_peak in list.files("E:\\20220921-WSL-5hmc\\data\\diff_peak_new",pattern = "*deseq2.txt$",full.names = TRUE)){
  BaseName <- paste(str_split(basename(diff_peak),"_",simplify = T)[1:3],collapse = "-")
  
  tmp_up_df <- read.table(diff_peak,header = T,stringsAsFactors = F) %>%
    filter(Fold>1 & p.value<0.05) %>%
    makeGRangesFromDataFrame()
  target_ls[paste0(BaseName,"_up")] <- tmp_up_df
  tmp_down_df <- read.table(diff_peak,header = T,stringsAsFactors = F) %>%
    filter(Fold<(-1) & p.value<0.05) %>%
    makeGRangesFromDataFrame()
  target_ls[paste0(BaseName,"_down")] <- tmp_down_df
}
PeakOverlap <-
  enrichPeakOverlap(queryPeak   = tet2_binding_sites_GR,
                  targetPeak    = target_ls,
                  TxDb          = txdb,
                  pAdjustMethod = "BH",
                  nShuffle      = 50,
                  chainFile     = NULL,
                  verbose       = FALSE)
write.table(as.data.frame(PeakOverlap),
            file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\diff_peak_tet2binding_OL.txt",
            sep = "\t",quote = F,row.names = F)

## tet2 binding and non-binding region foldchange (need to intersect peaks in server)

## annotate the diff_peak files with overlap status with tet2
library(stringr)
out_path <- "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq"
for(diff_peak in list.files("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno",
                            pattern = "*peak_deseq2.txt_anno$",full.names = TRUE)){
  BaseName <- paste(c(str_split(basename(diff_peak),"_",simplify = T)[1:3],"tet2_ol.txt"),collapse = "_")
  
  tmp_peak <- readPeakFile(diff_peak)
  overlap_index <- as.data.frame(GenomicRanges::findOverlaps(tet2_binding_sites_GR,tmp_peak))$subjectHits
  tmp_peak_df <- as.data.frame(tmp_peak) %>% 
    mutate(overlap_with_tet2="no")
  tmp_peak_df[overlap_index,"overlap_with_tet2"] <- "yes"
  write.table(tmp_peak_df,
              file = paste0(out_path,"\\",BaseName),
              quote = F, row.names = F,sep = "\t")
  }

## plot the fold changes of peaks that overlap or don't overlap with tet2 binding region
library(ggplot2)
library(gridExtra)
for(diff_peak in list.files("E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq",
                            pattern = "*tet2_ol.txt$",full.names = TRUE)){
  BaseName <- paste(str_split(basename(diff_peak),"_",simplify = T)[1:3],collapse = "_")
  annotated_diff_peak <- 
    read.table(diff_peak,header = T,
                                    sep = "\t",stringsAsFactors = F,
                                    quote = "\"")
  annotated_diff_peak$fold_group <- 
    ifelse(annotated_diff_peak$Fold>0,"Fold_gt0","Fold_lt0")
  
  p1 <- 
    annotated_diff_peak %>%
    filter(fold_group=="Fold_gt0") %>%
    ggplot(aes(x=overlap_with_tet2,y=Fold,fill=overlap_with_tet2))+
    geom_boxplot(width=0.6)+
    scale_y_continuous(limits = c(0,quantile(annotated_diff_peak[which(annotated_diff_peak$fold_group=="Fold_gt0"),"Fold"],0.9)))+
    scale_x_discrete(name="",breaks=c("no","yes"),labels=c("non-OL","OL"))+
    scale_fill_manual(name="overlap with \ntet2 ChIP-seq",
                      values = list("no"="#6F7EB3","yes"="#D05B5B"),
                      labels=c("NO","YES"))+
    labs(title="Fold > 0 ")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10,colour = "black"),
          axis.title = element_text(size = 12,colour = "black"),
          plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 13),
          legend.position = "NA")
  p2 <- 
    annotated_diff_peak %>%
    filter(fold_group=="Fold_lt0") %>%
    ggplot(aes(x=overlap_with_tet2,y=Fold,fill=overlap_with_tet2))+
    geom_boxplot(width=0.6)+
    scale_y_continuous(limits = c(quantile(annotated_diff_peak[which(annotated_diff_peak$fold_group=="Fold_lt0"),"Fold"],0.1),0))+
    scale_x_discrete(name="",breaks=c("no","yes"),labels=c("non-OL","OL"))+
    scale_fill_manual(name="overlap with \ntet2 ChIP-seq",
                      values = list("no"="#6F7EB3","yes"="#D05B5B"),
                      labels=c("NO","YES"))+
    labs(title="Fold < 0 ")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10,colour = "black"),
          axis.title = element_text(size = 12,colour = "black"),
          plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 13))
  
  p <- 
    grid.arrange(p1,p2,ncol=2,widths=c(1,1.5),top=BaseName)
  ggsave(p,filename = paste0("E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\tet2_binding_region_Fold_boxplot\\",BaseName,".pdf"),
         width = 8,height = 5)
  
  
  }

## annotate the overlap status with tet2 ChIP-seq data of each sample
library(dplyr)
library(ChIPseeker)
## initiate the final dataframe
peak_annotated <- data.frame()
out_path <- "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq"
for(diff_peak in list.files("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno",
                            pattern = "^[ADJUP].*H_peaks.narrowPeak_anno$",full.names = TRUE)){
  BaseName <- str_split(basename(diff_peak),"_",simplify = T)[1]
  
  tmp_peak <- readPeakFile(diff_peak)
  overlap_index <- as.data.frame(GenomicRanges::findOverlaps(tet2_binding_sites_GR,tmp_peak))$subjectHits
  tmp_peak_df <- as.data.frame(tmp_peak) %>% 
    mutate(overlap_with_tet2="no") %>%
    mutate(sample=!!BaseName) %>%
    dplyr::select(V7,annotation,ENSEMBL,SYMBOL,overlap_with_tet2,sample)
  tmp_peak_df[overlap_index,"overlap_with_tet2"] <- "yes"
  colnames(tmp_peak_df)[1] <- "fold"
  peak_annotated <- rbind(peak_annotated,tmp_peak_df)
  #write.table(tmp_peak_df,
  #            file = paste0(out_path,"\\",BaseName),
  #            quote = F, row.names = F,sep = "\t")
}

write.table(peak_annotated,
            file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\Sample_peak_tet2_bind_anno.txt",
            quote = F, sep = "\t", row.names = F)
## save the data in RDS format
saveRDS(peak_annotated,file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\Sample_peak_tet2_bind_anno.rds")

##plot the fold of peaks in each sample
peak_annotated <-
peak_annotated %>% mutate(stage=case_when(sample %in% c("U3H","U5H","JH") ~ "UDH",
                                          sample %in% c("A10H","A2H","AWH") ~ "ADH",
                                          sample %in% c("D10H","D20H","D3H") ~ "DCIS",
                                          sample %in% c("PT3H","PT5H","PT6H","PT7H") ~ "PT")) %>%
  mutate(stage=factor(stage,levels = c("UDH","ADH","DCIS","PT")))
library(ggplot2)
p <-
ggplot(peak_annotated,aes(x=stage,y=fold,fill=overlap_with_tet2))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Dark2")+
  labs(x="",y="5hmC fold enrichment")+
  scale_y_continuous(limits = c(0,10))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12,color = "black"),
        axis.title = element_text(size=14,color = "black"))
ggsave(p,
       filename = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\tet2_ol_peak_fold_boxplot.pdf",
       width = 10,height = 8)


## test before loop
test <- "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\A2H_peaks.narrowPeak_anno"
tmp_peak <- readPeakFile(test)
overlap_index <- as.data.frame(GenomicRanges::findOverlaps(tet2_binding_sites_GR,tmp_peak))$subjectHits


