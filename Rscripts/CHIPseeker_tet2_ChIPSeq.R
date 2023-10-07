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
# tet2_binding_sites <-
#   read.table("E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\hglft_genome_32e07_39b920.bed",
#            sep = "\t", header = F)
# colnames(tet2_binding_sites) <- c("seqnames","start","end","name","addin")
# tet2_binding_sites_GR <- makeGRangesFromDataFrame(tet2_binding_sites)
# saveRDS(tet2_binding_sites_GR,
#         file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\tet2_binding_sites_GR.rds")
tet2_binding_sites_GR <- readRDS("E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\tet2_binding_sites_GR.rds")
#tet2_binding_window <- 
#  makeBioRegionFromGranges(gr=tet2_binding_sites_GR,by="peak",type="body",
#                           upstream = rel(0.2),
#                           downstream = rel(0.2))

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
##step1: make Bioregions from tet2-binding peaks
## +-3kb around peak center
tet2_binding_for_covplot <- 
  GRanges(seqname=tet2_binding_sites_GR@seqnames,
          strand=tet2_binding_sites_GR@strand,
          ranges=IRanges(start=tet2_binding_sites_GR@ranges@start+as.integer(tet2_binding_sites_GR@ranges@width/2)-3000,
                         width = 6000)) %>%
  makeBioRegionFromGranges(type = "body",by="gene")

out_path <- "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\diffPeak_tet2_avg_plots"
dir.create(out_path)
for(diff_peak in list.files("E:\\20220921-WSL-5hmc\\data\\diff_peak_new",
                            pattern = "*deseq2.txt$",full.names = TRUE)){
  BaseName <- paste(str_split(basename(diff_peak),"_",simplify = T)[1:3],
                    collapse = "-")
  
  tmp_up_df <- read.table(diff_peak,header = T,stringsAsFactors = F) %>%
    filter(Fold>1 & p.value<0.05) %>%
    makeGRangesFromDataFrame()
  tmp_down_df <- read.table(diff_peak,header = T,stringsAsFactors = F) %>%
    filter(Fold<(-1) & p.value<0.05) %>%
    makeGRangesFromDataFrame()
  tagMatrix_up <- getTagMatrix(tmp_up_df,
                               windows = tet2_binding_for_covplot,
                               nbin = 300
                               )
  tagMatrix_down <- getTagMatrix(tmp_down_df, 
                                 windows = tet2_binding_for_covplot,
                                 nbin = 300)
  
  if(!is.null(ncol(tagMatrix_up)) & dim(tagMatrix_up)[1]>1){
    p <- plotPeakProf(tagMatrix_up,conf = 0.95)+
      labs(title = paste0(BaseName,"_Fold1"))+
      scale_x_continuous(breaks = c(0,75,150,225,300),
                         labels = c("-3kb","-1.5kb","center","+1.5kb","+3kb"))+
      theme(axis.text = element_text(size = 10,color = "black"),
            plot.title = element_text(hjust = 0.5,vjust = 0.5))
    ggsave(filename = paste0(BaseName,"_Fold1.pdf"),
           p,width = 5,height = 4)
  }
  if(!is.null(ncol(tagMatrix_down)) & dim(tagMatrix_down)[1]>1){
  ## plot 
  p <- plotPeakProf(tagMatrix_down,conf = 0.95)+
    labs(title = paste0(BaseName,"_FoldM1"))+
    scale_x_continuous(breaks = c(0,75,150,225,300),
                       labels = c("-3kb","-1.5kb","center","+1.5kb","+3kb"))+
    theme(axis.text = element_text(size = 10,color = "black"),
          plot.title = element_text(hjust = 0.5,vjust = 0.5))
  ggsave(filename = paste0(BaseName,"_FoldM1.pdf"),
         p,width = 5,height = 4)
  }
}

## test on overlap of tet2 ChIP-seq with DMRs and DHMRs
target_ls <- list()
for(diff_peak in list.files("E:\\20220921-WSL-5hmc\\data\\diff_peak_new",
                            pattern = "*deseq2.txt$",full.names = TRUE)){
  BaseName <- paste(str_split(basename(diff_peak),"_",simplify = T)[1:3],
                    collapse = "-")
  
  tmp_up_df <- read.table(diff_peak,header = T,stringsAsFactors = F) %>%
    filter(Fold>1 & p.value<0.05) %>%
    makeGRangesFromDataFrame()
  target_ls[paste0(BaseName,"_fold1")] <- tmp_up_df
  tmp_down_df <- read.table(diff_peak,header = T,stringsAsFactors = F) %>%
    filter(Fold<(-1) & p.value<0.05) %>%
    makeGRangesFromDataFrame()
  target_ls[paste0(BaseName,"_foldM1")] <- tmp_down_df
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
  BaseName <- paste(c(str_split(basename(diff_peak),"_",simplify = T)[1:3],"tet2_ol.txt"),
                    collapse = "_")
  
  tmp_peak <- readPeakFile(diff_peak)
  overlap_index <- as.data.frame(GenomicRanges::findOverlaps(tet2_binding_sites_GR,
                                                             tmp_peak))$subjectHits
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
saveRDS(peak_annotated,
        file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\Sample_peak_tet2_bind_anno.rds")

peak_annotated <-
        readRDS(file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\Sample_peak_tet2_bind_anno.rds")

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
  geom_boxplot(outlier.colour = "NA")+
  scale_fill_brewer(palette = "Dark2")+
  labs(x="",y="5hmC fold enrichment")+
  scale_y_continuous(limits = c(1,8))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12,color = "black"),
        axis.title = element_text(size=14,color = "black"))
ggsave(p,
       filename = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\tet2_ol_peak_fold_boxplot.pdf",
       width = 5.5,height = 4)


## plot the number of overlapped peaks between DhMRs and tet2-binding region
library(stringr)
library(dplyr)
data_dir <- "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq"
DhMRs_with_tet2 <- data.frame()
DhMRs_with_tet2_region <- data.frame()
for(i in list.files(data_dir,pattern = "*_H_tet2_ol.txt$",full.names = T)){
  comp_name <- 
    paste(str_split(basename(i),"_",simplify = T)[1:2],collapse = "_")
  DhMRs_file <-
    read.table(i,
               header = T,
               stringsAsFactors = F,
               sep = "\t",
               quote = "\"")
  Fold1_num_tet2 <-
    DhMRs_file %>%
      filter(Fold>1 & p.value<0.05 & overlap_with_tet2=="yes") %>%
      nrow()
  Fold1_num <-
    DhMRs_file %>%
    filter(Fold>1 & p.value<0.05)
  FoldM1_num_tet2 <-
    DhMRs_file %>%
    filter(Fold<(-1) & p.value<0.05 & overlap_with_tet2=="yes") %>%
    nrow()
  FoldM1_num <-
    DhMRs_file %>%
    filter(Fold<(-1) & p.value<0.05)
  tmp_df <- 
    data.frame(sample=comp_name,
               change=c("Fold1","FoldM1"),
               number=c(Fold1_num_tet2,FoldM1_num_tet2),
               DhMR_num=c(nrow(Fold1_num),nrow(FoldM1_num)),
               prop=c(Fold1_num_tet2/nrow(Fold1_num),FoldM1_num_tet2/nrow(FoldM1_num))
               )
  DhMRs_with_tet2 <- rbind(DhMRs_with_tet2,tmp_df)
  
  Fold1_freq <- Fold1_num %>%
    mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                        ifelse(grepl("^Exon",annotation),"Exon",annotation))) %>%
    group_by(group) %>%
    summarise(count=n()) %>%
    mutate(sample=comp_name, change="Fold1")
  
  FoldM1_freq <- FoldM1_num %>%
    mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                        ifelse(grepl("^Exon",annotation),"Exon",annotation))) %>%
    group_by(group) %>%
    summarise(count=n()) %>%
    mutate(sample=comp_name, change="FoldM1")
  
  tmp_freq_df <- rbind(Fold1_freq,FoldM1_freq)
  DhMRs_with_tet2_region <- rbind(DhMRs_with_tet2_region,tmp_freq_df)
}
write.table(DhMRs_with_tet2,
            file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\DhMRs_with_tet2_num.txt",
            row.names = F,
            quote = F,
            sep = "\t"
            )
write.table(DhMRs_with_tet2_region,
            file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\DhMRs_with_tet2_region_num.txt",
            row.names = F,
            quote = F,
            sep = "\t"
)

library(ggplot2)
library(ggsci)
p<- 
  ggplot(DhMRs_with_tet2,
         aes(x=factor(sample,levels = c("UDH_ADH","ADH_DCIS","DCIS_PT")),
             y=number,fill=change))+
    geom_bar(stat = "identity",position = "dodge",width = 0.6)+
    geom_text(aes(label=number),vjust=0,
              position = position_dodge(width=0.6))+
    scale_y_continuous(limits = c(0,950),expand = c(0,0))+
    labs(x="",y="Frequency",title = "DhMRs overlapped with tet2-BR")+
    scale_fill_manual(values = c("Fold1"="#0081a7","FoldM1"="#f28482"))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10,colour = "black"),
          plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 12))
ggsave(p,filename = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\DhMRs_tet2_num_barplot.pdf",
       width = 6,height = 4)

## plot for genomic region distribution
DhMRs_with_tet2_region$sample <- 
  factor(DhMRs_with_tet2_region$sample,
         levels = c("UDH_ADH","ADH_DCIS","DCIS_PT"))
p <- 
  ggplot(DhMRs_with_tet2_region,
         aes(y=change,
             x=count,fill=group))+
    geom_col(position = "fill")+
    facet_grid(sample~.)+
    labs(x="proportion",y="")+
    scale_x_continuous(expand = expansion(mult = c(0,0.01)))+
    scale_fill_nejm()+
    guides(fill=guide_legend(reverse=TRUE))+
    theme_classic()+
    theme(strip.background = element_rect(colour = NA),
          axis.text = element_text(size = 10,color="black"))
ggsave(p,
       filename = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\DhMRs_tet2_region_num_barplot.pdf",
       width = 8,height = 4)


## same analysis on DMRs
DMRs_with_tet2 <- data.frame()
DMRs_with_tet2_region <- data.frame()
for(i in list.files(data_dir,pattern = "*_M_tet2_ol.txt$",full.names = T)){
  comp_name <- 
    paste(str_split(basename(i),"_",simplify = T)[1:2],collapse = "_")
  DMRs_file <-
    read.table(i,
               header = T,
               stringsAsFactors = F,
               sep = "\t",
               quote = "\"")
  Fold1_num_tet2 <-
    DMRs_file %>%
    filter(Fold>1 & p.value<0.05 & overlap_with_tet2=="yes") %>%
    nrow()
  Fold1_num <-
    DMRs_file %>%
    filter(Fold>1 & p.value<0.05)
  FoldM1_num_tet2 <-
    DMRs_file %>%
    filter(Fold<(-1) & p.value<0.05 & overlap_with_tet2=="yes") %>%
    nrow()
  FoldM1_num <-
    DMRs_file %>%
    filter(Fold<(-1) & p.value<0.05)
  tmp_df <- 
    data.frame(sample=comp_name,
               change=c("Fold1","FoldM1"),
               number=c(Fold1_num_tet2,FoldM1_num_tet2),
               DMR_num=c(nrow(Fold1_num),nrow(FoldM1_num)),
               prop=c(Fold1_num_tet2/nrow(Fold1_num),FoldM1_num_tet2/nrow(FoldM1_num))
    )
  DMRs_with_tet2 <- rbind(DMRs_with_tet2,tmp_df)
  
  Fold1_freq <- Fold1_num %>%
    mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                        ifelse(grepl("^Exon",annotation),"Exon",annotation))) %>%
    group_by(group) %>%
    summarise(count=n()) %>%
    mutate(sample=comp_name, change="Fold1")
  
  FoldM1_freq <- FoldM1_num %>%
    mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                        ifelse(grepl("^Exon",annotation),"Exon",annotation))) %>%
    group_by(group) %>%
    summarise(count=n()) %>%
    mutate(sample=comp_name, change="FoldM1")
  
  tmp_freq_df <- rbind(Fold1_freq,FoldM1_freq)
  DMRs_with_tet2_region <- rbind(DMRs_with_tet2_region,tmp_freq_df)
}

write.table(DMRs_with_tet2,
            file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\DMRs_with_tet2_num.txt",
            row.names = F,
            quote = F,
            sep = "\t"
)
write.table(DMRs_with_tet2_region,
            file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\DMRs_with_tet2_region_num.txt",
            row.names = F,
            quote = F,
            sep = "\t"
)

p<- 
  ggplot(DMRs_with_tet2,
         aes(x=factor(sample,levels = c("UDH_ADH","ADH_DCIS")),
             y=number,fill=change))+
  geom_bar(stat = "identity",position = "dodge",width = 0.6)+
  geom_text(aes(label=number),vjust=0,
            position = position_dodge(width=0.6))+
  scale_y_continuous(limits = c(0,125),expand = c(0,0))+
  labs(x="",y="Frequency",title = "DMRs overlapped with tet2-BR")+
  scale_fill_manual(values = c("Fold1"="#0081a7","FoldM1"="#f28482"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 12))
ggsave(p,filename = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\DMRs_tet2_num_barplot.pdf",
       width = 6,height = 4)

## plot for genomic region distribution
DMRs_with_tet2_region$sample <- 
  factor(DMRs_with_tet2_region$sample,
         levels = c("UDH_ADH","ADH_DCIS"))
p <- 
  ggplot(DMRs_with_tet2_region,
         aes(y=change,
             x=count,fill=group))+
  geom_col(position = "fill",width = 0.7)+
  #coord_fixed(ratio = 1)+
  facet_grid(sample~.)+
  labs(x="proportion",y="")+
  scale_x_continuous(expand = expansion(mult = c(0,0.01)))+
  scale_fill_nejm()+
  guides(fill=guide_legend(reverse=TRUE))+
  theme_classic()+
  theme(strip.background = element_rect(colour = NA),
        axis.text = element_text(size = 10,color="black"),
        aspect.ratio = 0.3)
ggsave(p,
       filename = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\DMRs_tet2_region_num_barplot.pdf",
       width = 8,height = 3)


## plot the expression of tet2 in each 
RNAseq_df <- 
read.table("E:\\20220921-WSL-5hmc\\pre-info\\13059_2014_3265_MOESM4_ESM.txt",
           skip = 2,sep = "\t",comment.char = "",
           row.names = 1)
meta_infor <- read.table("E:\\20220921-WSL-5hmc\\pre-info\\13059_2014_3265_MOESM4_ESM.txt",
                         nrows = 2,
                         sep = "\t",comment.char = "")

tet2_expr <- 
  data.frame(value=unlist(RNAseq_df["TET2",]),
             group=unlist(meta_infor[2,-1])) %>%
  filter(!group=="nl")
tet2_expr$group <- factor(tet2_expr$group,
                          levels = c("normal","en","dcis","idc"))
ggplot(tet2_expr,aes(group,value,fill=group))+
  geom_violin()+
  geom_point()+
  guides()


## compare the expression levels of tet2-targeted RNAs and others

tet2_br <- 
  read.table("E:/20220921-WSL-5hmc/analysis/GSE153251_tet2_ChIPSeq/tet2_ChIPSeq_annotated.txt",
            header = T, sep = "\t", quote = "\"",
            stringsAsFactors = F,comment.char = "")
dim(tet2_br)
tet2_targets <- unique(tet2_br$ENSEMBL)
## 7766

## load gene expression data
data_dir <- "E:\\20220921-WSL-5hmc\\data\\RNAseq_data"
normalized_count_df <- 
  readRDS(paste0(data_dir,"//normalized_count_df.rds"))

plot_RNA_group_Tagcov <- 
  function(RNA_expr_nor, sample_name){
    ## RNA_expr_nor need to be a dataframe with:
    ## gene_id as rownames, expression value as column (named value)
    
    # get gene ids in each quantile group
    RNA_expr_nor <- 
      RNA_expr_nor %>%
      filter(!value==0) %>%
      mutate(group=ifelse(row.names(.) %in% tet2_targets,"target","not"))
    p <-
      ggplot(RNA_expr_nor,aes(log2(value),fill=group))+
        #geom_boxplot()
        geom_density(alpha=0.6)+
        scale_x_continuous(expand = c(0,0),name = "log2(normalized value)")+
        scale_fill_manual(values = c("not"="#f18d91","target"="#639ccf"))+
        scale_y_continuous(expand = c(0,0),limits = c(0,0.2))+
        theme_classic()+
        theme(panel.grid = element_blank())
    
    ggsave(p, 
           filename=paste0("E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\",
                           sample_name,"_tet2vsNo.pdf"),
           width=13, height=8) 
  }
for(i in colnames(normalized_count_df)){
  RNA_expr <-
    data.frame(value=normalized_count_df[,i],
               row.names=row.names(normalized_count_df))
    plot_RNA_group_Tagcov(RNA_expr, i)

}

## continuing
## consider the existence of 5hmC on transcripts
data_dir <- "E:\\20220921-WSL-5hmc\\data\\RNAseq_data"
normalized_count_df <- 
  readRDS(paste0(data_dir,"//normalized_count_df.rds"))

plot_RNA_group_Tagcov <- 
  function(RNA_expr_nor, hmC_file, sample_name){
    ## RNA_expr_nor need to be a dataframe with:
    ## gene_id as rownames, expression value as column (named value)
    
    # get gene ids in each quantile group
    hmC_peak <-
    read.table(hmC_file, 
               header = T,sep = "\t",
               comment.char = "",quote = "\"")
    RNA_expr_nor <- 
      RNA_expr_nor %>%
      filter(!value==0) %>%
      mutate(tet2_sign=ifelse(row.names(.) %in% tet2_targets,"target","not")) %>%
      mutate(hmC_sign=ifelse(row.names(.) %in% hmC_peak$ENSEMBL,"5hmC","no-5hmC"))
    
    p <-
      ggplot(RNA_expr_nor,aes(log2(value),
                              fill=interaction(tet2_sign,hmC_sign)))+
      #geom_boxplot()
      geom_density(alpha=0.7)+
      labs(fill="RNA status")+
      scale_x_continuous(expand = c(0,0),name = "log2(normalized value)")+
      scale_fill_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3"))+
      scale_y_continuous(expand = c(0,0),limits = c(0,0.2))+
      theme_classic()+
      theme(panel.grid = element_blank())
    
    ggsave(p, 
           filename=paste0("E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\",
                           sample_name,"_tet2vsNo_hmC.pdf"),
           width=13, height=8) 
  }
peak_dir <- "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\"
for(i in colnames(normalized_count_df)){
  RNA_expr <-
    data.frame(value=normalized_count_df[,i],
               row.names=row.names(normalized_count_df))
  plot_RNA_group_Tagcov(RNA_expr, 
                        hmC_file = paste0(peak_dir,i,"H_peaks.narrowPeak_anno"),
                        i)
  
}


