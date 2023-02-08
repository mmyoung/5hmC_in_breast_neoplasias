library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(GenomicRanges)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

out_path <- "E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno"
dir.create(out_path)
tmp <- read.table("E:\\20220921-WSL-5hmc\\data\\cfDNA\\breast_healthy_diff_peak_deseq2.txt",
                  header=T,row.names = 1) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
  
peakAnno.edb <- 
  annotatePeak(tmp, tssRegion = c(-2000, 2000),
               genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
               #ignoreUpstream = TRUE,ignoreDownstream = TRUE,
               TxDb = txdb, annoDb = "org.Hs.eg.db")
write.table(as.data.frame(peakAnno.edb),
            file = paste0(out_path,"\\","cfDNA_diff_5hmC_ChIPseeker_anno.txt"),
            quote = F,row.names = F,sep = "\t")

annotation_region_count <- 
  function(annotation_df, sample_id){
    peak_region_count <- data.frame()
    tmp_df <-
      annotation_df %>%
      mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                          ifelse(grepl("^Exon",annotation),"Exon",annotation))) %>%
      group_by(group) %>%
      summarise(n()) %>%
      as.data.frame() %>%
      mutate(comparison=!!sample_id)
    return(tmp_df)
  }

## sample peak annotation
library(stringr)
library(reshape2)
library(dplyr)
peak_region_count <- data.frame()
out_path <- "E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno"

for (sample_peak in list.files("E:\\20220921-WSL-5hmc\\data\\cfDNA",pattern = "*_peaks.narrowPeak",full.names = T)){
  BaseName <- paste(str_split(basename(sample_peak),pattern = "_",simplify = T)[1:2],collapse = "_")
  tmp <- readPeakFile(sample_peak,header=F)
  peakAnno.edb <- 
    annotatePeak(tmp, tssRegion = c(-2000, 2000),
                 genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
                 #ignoreUpstream = TRUE,ignoreDownstream = TRUE,
                 TxDb = txdb, annoDb = "org.Hs.eg.db")
  peakAnno.edb.df <- as.data.frame(peakAnno.edb)
  
  tmp_df  <- annotation_region_count(peakAnno.edb.df,BaseName)
  peak_region_count <- rbind(peak_region_count,tmp_df)
  
  write.table(peakAnno.edb.df,
              file = paste0(out_path,"\\",BaseName,"_anno.txt"),
              quote = F,row.names = F,sep = "\t")
}
write.table(peak_region_count,
            file = paste0(out_path,"\\01_peak_region_counts.txt"),
            quote = F,row.names = F,sep = "\t")

## calculate the average enrichment in each genomic region
peak_region_avg <- data.frame() 
for (diff_res in list.files(out_path,pattern = "^[bh].*anno.txt$",full.names = TRUE)){
  base_name <- paste(str_split(basename(diff_res),"_",simplify = T)[1:2],collapse = "_")
  diff_res_df <- read.table(diff_res,header = T,sep = "\t",comment.char = "",quote = "\"")
  tmp_df <-
    diff_res_df %>%
    mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                        ifelse(grepl("^Exon",annotation),"Exon",annotation))) %>%
    group_by(group) %>%
    summarise_at(vars(V7),list(name = mean)) %>%
    as.data.frame() %>%
    mutate(sample=!!base_name)
  peak_region_avg <- rbind(peak_region_avg,tmp_df)
}
write.table(peak_region_avg,
            file = paste0(out_path,"\\02_peak_region_avg.txt"),
            quote = F,row.names = F,sep = "\t")


## boxplot for the enrichment of all peaks in each sample
library(dplyr)
library(stringr)
cfDNA_dir <- "E:\\20220921-WSL-5hmc\\data\\cfDNA"
tumorDNA_dir <- "E:\\20220921-WSL-5hmc\\data\\sample_peak"

hmC_file_ls <- c(list.files(cfDNA_dir,pattern = "^[bh].*narrowPeak$",full.names = TRUE),
list.files(tumorDNA_dir,pattern = "^[UPJ].*H.*narrowPeak$",full.names = TRUE))

tmp_res <-
  lapply(hmC_file_ls,
       function(file_name){
         anno_file <- read.table(file_name,header = F,comment.char = "",
                                 quote = "\"",stringsAsFactors = F,
                                 sep="\t")
         return(anno_file['V7'])
       })


names(tmp_res) <- 
  unlist(lapply(hmC_file_ls,
       function(file_name){
         BaseName <- str_split(basename(file_name),"_",simplify = T)[1]
       }
))

tmp_res['breast']
names(tmp_res)

tmp <- do.call(rbind,lapply(tmp_res,function(x){data.frame(x)})) %>%
  mutate(sample=gsub("\\..*","",row.names(.))) %>%
  mutate(group=case_when(sample %in% c("JH","U3H","U5H") ~ "UDH",
                         grepl("^PT",sample) ~ "PT",
                         TRUE ~ sample))
unique(tmp$group)

library(ggplot2)
tmp$group <- factor(tmp$group,
                    levels = c("UDH","PT","healthy","breast"),
                    labels = c("UDH","PT","cf-Normal","cf-cancer"))
p <- 
ggplot(tmp,aes(x=group,y=log2(V7),fill=group))+
  geom_violin(adjust=2.5)+
  geom_boxplot(width=0.1)+
  scale_y_continuous(limits = c(0.5,3.8))+
  labs(x="",y="log2(5hmC fold enrichment)")+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "NA",
        axis.text = element_text(size = 12,color = "black"),
        axis.title = element_text(size=14,colour = "black"))
ggsave(p,filename = paste0(out_path,"\\cfDNA_UDH_PT_5hmC_violin_plot.pdf"),width = 5,height = 4)
## compare the mean value of cf-Normal and cf-cancer
mean(tmp[which(tmp$group=="cf-Normal"),"V7"],na.rm=T)-mean(tmp[which(tmp$group=="cf-cancer"),"V7"],na.rm = T)

