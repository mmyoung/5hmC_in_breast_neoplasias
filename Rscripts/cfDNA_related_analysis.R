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

## 2023-02-26
##overlap between DhMRs of tumor and cfDNA
library(GenomicRanges)
library(dplyr)
UDH_PT_Fold1 <- 
  read.table("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\UDH_PT_H_diff_peak_deseq2.txt_anno",
           header = T,stringsAsFactors = T, sep = "\t", 
           comment.char = "",quote = "\"") %>%
  filter(Fold>1 & p.value<0.05) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

UDH_PT_FoldM1 <- 
  read.table("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\UDH_PT_H_diff_peak_deseq2.txt_anno",
             header = T,stringsAsFactors = T, sep = "\t", 
             comment.char = "",quote = "\"") %>%
  filter(Fold<(-1) & p.value<0.05) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)


cfDNA_Fold1 <-
  read.table("E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno\\cfDNA_diff_5hmC_ChIPseeker_anno.txt",
           header = T,stringsAsFactors = T,
           sep = "\t", comment.char = "",quote = "\"") %>%
  filter(Fold>1 & p.value<0.05) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

cfDNA_FoldM1 <-
  read.table("E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno\\cfDNA_diff_5hmC_ChIPseeker_anno.txt",
             header = T,stringsAsFactors = T,
             sep = "\t", comment.char = "",quote = "\"") %>%
  filter(Fold<(-1) & p.value<0.05) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

## FindOverlaps
## hypo-hydroxymethylation
OL_Fold1 <- findOverlaps(UDH_PT_Fold1,cfDNA_FoldM1)
length(OL_Fold1)
T_hypo_peak <-
  cbind(as.data.frame(UDH_PT_Fold1)[as.data.frame(OL_Fold1)[["queryHits"]],],
        as.data.frame(cfDNA_FoldM1)[as.data.frame(OL_Fold1)[["subjectHits"]],])
write.table(T_hypo_peak,
            "E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno\\cf_IDC_common_hypo_5hmC.txt",
            row.names = F,
            sep = "\t",quote = F)

## hyper-hydroxymethylation
OL_FoldM1 <- findOverlaps(UDH_PT_FoldM1,cfDNA_Fold1)
length(OL_FoldM1)
T_hyper_peak <-
  cbind(as.data.frame(UDH_PT_FoldM1)[as.data.frame(OL_FoldM1)[["queryHits"]],],
        as.data.frame(cfDNA_Fold1)[as.data.frame(OL_FoldM1)[["subjectHits"]],])
write.table(T_hyper_peak,
            "E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno\\cf_IDC_common_hyper_5hmC.txt",
            row.names = F,
            sep = "\t",quote = F)


library(ChIPpeakAnno)
pdf("E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno\\cf_PT_DhMRs_venn.pdf",
    width = 5,height = 4)
makeVennDiagram(list(UDH_PT_Fold1,cfDNA_Fold1,UDH_PT_FoldM1,cfDNA_FoldM1),
                NameOfPeaks = c("PT_f1","cf_f1","PT_fM1","cf_fM1"))
dev.off()

## hyper-5hmC overlap
pdf("E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno\\cf_PT_hyper-DhMRs_venn.pdf",
    width = 5,height = 4)
makeVennDiagram(list(cfDNA_Fold1,UDH_PT_FoldM1),
                NameOfPeaks = c("cf_f1","PT_fM1"))
dev.off()

## hypo-5hmC overlap
pdf("E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno\\cf_PT_hypo-DhMRs_venn.pdf",
    width = 5,height = 4)
makeVennDiagram(list(cfDNA_FoldM1,UDH_PT_Fold1),
                NameOfPeaks = c("cf_fM1","PT_f1"))
dev.off()


## Plot the distribution of DhMRs near TSS
plot_up_down_avg <- 
  function(diff_peak,out_path,anno_name){
    library(ggplot2)
    promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
    Fold1_peak <- 
      read.table(diff_peak,
                 header = T,stringsAsFactors = F) %>% 
      filter(Fold>1 & p.value<0.05) %>%
      makeGRangesFromDataFrame()
    tagMatrix_up <- getTagMatrix(Fold1_peak, windows=promoter)
    Fold1_ss <- colSums(tagMatrix_up)
    Fold1_ss <- Fold1_ss/sum(Fold1_ss)
    freq_df_1 <- 
      data.frame(Peak_count_frequency=Fold1_ss,
                 pos=seq.int(length(Fold1_ss)),
                 group="Fold1")
    
    FoldM1_peak <- 
      read.table(diff_peak,
                 header = T,stringsAsFactors = F) %>% 
      filter(Fold<(-1) & p.value<0.05) %>%
      makeGRangesFromDataFrame()
    tagMatrix_down <- getTagMatrix(FoldM1_peak, windows=promoter)
    FoldM1_ss <- colSums(tagMatrix_down)
    FoldM1_ss <- FoldM1_ss/sum(FoldM1_ss)
    
    freq_df_M1 <- 
      data.frame(Peak_count_frequency=FoldM1_ss,
                 pos=seq.int(length(FoldM1_ss)),
                 group="FoldM1")
    
    data_for_plot <- rbind(freq_df_1,freq_df_M1)
    p <- 
      ggplot(data_for_plot,aes(pos,Peak_count_frequency,color=group))+
      geom_line()+
      geom_vline(xintercept = 3000,linetype=2)+
      scale_color_manual(values = c("#2a9d8f","#e76f51"))+
      scale_x_continuous(name = "",
                         breaks = c(0,1500,3000,4500,6000),
                         labels = c("-3000bp","-1500bp","TSS","1500bp","3000bp"))+
      scale_y_continuous(name = "peak count frequency",
                         labels = function(x) format(x, scientific = TRUE))+
      labs(title=anno_name)+
      theme_bw()+
      theme(panel.grid = element_blank(),
            axis.text = element_text(colour = "black",size = 10),
            plot.title = element_text(hjust = 0.5,vjust = 0.5,size=12))
    ggsave(p, filename = paste0(out_path,"\\",anno_name,"_Fold1_M1_avgPlot.pdf"),
           width = 5,height = 4)
  }
plot_up_down_avg(diff_peak="E:\\20220921-WSL-5hmc\\data\\cfDNA\\breast_healthy_diff_peak_deseq2.txt",
                 out_path="E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno",
                 anno_name="cfDNA")


## make heatmap for tumor and cfDNA DhMRs

T_hypo_plot_mat <-
  as.matrix(as.data.frame(UDH_PT_Fold1)[as.data.frame(OL_Fold1)[["queryHits"]],c("Conc_UDH_H","Conc_PT_H")])

cfDNA_hypo_plot_mat <- 
as.matrix(as.data.frame(cfDNA_FoldM1)[as.data.frame(OL_Fold1)[["subjectHits"]],c("Conc_breast","Conc_healthy")])

library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(c(0,4,8),c("#435BF7","#F1D8D8","#E26666"))
ht1 = Heatmap(T_hypo_plot_mat, name = "Tumor",
              cluster_rows = F,
              show_row_names = F,
              cluster_columns = F,
              column_labels = c("UDH","IDC"),
              column_names_rot = 0,
              col = col_fun)
ht2 = Heatmap(cfDNA_hypo_plot_mat[,c("Conc_healthy","Conc_breast")], name = "cfDNA",
              cluster_rows = F,
              show_row_names = F,
              cluster_columns = F,
              column_labels = c("healthy","tumor"),
              column_names_rot = 0,
              col = col_fun)
pdf("E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno\\tumor_cfDNA_hypo_5hmC_heatmap.pdf",
    width = 5,height =7)  
draw(ht1+ht2,
     column_title='hypo-5hmC')
dev.off()

### hyper
T_hyper_plot_mat <-
  as.matrix(as.data.frame(UDH_PT_FoldM1)[as.data.frame(OL_FoldM1)[["queryHits"]],c("Conc_UDH_H","Conc_PT_H")])

cfDNA_hyper_plot_mat <- 
  as.matrix(as.data.frame(cfDNA_Fold1)[as.data.frame(OL_FoldM1)[["subjectHits"]],c("Conc_breast","Conc_healthy")])

ht1 = Heatmap(T_hyper_plot_mat, name = "Tumor",
              cluster_rows = F,
              show_row_names = F,
              cluster_columns = F,
              column_labels = c("UDH","IDC"),
              column_names_rot = 0,
              col = col_fun)
ht2 = Heatmap(cfDNA_hyper_plot_mat[,c("Conc_healthy","Conc_breast")], name = "cfDNA",
              cluster_rows = F,
              show_row_names = F,
              cluster_columns = F,
              column_labels = c("healthy","tumor"),
              column_names_rot = 0,
              col = col_fun)
pdf("E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno\\tumor_cfDNA_hyper_5hmC_heatmap.pdf",
    width = 5,height =7)  
draw(ht1+ht2,
     column_title='hyper-5hmC')
dev.off()
