
library(GenomicRanges)
library(GenomicFeatures)
library(forcats)
library(ggplot2)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
## evaluate the correlation between DMR and DhMRs
## UDH vs ADH
DMR_dir <- "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\"

DMR_DhmR_plot <-
  function(x_file,y_file,out_file,prefix,x_name,y_name){
  DhMRs_peak <- 
    read.table(y_file,
               stringsAsFactors = F, header = T,
               sep = "\t",comment.char = "",quote = "\"")
  DMRs_peak <- 
    read.table(x_file,
               stringsAsFactors = F, header = T,
               sep = "\t",comment.char = "",quote = "\"")
  
  DMRs_peak$diff_sign <- 
    ifelse(DMRs_peak$Fold>1 & DMRs_peak$p.value<0.05,"down",
           ifelse(DMRs_peak$Fold <(-1) & DMRs_peak$p.value<0.05,"up","no-sig"))
  DMRs_peak$hmC_fold <- NA
  
  DMRs_peak$hmC_fold[as.data.frame(findOverlaps(makeGRangesFromDataFrame(DMRs_peak),makeGRangesFromDataFrame(DhMRs_peak)))[["queryHits"]]] <-
    DhMRs_peak[as.data.frame(findOverlaps(makeGRangesFromDataFrame(DMRs_peak),makeGRangesFromDataFrame(DhMRs_peak)))[["subjectHits"]],"Fold"]
  DMRs_peak$hmC_sign <-
    ifelse(DMRs_peak$hmC_fold>0,"pos","neg")
  
  #DMRs_peak$diff_sign <-
  #  factor(DMRs_peak$diff_sign,
  #         levels = c("up","down","no-sig"),
  #         labels = c("up","down","no-change"))
  
  for(i in c("neg","pos")){
    
    p <-
      DMRs_peak %>%
      filter(!is.na(hmC_fold)) %>%
      filter(hmC_sign==!!i) %>%
      filter(!diff_sign=="no-sig") %>%
      mutate(diff_sign=factor(diff_sign,levels = c("up","down"))) %>%
      ggplot(aes(diff_sign,-hmC_fold))+
      geom_violin(aes(fill=diff_sign))+
      geom_boxplot(width=0.1)+
      labs(title = prefix)+
      scale_y_continuous(name = y_name)+
      scale_x_discrete(name=x_name)+
      #scale_fill_manual(values = c("#ee6c4d","#3d5a80","#adb5bd"))+
      scale_fill_manual(values = c("#ee6c4d","#3d5a80"))
    
    (p + theme(panel.grid.major = element_blank(),
               panel.background = element_rect(fill="white",colour = "black"),
               #axis.line = element_blank(),
               panel.grid.minor = element_blank(),
               axis.text = element_text(size = 10,colour = "black"),
               plot.title = element_text(hjust = 0.5,vjust = 0.5),
               legend.position = "none") 
    ) %>%
      ggsave(paste0(out_file,prefix,"_",i,"_boxplot.pdf"),.,width = 5,height = 4)
  }
  
}

setwd("E:\\20220921-WSL-5hmc\\analysis\\5hmC_5mC_correlation")
DMR_DhmR_plot(paste0(DMR_dir,"UDH_ADH_M_diff_peak_deseq2.txt_anno"),
              paste0(DMR_dir,"UDH_ADH_H_diff_peak_deseq2.txt_anno"),
              "E:\\20220921-WSL-5hmc\\analysis\\5hmC_5mC_correlation\\DMR-DhMRs_",
              prefix="UDH_ADH",
              x_name="5mC change",
              y_name="log2(5hmC Fold)")
DMR_DhmR_plot(paste0(DMR_dir,"ADH_DCIS_M_diff_peak_deseq2.txt_anno"),
              paste0(DMR_dir,"ADH_DCIS_H_diff_peak_deseq2.txt_anno"),
              "E:\\20220921-WSL-5hmc\\analysis\\5hmC_5mC_correlation\\DMR-DhMRs_",
              prefix="DCIS_ADH",
              x_name="5mC change",
              y_name="log2(5hmC Fold)")
DMR_DhmR_plot(paste0(DMR_dir,"UDH_ADH_H_diff_peak_deseq2.txt_anno"),
              paste0(DMR_dir,"UDH_ADH_M_diff_peak_deseq2.txt_anno"),
              "E:\\20220921-WSL-5hmc\\analysis\\5hmC_5mC_correlation\\DhMR-DMRs_",
              prefix="ADH_UDH",
              x_name="5hmC change",
              y_name="log2(5mC Fold)")
DMR_DhmR_plot(paste0(DMR_dir,"ADH_DCIS_H_diff_peak_deseq2.txt_anno"),
              paste0(DMR_dir,"ADH_DCIS_M_diff_peak_deseq2.txt_anno"),
              "E:\\20220921-WSL-5hmc\\analysis\\5hmC_5mC_correlation\\DhMR-DMRs_",
              prefix="DCIS_ADH",
              x_name="5hmC change",
              y_name="log2(5mC Fold)")

## make the violin plots for DhMRs and DMRs together
# DMRs_peak_1 <- 
#   read.table(paste0(DMR_dir,"UDH_ADH_M_diff_peak_deseq2.txt_anno"),
#              stringsAsFactors = F, header = T,
#              sep = "\t",comment.char = "",quote = "\"") 
# DhMRs_peak_1 <- 
#   read.table(paste0(DMR_dir,"UDH_ADH_H_diff_peak_deseq2.txt_anno"),
#              stringsAsFactors = F, header = T,
#              sep = "\t",comment.char = "",quote = "\"")


DMRs_peak_1 <-
  read.table(paste0(DMR_dir,"ADH_DCIS_M_diff_peak_deseq2.txt_anno"),
             stringsAsFactors = F, header = T,
             sep = "\t",comment.char = "",quote = "\"")
DhMRs_peak_1 <-
  read.table(paste0(DMR_dir,"ADH_DCIS_H_diff_peak_deseq2.txt_anno"),
             stringsAsFactors = F, header = T,
             sep = "\t",comment.char = "",quote = "\"")

DMRs_peak_1 %>%
  dplyr::slice(as.data.frame(findOverlaps(makeGRangesFromDataFrame(DMRs_peak_1),makeGRangesFromDataFrame(DhMRs_peak_1)))[["queryHits"]]) %>%
  mutate(group="5mC") %>%
  mutate(other_Fold = -DhMRs_peak_1[as.data.frame(findOverlaps(makeGRangesFromDataFrame(DMRs_peak_1),makeGRangesFromDataFrame(DhMRs_peak_1)))[["subjectHits"]],"Fold"]) %>%
  mutate(other_P = DhMRs_peak_1[as.data.frame(findOverlaps(makeGRangesFromDataFrame(DMRs_peak_1),makeGRangesFromDataFrame(DhMRs_peak_1)))[["subjectHits"]],"p.value"]) %>%
  mutate(change=case_when(Fold<(-1) & p.value <0.05 ~"up",
                             Fold>1 & p.value <0.05 ~"down",
                             .default =  "no-sig")) %>%
  dplyr::filter(abs(other_Fold)>1 & other_P<0.05) %>%
  dplyr::select(group,other_Fold,change) -> data_for_plot1
table(data_for_plot1$change)
#down no-sig 
#332  16401 

DhMRs_peak_1 %>%
  dplyr::slice(as.data.frame(findOverlaps(makeGRangesFromDataFrame(DMRs_peak_1),makeGRangesFromDataFrame(DhMRs_peak_1)))[["subjectHits"]]) %>%
  mutate(group="5hmC") %>%
  mutate(change=case_when(Fold<(-1) & p.value <0.05 ~"up",
                             Fold>1 & p.value <0.05 ~"down",
                             .default =  "no-sig")) %>%
  mutate(other_Fold = -DMRs_peak_1[as.data.frame(findOverlaps(makeGRangesFromDataFrame(DMRs_peak_1),makeGRangesFromDataFrame(DhMRs_peak_1)))[["queryHits"]],"Fold"]) %>%
  mutate(other_P = DMRs_peak_1[as.data.frame(findOverlaps(makeGRangesFromDataFrame(DMRs_peak_1),makeGRangesFromDataFrame(DhMRs_peak_1)))[["queryHits"]],"p.value"]) %>%
  dplyr::filter(abs(other_Fold)>1 & other_P<0.05) %>%
  dplyr::select(group,other_Fold, change) -> data_for_plot2
table(data_for_plot2$change)
#  down no-sig 
#332   2503 
  
p <- 
rbind(data_for_plot1,data_for_plot2) %>%
  #filter(!change=="no-sig") %>%
  mutate(change_b=ifelse(change %in% c("up","down"),"change","no")) %>%
  ggplot(aes(x=group,y=abs(other_Fold)))+
  geom_violin(aes(fill=change_b),
              position = position_dodge(width = 0.6),
              width=0.5)+
  geom_boxplot(aes(fill=change_b),width=0.1,
               position=position_dodge(width = 0.6),
               outlier.colour = "NA")+
  geom_segment(aes(x=0.75,y=5,xend=1.25,yend=5))+
  annotate("text",x = 1,y = 5.3,label="5mC change")+
  geom_segment(aes(x=1.75,y=5,xend=2.25,yend=5))+
  annotate("text",x = 2,y = 5.3,label="5hmC change")+
  #labs(title = "ADH vs UDH")+
  labs(title = "DCIS vs ADH")+
  scale_y_continuous(limits = c(1,6),name = "abs(log2FoldChange)")+
  scale_x_discrete(name="")+
  scale_fill_manual(values = c("no"="#adb5bd","change"="#ee6c4d"))+
  #scale_fill_manual(values = c("up"="#ee6c4d","down"="#3d5a80","no-sig"="#adb5bd"))+
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill="white",colour = "black"),
        #axis.line = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12,colour = "black"),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))
#ggsave(p, filename = "E:\\20220921-WSL-5hmc\\analysis\\5hmC_5mC_correlation\\ADH_UDH_DMR_DhMR_corr_violin.pdf",
ggsave(p, filename = "E:\\20220921-WSL-5hmc\\analysis\\5hmC_5mC_correlation\\DCIS_ADH_DMR_DhMR_corr_violin.pdf",
      width = 5,height = 4)

## coverage plot of 5hmC in DMR regions
library(ChIPseeker)
## ADH vs UDH

DMR_DhMR_avg_plot <-
  function(prefix,out_file,title){
    UDH_ADH_DMR <- 
    read.table(paste0(DMR_dir,prefix,"_M_diff_peak_deseq2.txt_anno"),
               stringsAsFactors = F, header = T,
               sep = "\t",comment.char = "",quote = "\"") %>%
    filter(abs(Fold) >1 & p.value<0.05) %>%
    mutate(type="gene") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    makeBioRegionFromGranges(by = "peak",
                             type = "body",
                             upstream = 3000, downstream = 3000)
  
  UDH_ADH_DhMR <- 
    read.table(paste0(DMR_dir,prefix,"_H_diff_peak_deseq2.txt_anno"),
               stringsAsFactors = F, header = T,
               sep = "\t",comment.char = "",quote = "\"") %>%
    filter(abs(Fold) >1 & p.value<0.05) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  Tag_mat <- 
    getTagMatrix(UDH_ADH_DhMR,windows = UDH_ADH_DMR,
                 nbin=200)
  (plotAvgProf(Tag_mat,xlim = c(1,200),
              xlab="DMRs region (5'->3')", 
              ylab = "Read Count Frequency")+
    labs(title = paste0("DhMRs flanking DMRs (",title,")"))+
    scale_x_continuous(breaks = c(0,50,100,150,200),
                       labels = c("-3kb","start","center","end","+3kb"))+
    theme(plot.title = element_text(hjust = 0.5,vjust = 0.5),
          axis.text = element_text(size = 10, colour = "black"))) %>%
    ggsave(filename = out_file,.,width = 5,height = 4)
  
  }

DMR_DhMR_avg_plot(prefix="UDH_ADH",
                  out_file="E:\\20220921-WSL-5hmc\\analysis\\5hmC_5mC_correlation\\ADHvsUDH_DMR_DhMR_avg_plot.pdf",
                  title="ADH vs UDH")
DMR_DhMR_avg_plot(prefix="ADH_DCIS",
                  out_file="E:\\20220921-WSL-5hmc\\analysis\\5hmC_5mC_correlation\\DCISvsADH_DMR_DhMR_avg_plot.pdf",
                  title="DCIS vs ADH")

## coverage of 5mC flanking DhMRs
DMR_DhMR_avg_plot <-
  function(prefix,out_file,title){
    UDH_ADH_DMR <- 
      read.table(paste0(DMR_dir,prefix,"_H_diff_peak_deseq2.txt_anno"),
                 stringsAsFactors = F, header = T,
                 sep = "\t",comment.char = "",quote = "\"") %>%
      filter(abs(Fold) >1 & p.value<0.05) %>%
      mutate(type="gene") %>%
      makeGRangesFromDataFrame(keep.extra.columns = T) %>%
      makeBioRegionFromGranges(by = "peak",
                               type = "body",
                               upstream = 3000, downstream = 3000)
    
    UDH_ADH_DhMR <- 
      read.table(paste0(DMR_dir,prefix,"_M_diff_peak_deseq2.txt_anno"),
                 stringsAsFactors = F, header = T,
                 sep = "\t",comment.char = "",quote = "\"") %>%
      filter(abs(Fold) >1 & p.value<0.05) %>%
      makeGRangesFromDataFrame(keep.extra.columns = T)
    
    Tag_mat <- 
      getTagMatrix(UDH_ADH_DhMR,windows = UDH_ADH_DMR,
                   nbin=200)
    (plotAvgProf(Tag_mat,xlim = c(1,200),
                 xlab="DhMRs region (5'->3')", 
                 ylab = "Read Count Frequency")+
        labs(title = paste0("DMRs flanking DhMRs (",title,")"))+
        scale_x_continuous(breaks = c(0,50,100,150,200),
                           labels = c("-3kb","start","center","end","+3kb"))+
        theme(plot.title = element_text(hjust = 0.5,vjust = 0.5),
              axis.text = element_text(size = 10, colour = "black"))) %>%
      ggsave(filename = out_file,.,width = 5,height = 4)
    
  }

DMR_DhMR_avg_plot(prefix="UDH_ADH",
                  out_file="E:\\20220921-WSL-5hmc\\analysis\\5hmC_5mC_correlation\\ADHvsUDH_DhMR_DMR_avg_plot.pdf",
                  title="ADH vs UDH")
DMR_DhMR_avg_plot(prefix="ADH_DCIS",
                  out_file="E:\\20220921-WSL-5hmc\\analysis\\5hmC_5mC_correlation\\DCISvsADH_DhMR_DMR_avg_plot.pdf",
                  title="DCIS vs ADH")

