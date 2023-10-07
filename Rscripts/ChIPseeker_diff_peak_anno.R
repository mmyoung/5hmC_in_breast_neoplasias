library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(org.Hs.eg.db)
library(reshape2)
library(dplyr)
library(ChIPpeakAnno)
library(GenomicRanges)
library(ggplot2)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## restore the diff peaks in clean format
new_peak_path <- "E:\\20220921-WSL-5hmc\\data\\diff_peak_new"
for(diff_peak in list.files("E:\\20220921-WSL-5hmc\\data\\diff_peak",pattern = "*deseq2.txt$",full.names = TRUE)){
  BaseName <- basename(diff_peak)
  peak <- read.table(diff_peak,row.names = 1)
  write.table(peak,file = paste0(new_peak_path,"\\",BaseName),quote = F, row.names = F,sep = "\t")
}

## annotate the diff peak analysis results
out_path <- "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno"
for(sample_peak in list.files("E:\\20220921-WSL-5hmc\\data\\diff_peak_new",pattern = "*.txt$",full.names = TRUE)){
  #peak <- read.table(sample_peak,row.names = 1)
  tmp <- readPeakFile(sample_peak,header=T)
  BaseName <- basename(sample_peak)
  peakAnno.edb <- 
    annotatePeak(tmp, tssRegion = c(-2000, 2000),
                 genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
                 #ignoreUpstream = TRUE,ignoreDownstream = TRUE,
                 TxDb = txdb, annoDb = "org.Hs.eg.db")
  write.table(as.data.frame(peakAnno.edb),file = paste0(out_path,"\\",BaseName,"_anno"),quote = F,row.names = F,sep = "\t")
}

table(str_split(as.data.frame(peakAnno.edb)[,"annotation"]," \\(",simplify = T)[,1])

## count the numbers for each region in differential methylated regions
peak_region_count <- data.frame() 
for (diff_res in list.files("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno",
                            pattern = "anno$",full.names = TRUE)){
  base_name <- paste(str_split(basename(diff_res),"_",simplify = T)[1:3],collapse = "-")
  diff_res_df <- read.table(diff_res,header = T,sep = "\t",comment.char = "",quote = "\"")
  tmp_df <-
    diff_res_df %>%
    mutate(diff_status=ifelse(p.value<0.05 & Fold>1,"P<0.05&Fold>1",
                              ifelse(p.value<0.05 & Fold<(-1), "P<0.05&Fold<-1","NS")))%>%
    filter(!diff_status=="NS") %>%
    mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                        ifelse(grepl("^Exon",annotation),"Exon",annotation))) %>%
    group_by(diff_status,group) %>%
    summarise(n()) %>%
    as.data.frame() %>%
    mutate(comparison=!!base_name)
  peak_region_count <- rbind(peak_region_count,tmp_df)
}

##calculate the number of DhMRs and DMRs in each stage
peak_region_count <- data.frame()
for (diff_res in list.files("E:\\20220921-WSL-5hmc\\data\\diff_peak_new",
                            pattern = "diff_peak_deseq2.txt$",full.names = TRUE)){
  base_name <- paste(str_split(basename(diff_res),"_",simplify = T)[1:3],collapse = "-")
  diff_res_df <- read.table(diff_res,header = T,sep = "\t",comment.char = "",quote = "\"")
  tmp_df <-
    diff_res_df %>%
    mutate(diff_status=ifelse(p.value<0.05 & Fold>1,"P<0.05&Fold>1",
                              ifelse(p.value<0.05 & Fold<(-1), "P<0.05&Fold<-1","NS")))%>%
    filter(!diff_status=="NS") %>%
    group_by(diff_status) %>%
    summarise(n()) %>%
    as.data.frame() %>%
    mutate(comparison=!!base_name)
  peak_region_count <- rbind(peak_region_count,tmp_df)
}
write.table(peak_region_count,
            file = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\diff_peak_counts.txt",
            quote = F,row.names = F,sep = "\t")

DEG_count <- data.frame()
for (diff_res in list.files("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno",
                            pattern = "txt_anno$",full.names = TRUE)){
  base_name <- paste(str_split(basename(diff_res),"_",simplify = T)[1:3],collapse = "-")
  diff_res_df <- read.table(diff_res,header = T,sep = "\t",comment.char = "",quote = "\"")
  tmp_df <-
    diff_res_df %>%
    mutate(diff_status=ifelse(p.value<0.05 & Fold>1,"P<0.05&Fold>1",
                              ifelse(p.value<0.05 & Fold<(-1), "P<0.05&Fold<-1","NS")))%>%
    filter(!diff_status=="NS") %>%
    group_by(diff_status) %>%
    filter(!ENSEMBL=="NA") %>% 
    filter(!duplicated(ENSEMBL)) %>%
    summarise(n()) %>%
    as.data.frame() %>%
    mutate(comparison=!!base_name)
  DEG_count <- rbind(DEG_count,tmp_df)
}

write.table(DEG_count,
            file = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\DMG_counts.txt",
            quote = F,row.names = F,sep = "\t")



## plot the average plot for differential peaks
setwd("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno")
for(diff_peak in list.files("E:\\20220921-WSL-5hmc\\data\\diff_peak_new",pattern = "*deseq2.txt$",full.names = TRUE)){
  BaseName <- paste(str_split(basename(diff_peak),"_",simplify = T)[1:3],collapse = "-")
  
  tmp_up_df <- read.table(diff_peak,header = T,stringsAsFactors = F) %>%
    filter(Fold>1 & p.value<0.05) %>%
    makeGRangesFromDataFrame()
  tmp_down_df <- read.table(diff_peak,header = T,stringsAsFactors = F) %>%
    filter(Fold<(-1) & p.value<0.05) %>%
    makeGRangesFromDataFrame()
  promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
  tagMatrix_up <- getTagMatrix(tmp_up_df, windows=promoter)
  tagMatrix_down <- getTagMatrix(tmp_down_df, windows=promoter)
  
  if(!is.null(ncol(tagMatrix))){
    out_file_name <- paste(c(BaseName,"_UP_","tagAvgProf"),collapse = "")
    pdf(paste0(out_file_name,".pdf"),width = 10,height = 8)
    plot(plotAvgProf(tagMatrix_up, xlim=c(-3000, 3000),
                xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency"))
    dev.off()
    ## plot 
    out_file_name <- paste(c(BaseName,"_DOWN_","tagAvgProf"),collapse = "")
    pdf(paste0(out_file_name,".pdf"),width = 10,height = 8)
    plot(plotAvgProf(tagMatrix_down, xlim=c(-3000, 3000),
                     xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency"))
    dev.off()
    
  }
}
## plot the up and down peaks together
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

for(i in list.files("E:\\20220921-WSL-5hmc\\data\\diff_peak_new",
                    pattern = "*.deseq2.txt$",
                    full.names = T)){
  BaseName <- paste(str_split(basename(i),"_",simplify = T)[1:3],collapse = "_")
  plot_up_down_avg(diff_peak = i,
                   out_path = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\DMRs_avg_plots",
                   anno_name = BaseName)
}


############################### make plot for the counts of differential peaks in each region #############################################
library(stringr)
diff_peak_region_counts <- 
read.table("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\diff_peak_region_counts.txt",
           header=T, stringsAsFactors = F, comment.char = "",sep = "\t",quote = "\"")
diff_peak_region_counts$methyl_group <- str_split(diff_peak_region_counts$comparison,"-",simplify = T)[,3]
diff_peak_region_counts$comparison_group <- apply(str_split(diff_peak_region_counts$comparison,"-",simplify = T)[,c(2,1)], 1, function(x) {paste(x,collapse="vs")})

library(ggplot2)
library(RColorBrewer)
library(ggsci)
## methylation 
p <- 
ggplot(subset(diff_peak_region_counts,methyl_group=="M"),aes(x=comparison_group,y=n..,fill=group))+
  geom_bar(position = "fill",stat = "identity",alpha=0.8)+
  facet_wrap(.~diff_status)+
  scale_y_continuous(expand = c(0,0.03))+
  labs(x="",y="proportion",title="differential 5mC distribution")+
  scale_fill_nejm()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14))

  ggsave(p,filename = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\5mC_genomic_region_distri_barplot.pdf",width = 10,height = 7)


## hydroxymethylation
p <- 
  ggplot(subset(diff_peak_region_counts,methyl_group=="H"),aes(x=comparison_group,y=n..,fill=group))+
  geom_bar(position = "fill",stat = "identity",alpha=0.8)+
  facet_wrap(.~diff_status)+
  scale_y_continuous(expand = c(0,0.03))+
  labs(x="",y="proportion",title="differential 5hmC distribution")+
  scale_fill_nejm()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5),
        axis.text = element_text(size=10,angle = 30),
        axis.title = element_text(size=12))
  
  ggsave(p,filename = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\5hmC_genomic_region_distri_barplot.pdf",
         width = 10,height = 7)
  

##overlap between DhMRs and DMRs, and correlation between the foldchanges
data_dir <- "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\"
library(ggrastr)
library(ChIPpeakAnno)
DMR_DhMR_ol  <- function(pair_name){
  pair_ls <- str_split(pair_name,"_",simplify = T)
  DMR_df  <-
    read.table(paste0(data_dir,pair_name,"_M_diff_peak_deseq2.txt_anno"),
               stringsAsFactors = F, header = T,
               sep = "\t",comment.char = "",quote = "\"") %>%
    mutate(mC_sign=ifelse(abs(Fold)>1 & p.value<0.05,"diff","no-sig"))%>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
  DhMR_df  <-
    read.table(paste0(data_dir,pair_name,"_H_diff_peak_deseq2.txt_anno"),
               stringsAsFactors = F, header = T,
               sep = "\t",comment.char = "",quote = "\"") %>%
    mutate(hmC_sign=ifelse(abs(Fold)>1 & p.value<0.05,"diff","no-sig"))%>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  pdf(file=paste0(data_dir,pair_name,"_DMR_DhMR_venn.pdf"),
      width = 5,height = 4)
  makeVennDiagram(list(DMR_df,DhMR_df),
                  NameOfPeaks = c("DMR","DhMR"),
                  scaled=F,
                  main=paste0(pair_ls[[1]]," -> ",pair_ls[[2]]))
  dev.off()  
  
  data_for_point =
    data.frame(DhMR_fold=DhMR_df[as.data.frame(findOverlaps(DhMR_df,DMR_df))[["queryHits"]],]$Fold,
               hmC_sign=DhMR_df[as.data.frame(findOverlaps(DhMR_df,DMR_df))[["queryHits"]],]$hmC_sign,
               DMR_fold=DMR_df[as.data.frame(findOverlaps(DhMR_df,DMR_df))[["subjectHits"]],]$Fold,
               mC_sign=DMR_df[as.data.frame(findOverlaps(DhMR_df,DMR_df))[["subjectHits"]],]$mC_sign
               
      )
  
  ## save the DhMRs and DMRs corresponding genes
   data.frame(gene=DhMR_df[as.data.frame(findOverlaps(DhMR_df,DMR_df))[["queryHits"]],]$SYMBOL,
               DhMR_fold=DhMR_df[as.data.frame(findOverlaps(DhMR_df,DMR_df))[["queryHits"]],]$Fold,
               hmC_sign=DhMR_df[as.data.frame(findOverlaps(DhMR_df,DMR_df))[["queryHits"]],]$hmC_sign,
               DMR_fold=DMR_df[as.data.frame(findOverlaps(DhMR_df,DMR_df))[["subjectHits"]],]$Fold,
               mC_sign=DMR_df[as.data.frame(findOverlaps(DhMR_df,DMR_df))[["subjectHits"]],]$mC_sign
               
    ) %>%
    filter(mC_sign=="diff" & hmC_sign=="diff") %>%
    write.table(file = paste0("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\",pair_name,"_DMR_DhMR_gene.txt"),
                quote = F,sep = "\t",row.names = F)
  
  anno_1=sum(-data_for_point$DhMR_fold<0 & -data_for_point$DMR_fold>0 & data_for_point$mC_sign=="diff" & data_for_point$hmC_sign=="diff")
  anno_2=sum(-data_for_point$DhMR_fold>0 & -data_for_point$DMR_fold>0 & data_for_point$mC_sign=="diff" & data_for_point$hmC_sign=="diff")
  anno_3=sum(-data_for_point$DhMR_fold<0 & -data_for_point$DMR_fold<0 & data_for_point$mC_sign=="diff" & data_for_point$hmC_sign=="diff")
  anno_4=sum(-data_for_point$DhMR_fold>0 & -data_for_point$DMR_fold<0 & data_for_point$mC_sign=="diff" & data_for_point$hmC_sign=="diff")
p <-
  ggplot(data_for_point,aes(-DhMR_fold,-DMR_fold))+
    geom_point_rast(color="grey")+
    geom_point_rast(data = subset(data_for_point,mC_sign=="diff" & hmC_sign=="diff"),
               aes(-DhMR_fold,-DMR_fold),color="skyblue")+
  geom_vline(xintercept = c(-1,1),linetype="dashed")+
  geom_hline(yintercept = c(-1,1),linetype="dashed")+
    labs(x=paste0("DhMR log2(",pair_ls[[2]],"/",pair_ls[[1]],")"),
         y=paste0("DMR log2(",pair_ls[[2]],"/",pair_ls[[1]],")"), 
         title = paste0(pair_ls[[1]]," -> ", pair_ls[[2]]))+
    scale_x_continuous(limits = c(-5,5))+
    scale_y_continuous(limits = c(-6,6))+
    annotate("text",x=-4.5,y=4.5,label=anno_1)+
    annotate("text",x=4.5,y=4.5,label=anno_2)+
    annotate("text",x=-4.5,y=-4.5,label=anno_3)+
    annotate("text",x=4.5,y=-4.5,label=anno_4)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5,vjust = 0.5),
          axis.title = element_text(size=11,color="black"),
          axis.text = element_text(size=10,color = "black"))
ggsave(p,filename = paste0(data_dir,pair_name,"_DMR_DhMR_point.pdf"),
       width = 5, height = 4)
}

lapply(c("UDH_ADH","ADH_DCIS"),DMR_DhMR_ol)

## Plot avgplot around TES regions
plot_up_down_avg <- 
  function(diff_peak,out_path,anno_name){
    library(ggplot2)
    promoter <- 
      getBioRegion(TxDb = txdb,
                   upstream = 3000,
                   downstream = 3000,
                   by = "gene",
                   type = "end_site")
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
                         labels = c("-3000bp","-1500bp","TES","1500bp","3000bp"))+
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

for(i in list.files("E:\\20220921-WSL-5hmc\\data\\diff_peak_new",
                    pattern = "*.deseq2.txt$",
                    full.names = T)){
  BaseName <- paste(str_split(basename(i),"_",simplify = T)[1:3],collapse = "_")
  plot_up_down_avg(diff_peak = i,
                   out_path = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\DMRs_TES_avg_plots",
                   anno_name = BaseName)
}

