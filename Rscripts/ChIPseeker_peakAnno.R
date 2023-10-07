library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(dplyr)
library(ggrastr)
library(ggplot2)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## annotate the diff peak analysis results
out_path <- "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno"
for(sample_peak in list.files("E:\\20220921-WSL-5hmc\\data\\sample_peak",pattern = "*.narrowPeak$",full.names = TRUE)){
  #peak <- read.table(sample_peak,row.names = 1)
  tmp <- readPeakFile(sample_peak,header=F)
  BaseName <- basename(sample_peak)
  peakAnno.edb <- 
    annotatePeak(tmp, tssRegion = c(-2000, 2000),
                 genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
                 #ignoreUpstream = TRUE,ignoreDownstream = TRUE,
                 TxDb = txdb, annoDb = "org.Hs.eg.db")
  write.table(as.data.frame(peakAnno.edb),file = paste0(out_path,"\\",BaseName,"_anno"),quote = F,row.names = F,sep = "\t")
}

## count the numbers for each region in differential methylated regions
peak_region_count <- data.frame() 
for (diff_res in list.files(out_path,pattern = "anno$",full.names = TRUE)){
  base_name <- str_split(basename(diff_res),"_",simplify = T)[1]
  diff_res_df <- read.table(diff_res,header = T,sep = "\t",comment.char = "",quote = "\"")
  tmp_df <-
    diff_res_df %>%
    mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                        ifelse(grepl("^Exon",annotation),"Exon",annotation))) %>%
    group_by(group) %>%
    summarise(n()) %>%
    as.data.frame() %>%
    mutate(comparison=!!base_name)
  peak_region_count <- rbind(peak_region_count,tmp_df)
}
write.table(peak_region_count,
            file = paste0(out_path,"\\01_peak_region_counts.txt"),
            quote = F,row.names = F,sep = "\t")


## calculate the average enrichment in each genomic region
peak_region_avg <- data.frame() 
for (diff_res in list.files(out_path,pattern = "anno$",full.names = TRUE)){
  base_name <- str_split(basename(diff_res),"_",simplify = T)[1]
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
            file = paste0(out_path,"\\01_peak_region_avg.txt"),
            quote = F,row.names = F,sep = "\t")

##compare the 5hmC and 5mC level in each sample
out_path <- "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno"
sample_5mC_df <- data.frame()
for (diff_res in list.files(out_path,pattern = "M_peaks.narrowPeak_anno$",full.names = TRUE)){
  base_name <- str_split(basename(diff_res),"_",simplify = T)[1]
  diff_res_df <- read.table(diff_res,header = T,sep = "\t",comment.char = "",quote = "\"")
  tmp_df <- data.frame(SignalValue=diff_res_df[,"V7"],Sample=base_name)
  sample_5mC_df <- rbind(sample_5mC_df,tmp_df)
}
sample_5mC_df <- 
sample_5mC_df %>%
  mutate(group=case_when(Sample %in% c("A10M","A2M","AWM") ~ "ADH",
                         Sample %in% c("D10M","D20M","D3M") ~ "DCIS",
                         Sample %in% c("JM","U3M","U5M") ~ "UDH"
  )) %>%
  mutate(group=factor(group,levels=c("UDH","ADH","DCIS")))
library(ggplot2)
p_5mC <- 
ggplot(sample_5mC_df,aes(x=group,y=log1p(SignalValue),fill=group))+
  #geom_violin(width=0.7)+
  geom_boxplot(outlier.color = "NA")+
  labs(title="5mC level",x="")+
  scale_y_continuous(limits = c(0.75,2.5))+
  scale_fill_manual(values = c("UDH"="#1f70a9",
                               "ADH"="#83639f",
                               "DCIS"="#c22f2f"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 12),
        axis.text = element_text(size = 11,colour = "black"),
        legend.position = "None")
p_5mC
ggsave(p_5mC,
       filename = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\5mC_stage_boxplot.pdf",
       width = 5,height = 4)


sample_5hmC_df <- data.frame()
for (diff_res in list.files(out_path,pattern = "H_peaks.narrowPeak_anno$",full.names = TRUE)){
  base_name <- str_split(basename(diff_res),"_",simplify = T)[1]
  diff_res_df <- read.table(diff_res,header = T,sep = "\t",comment.char = "",quote = "\"")
  tmp_df <- data.frame(SignalValue=diff_res_df[,"V7"],Sample=base_name)
  sample_5hmC_df <- rbind(sample_5hmC_df,tmp_df)
}
sample_5hmC_df <- 
  sample_5hmC_df %>%
  mutate(group=case_when(Sample %in% c("A10H","A2H","AWH") ~ "ADH",
                         Sample %in% c("D10H","D20H","D3H") ~ "DCIS",
                         Sample %in% c("JH","U3H","U5H") ~ "UDH",
                         Sample %in% c("PT3H","PT5H","PT6H","PT7H") ~ "PT"
  )) %>%
  mutate(group=factor(group,levels=c("UDH","ADH","DCIS","PT")))
p_5hmC <- 
ggplot(sample_5hmC_df,aes(x=group,y=log1p(SignalValue),fill=group))+
  geom_boxplot(width=0.7,outlier.color = "NA")+
  scale_y_continuous(limits = c(1,2.5))+
  labs(title="5hmC level",x="")+
  scale_fill_manual(values = c("UDH"="#1f70a9",
                               "ADH"="#83639f",
                               "DCIS"="#c22f2f",
                               "PT"="#ffd166"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 12),
        axis.text = element_text(size = 11,colour = "black"),
        legend.position = "None")
ggsave(p_5hmC,
       filename = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\5hmC_stage_boxplot.pdf",
       width = 5,height = 4)
# library(ggpubr)
# p <- ggarrange(p_5mC,p_5hmC,ncol = 2)
# ggsave(p,filename = "E:\\20220921-WSL-5hmc\\analysis\\5hmC_5mC_boxplot.pdf",width = 10,height = 7)

## plot the average plot for all peaks
dir.create("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\average_plot")
setwd("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\average_plot")
for(diff_peak in list.files("E:\\20220921-WSL-5hmc\\data\\sample_peak",full.names = TRUE)){
  BaseName <- basename(diff_peak)
  tmp_peak <- readPeakFile(diff_peak)
  promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
  tagMatrix <- getTagMatrix(tmp_peak, windows=promoter)
  if(!is.null(ncol(tagMatrix))){
    out_file_name <- paste(c(str_split(BaseName,"_",simplify = T)[1],"tagAvgProf"),collapse = "_")
    pdf(paste0(out_file_name,".pdf"),width = 10,height = 8)
    plot(plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
                     xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency"))
    dev.off()
  }
}


## TagHeatmap
get_tagMatrix <-
function(loaded_peak){
  TagMat <- 
    getTagMatrix(peak = loaded_peak, upstream = 3000, downstream = 3000,
                 by = "gene", type = "body", nbin = 800,
                 TxDb = txdb, weightCol = "V5")
  ## need to make sure it's colMeans or colSums/sum(colSums)
  return(colMeans(TagMat))
}

data_path <- "E:\\20220921-WSL-5hmc\\data\\sample_peak"
ADH_peak_merged <-
  do.call(c,lapply(c("A2H","A10H","AWH"),
       function(x){
         readPeakFile(paste0(data_path,"\\",x,"_peaks.narrowPeak"))
       }))
length(ADH_peak_merged)
class(ADH_peak_merged)

UDH_peak_merged <-
  do.call(c,lapply(c("U3H","U5H","JH"),
                   function(x){
                     readPeakFile(paste0(data_path,"\\",x,"_peaks.narrowPeak"))
                   }))
DCIS_peak_merged <-
  do.call(c,lapply(c("D3H","D10H","D20H"),
                   function(x){
                     readPeakFile(paste0(data_path,"\\",x,"_peaks.narrowPeak"))
                   }))
PT_peak_merged <-
  do.call(c,lapply(c("PT3H","PT5H","PT6H","PT7H"),
                   function(x){
                     readPeakFile(paste0(data_path,"\\",x,"_peaks.narrowPeak"))
                   }))


## coordinate corresponds to x-labels
## 0,-3kb; 200,TSS;400,25%;600,50%;800,75%;1000,TTS;1200,+3kb
res <-cbind(get_tagMatrix(UDH_peak_merged),
        get_tagMatrix(ADH_peak_merged),
        get_tagMatrix(DCIS_peak_merged),
        get_tagMatrix(PT_peak_merged))
saveRDS(res,
        file = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\stage_peak_genebody_tagMatrix.rds")
data_for_plot <-
  t(res) %>%
  as.data.frame()
row.names(data_for_plot) <- 
  c("UDH","ADH","DCIS","PT")

library(pheatmap)

pheatmap(as.matrix(data_for_plot/100),
        cluster_cols = F,cluster_rows = F,
        scale = "row",
        show_colnames = F,
        )

rowSums(data_for_plot/100)
plot(res[,1]/100,type="l")
data_for_lineplot <- as.data.frame(res)
colnames(data_for_lineplot) <- c("UDH","ADH","DCIS","PT")
data_for_lineplot$xaxis <- seq(1,dim(data_for_lineplot)[1],by=1)
head(data_for_lineplot)
data_for_lineplot <- reshape2::melt(data_for_lineplot,id="xaxis")
head(data_for_lineplot)
library(ggplot2)
ggplot(data_for_lineplot,aes(x=xaxis,y=value,color=variable))+
  geom_line()

library(ComplexHeatmap)
library(circlize)
tag_mat_res  <-
  readRDS("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\stage_peak_genebody_tagMatrix.rds")
tag_mat_res_t <- t(tag_mat_res)
row.names(tag_mat_res_t) <- c("UDH","ADH","DCIS","PT")
colnames(tag_mat_res_t) <- seq(1,dim(tag_mat_res_t)[2],by=1)
col_fun <- colorRamp2(c(8,40,61),c("#481567FF","#238A8DFF","#DCE319FF"))

ha = HeatmapAnnotation(peak_freq=anno_lines(tag_mat_res,
                             gp = gpar(col = c("#1b9e77","#d95f02","#7570b3","#e7298a")), 
                             add_points = F,
                             size = unit(3, "mm"),
                             height = unit(3, "cm"),
                             
                      )
)
lgd <- 
  Legend(labels = c("UDH","ADH","DCIS","PT"),
         title = "stage",
       legend_gp = gpar(col=c("#1b9e77","#d95f02","#7570b3","#e7298a")),
       type = "lines",
       background = NA
       )
ht <- 
  Heatmap(tag_mat_res_t, 
        name = "mat",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col=col_fun,
        show_column_names = T,
        top_annotation = ha,
        border = T,
        column_labels = c(rep("",39),
                          "-3kb",rep("",199),
                          "TSS",rep("",799),
                          "TES",rep("",199),
                          "+3kb",rep("",40)),
        column_names_rot = 0,
        column_names_centered = T
        
        )

pdf("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\AllStage_5hmC_heamtap.pdf",
    width = 8,height = 6)
draw(ht,annotation_legend_list = list(lgd))
dev.off()

## plot the 5hmC levels in each genomic region separately
library(dplyr)
library(stringr)
peak_FE_df <- data.frame()
for (diff_res in list.files(out_path,pattern = "anno$",full.names = TRUE)){
  base_name <- str_split(basename(diff_res),"_",simplify = T)[1]
  tmp_df <- read.table(diff_res,header = T,sep = "\t",comment.char = "",quote = "\"") %>%
    mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                        ifelse(grepl("^Exon",annotation),"Exon",annotation))) %>%
    mutate(sample=!!base_name) %>%
    select(group,V7,sample)
    
  peak_FE_df <- rbind(peak_FE_df,tmp_df)
}

peak_FE_df <- 
  peak_FE_df %>%
  mutate(stage=case_when(grepl(paste(c("A10","A2","AW"),collapse = "|"),sample) ~ "ADH",
                         grepl(paste(c("D10","D20","D3"),collapse = "|"),sample) ~ "DCIS",
                         grepl(paste(c("J","U3","U5"),collapse = "|"),sample) ~ "UDH",
                         grepl(paste(c("PT3","PT5","PT6","PT7"),collapse = "|"),sample) ~ "IDC"
  )) %>%
  mutate(stage=factor(stage,levels=c("UDH","ADH","DCIS","IDC")))
saveRDS(peak_FE_df,
        file = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\stage_sample_region_FE.rds")
peak_FE_df <- readRDS("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\stage_sample_region_FE.rds")

library(ggplot2)
plot_box_compare <-
  function(tmp_region){
    p<-
  ggplot(peak_FE_df[peak_FE_df$group==tmp_region &
                      grepl("H",peak_FE_df$sample),],
         aes(x=stage,y=V7,fill=stage))+
    geom_boxplot()+
    labs(title=tmp_region)+
    #facet_wrap(group~.,scales = "free")+
    scale_y_continuous(name = "Fold Enrichment",
                       limits = c(1,7),breaks = c(1,3,5,7))+
    scale_x_discrete(name="")+
    scale_fill_manual(values = c("#1f70a9","#83639f","#c22f2f","#ffd166"))+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.position = "NA",
          axis.text = element_text(size = 9,colour = "black"),
          plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 12))
      ggsave(p,filename = paste0("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\Sample_compare_region\\",
                               str_replace(tmp_region,"<|=",""),".pdf"),
             width = 5,height = 4)
}
lapply(unique(peak_FE_df$group),plot_box_compare)

##5mC
library(ggplot2)
library(stringr)
plot_box_compare <-
  function(tmp_region){
    p<-
      ggplot(peak_FE_df[peak_FE_df$group==tmp_region &
                          grepl("M",peak_FE_df$sample),],
             aes(x=stage,y=V7,fill=stage))+
      geom_boxplot()+
      labs(title=tmp_region)+
      #facet_wrap(group~.,scales = "free")+
      scale_y_continuous(name = "Fold Enrichment",
                         limits = c(1,7),breaks = c(1,3,5,7))+
      scale_x_discrete(name="")+
      scale_fill_manual(values = c("#1f70a9","#83639f","#c22f2f"))+
      theme_bw()+
      theme(panel.grid = element_blank(),legend.position = "NA",
            axis.text = element_text(size = 9,colour = "black"),
            plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 12))
    ggsave(p,filename = paste0("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\sample_5mC_compare_region\\",
                               str_replace(tmp_region,"<|=",""),".pdf"),
           width = 5,height = 4)
  }
lapply(unique(peak_FE_df$group),plot_box_compare)


source("E:\\20220921-WSL-5hmc\\src\\.make_merged_5mC_peak_distribution.R")


