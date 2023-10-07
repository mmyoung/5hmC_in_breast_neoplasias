## make heatmap for differential methylation peaks
## but the dataset is too large to make the plot
## So just filter specific peaks in each stage
setwd("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap")
raw_counts_df <- read.delim("raw_counts_mat.tab")
dim(raw_counts_df)
#[1] 108640     16
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(dplyr)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#### Only need to run this in the first time
# data_for_plot <- 
#   raw_counts_df[,4:16]
# colnames(data_for_plot) <- 
# str_split(str_replace_all(colnames(raw_counts_df[4:16]),"X\\.",""),"_",simplify = T)[,1]
# data_for_plot<-
#   data_for_plot[,c("JH","U3H","U5H","A10H","A2H","AWH","D10H","D20H","D3H","PT3H","PT5H","PT6H","PT7H")]
# SF <- read.table("scalingFactor.tab",header = T)
# row.names(SF) <- str_split(SF$sample,"_",simplify = T)[,1]
# SF <- SF[c("JH","U3H","U5H","A10H","A2H","AWH","D10H","D20H","D3H","PT3H","PT5H","PT6H","PT7H"),]
# 
# ## 1. multiply each row by the scaling factor
# new_data <- sweep(as.matrix(data_for_plot),2,SF$scalingFactor,"*") %>% 
#   as.data.frame() 
# ## 2. average among replicates
# mean_mat =
#   data.frame(UDH_mean=apply(new_data[,c("JH","U3H","U5H")],1,mean),
#              ADH_mean=apply(new_data[,c("A10H","A2H","AWH")],1,mean),
#              DCIS_mean=apply(new_data[,c("D10H","D20H","D3H")],1,mean),
#              PT_mean=apply(new_data[,c("PT3H","PT5H","PT6H","PT7H")],1,mean)
#   )
# 
# saveRDS(mean_mat,
#         file = "E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\mean_mat.rds")
# saveRDS(data.frame(raw_counts_df[,1:3],
#                    row.names = paste0("peak",seq.int(dim(raw_counts_df)[1]))),
#         file = "E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\raw_peak_coord.rds")
peak_id<-
  readRDS("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\raw_peak_coord.rds")
colnames(peak_id) <- c("chr","start","end")
mean_mat <- 
  readRDS("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\mean_mat.rds")
## specifically high in one stage
## standard: max/second max >2

row.names(mean_mat) <- 
  paste0("peak",seq.int(dim(mean_mat)[1]))
mean_mat$ratio <-
  unlist(apply(mean_mat,1,function(x){
    max(x)/sort(x,decreasing = T)[2]
  }))
stage_spec_high <- 
  mean_mat %>%
    filter(ratio>2 & !is.infinite(ratio))
dim(stage_spec_high)
## 5871 5

## scale for plot
#scale each value by (x-min(x))/(max(x)-min(x))
run_scale = function(row_ls){
  #return((row_ls-median(row_ls))/mad(row_ls))
  (row_ls-min(row_ls))/(max(row_ls)-min(row_ls))
  }

stage_spec_high_plot = t(apply(stage_spec_high[,1:4],1,run_scale))
summary(stage_spec_high_plot)
## make plot
new_col = colorRamp2(c(-0.2,0.3,1.5),c("#435BF7","#F1D8D8","#E26666"))

pdf("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\stage_spec_peak.pdf",
    width = 5,height =8 )     
Heatmap(stage_spec_high_plot, 
        
        col = new_col,
        cluster_columns = F,
        cluster_rows = T,
        show_row_dend = F,
        show_row_names = F,
        clustering_method_rows = "ward.D",
        column_labels = c("UDH","ADH","DCIS","PT"),
        column_names_rot = 0,
        heatmap_legend_param = list(title="matrix")
)
dev.off()

stage_spec_peak_anno <- 
  cbind(peak_id[row.names(stage_spec_high_plot),],
        stage_spec_high_plot) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    ChIPseeker::annotatePeak(tssRegion = c(-2000, 2000),
                             genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
                             TxDb = txdb, annoDb = "org.Hs.eg.db")%>%
  as.data.frame() %>%
  write.table(file="E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\stage_spec_high_peak.txt",
              quote = F,sep="\t",row.names = F)

## check if the number of stage-specific peaks are right
## down
UDH_ADH <- 
read.table("E:/20220921-WSL-5hmc/data/diff_peak_new/UDH_ADH_H_diff_peak_deseq2.txt",
           header = T, stringsAsFactors = F) %>%
  filter(Fold<(-1) & p.value<0.05) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
UDH_DCIS <- 
  read.table("E:/20220921-WSL-5hmc/data/diff_peak_new/UDH_DCIS_H_diff_peak_deseq2.txt",
             header = T, stringsAsFactors = F) %>%
  filter(Fold<(-1) & p.value<0.05) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
UDH_PT <- 
  read.table("E:/20220921-WSL-5hmc/data/diff_peak_new/UDH_PT_H_diff_peak_deseq2.txt",
             header = T, stringsAsFactors = F) %>%
  filter(Fold<(-1) & p.value<0.05) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
findOverlaps(UDH_ADH[as.data.frame(findOverlaps(UDH_ADH,UDH_DCIS))[["queryHits"]],],
UDH_PT)
##1

##up
UDH_ADH <- 
  read.table("E:/20220921-WSL-5hmc/data/diff_peak_new/UDH_ADH_H_diff_peak_deseq2.txt",
             header = T, stringsAsFactors = F) %>%
  filter(Fold>1& p.value<0.05) %>%
  makeGRangesFromDataFrame()
UDH_DCIS <- 
  read.table("E:/20220921-WSL-5hmc/data/diff_peak_new/UDH_DCIS_H_diff_peak_deseq2.txt",
             header = T, stringsAsFactors = F) %>%
  filter(Fold>1 & p.value<0.05) %>%
  makeGRangesFromDataFrame()
UDH_PT <- 
  read.table("E:/20220921-WSL-5hmc/data/diff_peak_new/UDH_PT_H_diff_peak_deseq2.txt",
             header = T, stringsAsFactors = F) %>%
  filter(Fold>1 & p.value<0.05) %>%
  makeGRangesFromDataFrame()
length(findOverlaps(UDH_ADH[as.data.frame(findOverlaps(UDH_ADH,UDH_DCIS))[["queryHits"]],],
             UDH_PT))
##1456

###up in DCIS stage
ADH_DCIS <- 
  read.table("E:/20220921-WSL-5hmc/data/diff_peak_new/UDH_ADH_H_diff_peak_deseq2.txt",
             header = T, stringsAsFactors = F) %>%
  filter(Fold<(-1) & p.value<0.05) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
UDH_DCIS <- 
  read.table("E:/20220921-WSL-5hmc/data/diff_peak_new/UDH_DCIS_H_diff_peak_deseq2.txt",
             header = T, stringsAsFactors = F) %>%
  filter(Fold<(-1) & p.value<0.05) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
length(as.data.frame(findOverlaps(ADH_DCIS,UDH_DCIS))[["queryHits"]])

## specifically low
mean_mat$low_ratio <-
  unlist(apply(mean_mat,1,function(x){
    sort(x,decreasing = F)[2]/min(x)
  }))
stage_continue_low <- 
  mean_mat %>%
  filter(low_ratio>2 & !is.infinite(low_ratio))
dim(stage_continue_low)
# 1793 5
stage_continue_low_plot = t(apply(stage_continue_low[,1:4],1,run_scale))

new_col = colorRamp2(c(-0.2,0.3,1.5),c("#435BF7","#F1D8D8","#E26666"))
pdf("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\stage_spec_low_peak.pdf",
    width = 5,height =8 ) 
Heatmap(stage_continue_low_plot, 
        col = new_col,
        cluster_columns = F,
        cluster_rows = T,
        show_row_dend = F,
        show_row_names = F,
        clustering_method_rows = "ward.D",
        column_labels = c("UDH","ADH","DCIS","PT"),
        column_names_rot = 0,
        heatmap_legend_param = list(title="matrix")
)
dev.off()
## save stage-specific low peaks
stage_spec_peak_anno <- 
  cbind(peak_id[row.names(stage_continue_low_plot),],
        stage_continue_low_plot) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  ChIPseeker::annotatePeak(tssRegion = c(-2000, 2000),
                           genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
                           TxDb = txdb, annoDb = "org.Hs.eg.db")%>%
  as.data.frame() %>%
  write.table(file="E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\stage_spec_low_peak.txt",
              quote = F,sep="\t",row.names = F)


## continuously high
stage_continue_high <- 
  mean_mat %>%
  filter(ADH_mean > UDH_mean & 
           DCIS_mean > ADH_mean &
           PT_mean > DCIS_mean)
dim(stage_continue_high)
# 2952 5
stage_continue_high_plot = t(apply(stage_continue_high[,1:4],1,run_scale))
pdf("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuous_high_peak.pdf",
    width = 5,height =6 )  
Heatmap(stage_continue_high_plot, 
        col = new_col,
        cluster_columns = F,
        cluster_rows = T,
        show_row_dend = F,
        show_row_names = F,
        clustering_method_rows = "ward.D",
        column_labels = c("UDH","ADH","DCIS","PT"),
        column_names_rot = 0,
        heatmap_legend_param = list(title="matrix")
)
dev.off()

##save the peaks
peak_id[row.names(stage_continue_high),] %>%
  makeGRangesFromDataFrame(seqnames.field = "X..chr.",
                           start.field = "X.start.",
                           end.field = "X.end.") %>%
  annotatePeak(tssRegion = c(-2000, 2000),
               genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
               TxDb = txdb, annoDb = "org.Hs.eg.db") %>%
  as.data.frame() %>%
  write.table(file="E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuously_high_peak.txt",
              quote = F,sep="\t",row.names = F)

## continuously low
stage_continue_low <- 
  mean_mat %>%
  filter(ADH_mean < UDH_mean & 
           DCIS_mean < ADH_mean &
           PT_mean < DCIS_mean)
dim(stage_continue_low)
# 2037 5
stage_continue_low_plot = t(apply(stage_continue_low[,1:4],1,run_scale))
pdf("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuous_low_peak.pdf",
    width = 5,height =6 )
Heatmap(stage_continue_low_plot, 
        col = new_col,
        cluster_columns = F,
        cluster_rows = T,
        show_row_dend = F,
        show_row_names = F,
        clustering_method_rows = "ward.D",
        column_labels = c("UDH","ADH","DCIS","PT"),
        column_names_rot = 0,
        heatmap_legend_param = list(title="matrix")
)
dev.off()

##save the peaks
peak_id[row.names(stage_continue_low),] %>%
  makeGRangesFromDataFrame(seqnames.field = "X..chr.",
                           start.field = "X.start.",
                           end.field = "X.end.") %>%
  annotatePeak(tssRegion = c(-2000, 2000),
               genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
               TxDb = txdb, annoDb = "org.Hs.eg.db") %>%
  as.data.frame() %>%
  write.table(file="E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuously_low_peak.txt",
              quote = F,sep="\t",row.names = F)

## 3. scale each value by (x-median(x))/mad(x)
run_scale = function(row_ls){
 return((row_ls-median(row_ls))/mad(row_ls))
}

scaled_data = t(apply(mean_mat,1,run_scale))
saveRDS(scaled_data,
        file = "E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\scaled_counts.rds")
saveRDS(raw_counts_df[,1:3],
        file = "E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\scaled_counts_peak.rds")
summary(apply(scaled_data,1,sd))

scaled_data <- 
  readRDS(file = "E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\scaled_counts.rds")




peak_group = c(rep("ADH_DCIS_H_fold1",23428),rep("ADH_DCIS_H_foldM1",19695),
rep("DCIS_PT_H_fold1",15171),rep("DCIS_PT_H_foldM1",29648),
rep("UDH_ADH_H_fold1",20353),rep("UDH_ADH_H_foldM1",345))

## keep the peaks with consistent tendency between diff results with counts annotation
keep_peak_idx =
  ((scaled_data[,"ADH_mean"] - scaled_data[,"DCIS_mean"]) >0 & peak_group=="ADH_DCIS_H_fold1") +
  ((scaled_data[,"ADH_mean"] - scaled_data[,"DCIS_mean"]) <0 & peak_group=="ADH_DCIS_H_foldM1") +
  ((scaled_data[,"DCIS_mean"] - scaled_data[,"PT_mean"]) >0 & peak_group=="DCIS_PT_H_fold1") +
  ((scaled_data[,"DCIS_mean"] - scaled_data[,"PT_mean"]) <0 & peak_group=="DCIS_PT_H_foldM1") +
  ((scaled_data[,"UDH_mean"] - scaled_data[,"ADH_mean"]) >0 & peak_group=="UDH_ADH_H_fold1") +
  ((scaled_data[,"UDH_mean"] - scaled_data[,"ADH_mean"]) <0 & peak_group=="UDH_ADH_H_foldM1")
sum(keep_peak_idx[!is.na(keep_peak_idx)])
# 54839
keep_data <- scaled_data[keep_peak_idx==1,]
##
select_idx=apply(keep_data,1,function(x){
  ((sort(x,decreasing = T)[1] - sort(x,decreasing = T)[2])/abs(sort(x,decreasing = T)[1]) >0.8) |
    ((sort(x,decreasing = F)[2] - sort(x,decreasing = F)[1])/abs(sort(x,decreasing = F)[1]) >0.8)  
})
sum(select_idx[!is.na(select_idx)])
#30139

## SD>5
keep_data <- keep_data[apply(keep_data,1,sd)>2,]
dim(keep_data)
#12012
## Remove all NA rows
keep_data[is.na(keep_data)] <- 0
keep_data <- keep_data[!rowSums(keep_data==0),]
dim(keep_data)
# 2187 4

library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(c(-50,-10,50),c("blue","white","red"))
p=Heatmap(keep_data, 
        #col = col_fun,
        cluster_columns = F,
        show_row_dend = F,
        show_row_names = F,
        clustering_method_rows = "ward.D",
        row_km = 8
        )
row_order(p)

## make heatmaps for sequencially changed peaks
## continuously changed peaks filteration
scaled_data
continue_high <- 
  function(row_ls){
    (row_ls[1] < row_ls[2]) & (row_ls[2] < row_ls[3]) & (row_ls[3] < row_ls[4])
    }
up_peaks <- 
  na.omit(scaled_data[apply(scaled_data,1,continue_high),])
dim(up_peaks)
# 2952 4
#new_col = colorRamp2(c(-10,1,10),c("blue","#EFF0E1","#B2182B"))
new_col = colorRamp2(c(-10,5,10),c("blue","#EFF0E1","#B2182B"))

pdf("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuous_up.pdf",
    width = 5,height =8 )     
Heatmap(up_peaks, 
        
        col = new_col,
        cluster_columns = F,
        show_row_dend = F,
        show_row_names = F,
        clustering_method_rows = "ward.D",
        column_labels = c("UDH","ADH","DCIS","PT"),
        column_names_rot = 0,
        heatmap_legend_param = list(title="matrix")
        )
dev.off()

## annotate up_peaks
library(ChIPseeker)
library(GenomicRanges)
library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
up_peak_idx <- apply(scaled_data,1,continue_high)
up_peak_gr <- 
  na.omit(raw_counts_df[up_peak_idx,1:3]) %>% 
  makeGRangesFromDataFrame(seqnames.field = "X..chr.",
                           start.field = "X.start.",
                           end.field = "X.end.")
up_peak_anno <- annotatePeak(up_peak_gr,
                             tssRegion = c(-2000, 2000),
                             genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
                             TxDb = txdb, annoDb = "org.Hs.eg.db")
enrichKEGG(as.data.frame(up_peak_anno)[,"ENSEMBL"],
           organism = "hsa"
         )

#continuously down peaks
continue_down <- 
  function(row_ls){
    (row_ls[1] > row_ls[2]) & (row_ls[2] > row_ls[3]) & (row_ls[3] > row_ls[4])
  }
down_peaks <- 
  na.omit(scaled_data[apply(scaled_data,1,continue_down),])
dim(down_peaks)
# 2037 4
new_col = colorRamp2(c(-10,1,10),c("blue","#EFF0E1","#B2182B"))
pdf("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuous_down.pdf",
    width = 5,height =8 )     
Heatmap(down_peaks, 
        col = new_col,
        cluster_columns = F,
        show_row_dend = F,
        show_row_names = F,
        clustering_method_rows = "ward.D",
        column_labels = c("UDH","ADH","DCIS","PT"),
        column_names_rot = 0,
        heatmap_legend_param = list(title="matrix")
)
dev.off()
## annotate continuously down peaks
down_peak_idx <- apply(scaled_data,1,continue_down)
down_peak_gr <- 
  na.omit(raw_counts_df[down_peak_idx,1:3]) %>% 
  makeGRangesFromDataFrame(seqnames.field = "X..chr.",
                           start.field = "X.start.",
                           end.field = "X.end.")
down_peak_anno <- annotatePeak(down_peak_gr,
                             tssRegion = c(-2000, 2000),
                             genomicAnnotationPriority = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Intergenic"),
                             TxDb = txdb, annoDb = "org.Hs.eg.db")
enrichGO(as.data.frame(down_peak_anno)[,"ENSEMBL"],
         keyType = "ENSEMBL",ont = "BP",
         OrgDb = org.Hs.eg.db)

## 2023-4-4
low_sample_cov <- 
read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuous_low_avg.txt",
           stringsAsFactors = F,header = F) %>%
  mutate(group=case_when(V3 %in% c("A10H","A2H","AWH") ~ "ADH",
                         V3 %in% c("D10H","D20H","D3H") ~ "DCIS",
                         V3 %in% c("JH","U3H","U5H") ~ "UDH",
                         V3 %in% c("PT3H","PT5H","PT6H","PT7H") ~ "PT"
  ))
library(dplyr)
head(low_sample_cov)
low_sample_cov <- 
low_sample_cov %>%
  group_by(V1,group) %>%
  summarise_at(vars(V2),list(name=mean))
head(low_sample_cov)

library(ggplot2)
ggplot(low_sample_cov,aes(V1,name,col=group))+
  geom_line()+
  labs(title="Continuously low peaks")+
  scale_x_continuous(name = "position",
                     breaks = c(0,300,400,700),
                     labels = c("-3kb","Start","End","+3kb"))+
  scale_y_continuous(name = "normalizaed coverage")+
  scale_color_manual(values = c("UDH"="#1f70a9",
                                "ADH"="#83639f",
                                "DCIS"="#c22f2f",
                                "PT"="#ffd166"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))


read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuous_high_avg.txt",
           stringsAsFactors = F,header = F) %>%
  mutate(group=case_when(V3 %in% c("A10H","A2H","AWH") ~ "ADH",
                         V3 %in% c("D10H","D20H","D3H") ~ "DCIS",
                         V3 %in% c("JH","U3H","U5H") ~ "UDH",
                         V3 %in% c("PT3H","PT5H","PT6H","PT7H") ~ "PT"
  )) %>%
  group_by(V1,group) %>%
  summarise_at(vars(V2),list(name=mean)) %>%
  ggplot(aes(V1,name,col=group))+
  geom_line()+
  scale_x_continuous(name = "position",
                     breaks = c(0,300,400,700),
                     labels = c("-3kb","Start","End","+3kb"))+
  scale_y_continuous(name = "normalizaed coverage")+
  scale_color_manual(values = c("UDH"="#1f70a9",
                                "ADH"="#83639f",
                                "DCIS"="#c22f2f",
                                "PT"="#ffd166"))+
  labs(title="Continuously high peaks")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))
## results: the plots doesn't look satisfactory

## 2023-4-24
## Try to divide the IP coverage by Input coverage
library(ggplot2)
#dealt_data <- 
cbind(read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuous_high_avg_IP.txt",
           stringsAsFactors = F,header = F),
      read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuous_high_avg_Input.txt",
            stringsAsFactors = F,header = F)) %>%
  `colnames<-` (c("Coord_IP","IP_cov","IP_s","Coord_Input","Input_cov","Input_s")) %>%
  mutate(group=case_when(IP_s %in% c("A10H","A2H","AWH") ~ "ADH",
                         IP_s %in% c("D10H","D20H","D3H") ~ "DCIS",
                         IP_s %in% c("JH","U3H","U5H") ~ "UDH",
                         IP_s %in% c("PT3H","PT5H","PT6H","PT7H") ~ "PT")) %>%
  group_by(Coord_IP,group) %>%
  summarise_at(vars(IP_cov,Input_cov),mean, na.rm=T) %>%
  mutate(fold=IP_cov-Input_cov) %>%
  ggplot(aes(Coord_IP,fold,col=group))+
  geom_line()+
  scale_x_continuous(name = "position",
                     breaks = c(0,300,400,700),
                     labels = c("-3kb","Start","End","+3kb"))+
  scale_y_continuous(name = "normalizaed coverage")+
  scale_color_manual(values = c("UDH"="#1f70a9",
                                "ADH"="#83639f",
                                "DCIS"="#c22f2f",
                                "PT"="#ffd166"))+
  labs(title="Continuously high peaks")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))

cbind(read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuous_low_avg_IP.txt",
                 stringsAsFactors = F,header = F),
      read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuous_low_avg_Input.txt",
                 stringsAsFactors = F,header = F)) %>%
  `colnames<-` (c("Coord_IP","IP_cov","IP_s","Coord_Input","Input_cov","Input_s")) %>%
  mutate(group=case_when(IP_s %in% c("A10H","A2H","AWH") ~ "ADH",
                         IP_s %in% c("D10H","D20H","D3H") ~ "DCIS",
                         IP_s %in% c("JH","U3H","U5H") ~ "UDH",
                         IP_s %in% c("PT3H","PT5H","PT6H","PT7H") ~ "PT")) %>%
  group_by(Coord_IP,group) %>%
  summarise_at(vars(IP_cov,Input_cov),mean, na.rm=T) %>%
  mutate(fold=IP_cov/Input_cov) %>%
  ggplot(aes(Coord_IP,fold,col=group))+
  geom_line()+
  scale_x_continuous(name = "position",
                     breaks = c(0,300,400,700),
                     labels = c("-3kb","Start","End","+3kb"))+
  scale_y_continuous(name = "normalizaed coverage")+
  scale_color_manual(values = c("UDH"="#1f70a9",
                                "ADH"="#83639f",
                                "DCIS"="#c22f2f",
                                "PT"="#ffd166"))+
  labs(title="Continuously low peaks")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))
