## make heatmap for differential methylation peaks
## but the dataset is too large to make the plot
## So just filter specific peaks in each stage
setwd("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap")
raw_counts_df <- read.delim("raw_counts_mat.tab")
library(stringr)
data_for_plot <- 
  raw_counts_df[,4:16]
colnames(data_for_plot) <- 
str_split(str_replace_all(colnames(raw_counts_df[4:16]),"X\\.",""),"_",simplify = T)[,1]
data_for_plot<-
  data_for_plot[,c("JH","U3H","U5H","A10H","A2H","AWH","D10H","D20H","D3H","PT3H","PT5H","PT6H","PT7H")]
SF <- read.table("scalingFactor.tab",header = T)
row.names(SF) <- str_split(SF$sample,"_",simplify = T)[,1]
SF <- SF[c("JH","U3H","U5H","A10H","A2H","AWH","D10H","D20H","D3H","PT3H","PT5H","PT6H","PT7H"),]

## 1. multiply each row by the scaling factor
new_data <- sweep(as.matrix(data_for_plot),2,SF$scalingFactor,"*") %>% 
  as.data.frame() 
## 2. average among replicates
mean_mat =
  data.frame(UDH_mean=apply(new_data[,c("JH","U3H","U5H")],1,mean),
             ADH_mean=apply(new_data[,c("A10H","A2H","AWH")],1,mean),
             DCIS_mean=apply(new_data[,c("D10H","D20H","D3H")],1,mean),
             PT_mean=apply(new_data[,c("PT3H","PT5H","PT6H","PT7H")],1,mean)
  )
## 3. scale each value by (x-median(x))/mad(x)
run_scale = function(row_ls){
 return((row_ls-median(row_ls))/mad(row_ls))
}

scaled_data = t(apply(mean_mat,1,run_scale))
saveRDS(scaled_data,
        file = "E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\scaled_counts.rds")
summary(apply(scaled_data,1,sd))

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

