library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(dplyr)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


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
  do.call(c,lapply(c("A2M","A10M","AWM"),
       function(x){
         readPeakFile(paste0(data_path,"\\",x,"_peaks.narrowPeak"))
       }))
length(ADH_peak_merged)
class(ADH_peak_merged)

UDH_peak_merged <-
  do.call(c,lapply(c("U3M","U5M","JM"),
                   function(x){
                     readPeakFile(paste0(data_path,"\\",x,"_peaks.narrowPeak"))
                   }))
DCIS_peak_merged <-
  do.call(c,lapply(c("D3M","D10M","D20M"),
                   function(x){
                     readPeakFile(paste0(data_path,"\\",x,"_peaks.narrowPeak"))
                   }))


## coordinate corresponds to x-labels
## 0,-3kb; 200,TSS;400,25%;600,50%;800,75%;1000,TTS;1200,+3kb
res <-cbind(get_tagMatrix(UDH_peak_merged),
        get_tagMatrix(ADH_peak_merged),
        get_tagMatrix(DCIS_peak_merged))
saveRDS(res,
        file = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\stage_5mCpeak_genebody_tagMatrix.rds")
data_for_plot <-
  t(res) %>%
  as.data.frame()
row.names(data_for_plot) <- 
  c("UDH","ADH","DCIS")


library(ComplexHeatmap)
library(circlize)
tag_mat_res  <-
  readRDS("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\stage_5mCpeak_genebody_tagMatrix.rds")
tag_mat_res_t <- t(tag_mat_res)
row.names(tag_mat_res_t) <- c("UDH","ADH","DCIS")
colnames(tag_mat_res_t) <- seq(1,dim(tag_mat_res_t)[2],by=1)
col_fun <- colorRamp2(c(8,40,61),c("#481567FF","#238A8DFF","#DCE319FF"))

ha = HeatmapAnnotation(peak_freq=anno_lines(tag_mat_res,
                             gp = gpar(col = c("#1b9e77","#d95f02","#7570b3")), 
                             add_points = F,
                             size = unit(3, "mm"),
                             height = unit(3, "cm"),
                             
                      )
)
lgd <- 
  Legend(labels = c("UDH","ADH","DCIS"),
         title = "stage",
       legend_gp = gpar(col=c("#1b9e77","#d95f02","#7570b3")),
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

pdf("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\AllStage_5mC_heamtap.pdf",
    width = 8,height = 6)
draw(ht,annotation_legend_list = list(lgd))
dev.off()