library(DESeq2)
library(stringr)
HiSeqV2 <- 
  read.table("D:\\finished_project\\20191213-new_project\\数据分析\\HiSeqV2",
             sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
expr_df_rsem <- round(2^HiSeqV2-1,digits = 0)

clinical_info <- 
  read.table("D:\\finished_project\\20191213-new_project\\数据分析\\BRCA_clinicalMatrix",
             sep = "\t",header = T,stringsAsFactors = F)
clinical_info$sampleID <- str_replace_all(clinical_info$sampleID,"-",".")

expr_df_rsem <-
  expr_df_rsem[,intersect(colnames(HiSeqV2),
                 clinical_info$sampleID[clinical_info$sample_type %in% c("Primary Tumor","Solid Tissue Normal")])]
colData_df <- 
  data.frame(condition=clinical_info[match(colnames(expr_df_rsem),clinical_info$sampleID),"sample_type"],
           row.names = colnames(expr_df_rsem))
colData_df$condition <- factor(colData_df$condition,
                               levels = c("Primary Tumor","Solid Tissue Normal"),
                               labels = c("tumor","normal"))

dds_tmp <- DESeqDataSetFromMatrix(countData=expr_df_rsem,
                                  colData=colData_df, design=~condition)
dds_tmp <- DESeq(dds_tmp)
res_tmp <- results(dds_tmp,contrast=c("condition","tumor","normal"))
na.omit(row.names(res_tmp)[res_tmp$log2FoldChange >1 & res_tmp$pvalue <0.05])


## wilcoxon test for gene expression
table(clinical_info$sample_type)
#Metastatic       Primary Tumor Solid Tissue Normal 
#7                1101                 139 
T_sample <- clinical_info$sampleID[clinical_info$sample_type=="Primary Tumor"]
N_sample <- clinical_info$sampleID[clinical_info$sample_type=="Solid Tissue Normal"]

res <- 
do.call(rbind,apply(HiSeqV2,1,
      function(x) {
  return(data.frame(TN_diff = mean(unlist(x[T_sample]),na.rm=T)/mean(unlist(x[N_sample]),na.rm=T),
           Pvalue = wilcox.test(x[T_sample],x[N_sample])$p.value))
      }
))

Up_gene <- row.names(res)[which(res$Pvalue<0.05 & res$TN_diff >2)]
down_gene <- row.names(res)[which(res$Pvalue<0.05 & res$TN_diff <0.5)]

saveRDS(list(Up_gene,down_gene),
        file = "E:\\20220921-WSL-5hmc\\analysis\\BRCA_NOR_diff_gene\\BRCAvsNOR_diff_gene.rds")


## overlap between continuously changed peak targeted genes with tumor-related genes
T_related_genes <- readRDS("E:\\20220921-WSL-5hmc\\analysis\\BRCA_NOR_diff_gene\\BRCAvsNOR_diff_gene.rds")
continuous_high <- 
  read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuously_high_peak.txt",
           sep="\t", header = T, quote = "\"",stringsAsFactors = F)
continuous_low <- 
  read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuously_low_peak.txt",
             sep="\t", header = T, quote = "\"",stringsAsFactors = F)

##tumor activated gene with continuously high 5hmC gene
library(clusterProfiler)
library(org.Hs.eg.db)
TAG_5hmC_gene_ls <- 
  intersect(T_related_genes[[1]],
            continuous_high$SYMBOL[!is.na(continuous_high$SYMBOL)])
TAG_5hmC_GO <-
  enrichGO(TAG_5hmC_gene_ls,
         OrgDb = org.Hs.eg.db,
         ont = "BP",
         keyType = "SYMBOL")

TRG_5hmC_gene_ls <- 
  intersect(T_related_genes[[2]],
            continuous_low$SYMBOL[!is.na(continuous_low$SYMBOL)])
TRG_5hmC_GO <-
  enrichGO(TRG_5hmC_gene_ls,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           keyType = "SYMBOL")

library(VennDiagram)
library(scales)

venn.diagram(
  x = list(
   TAG = T_related_genes[[1]],
   TRG = T_related_genes[[2]],
   A_5hmC = continuous_high$SYMBOL[!is.na(continuous_high$SYMBOL)]
  ),
  filename = 'E:\\20220921-WSL-5hmc\\analysis\\BRCA_NOR_diff_gene\\continuous_high_venn.tiff',
  output = TRUE ,
  imagetype="tiff" ,
  height = 480 , 
  width = 480 ,
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  rotation = 1
)

venn.diagram(
  x = list(
    TAG = T_related_genes[[1]],
    TRG = T_related_genes[[2]],
    A_5hmC = continuous_low$SYMBOL[!is.na(continuous_low$SYMBOL)]
  ),
  filename = 'E:\\20220921-WSL-5hmc\\analysis\\BRCA_NOR_diff_gene\\continuous_low_venn.tiff',
  output = TRUE ,
  imagetype="tiff" ,
  height = 480 , 
  width = 480 ,
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  rotation = 1
)

