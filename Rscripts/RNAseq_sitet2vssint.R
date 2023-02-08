sitet2_RNAseq_FPKM_df <- 
read.table("E:\\20220921-WSL-5hmc\\data\\GSE153250_rnaseq_RPKM.tsv",
  header = T,stringsAsFactors = F,row.names = 1)
sitet2_RNAseq_FPKM_df$gene_id <- 
  gsub("\\..*","",row.names(sitet2_RNAseq_FPKM_df))

wilcoxon_res<- 
  do.call(rbind,apply(sitet2_RNAseq_FPKM_df[,grep("sitet2|sint",colnames(sitet2_RNAseq_FPKM_df))],1,
        function(x){
          tmp_test <- wilcox.test(x[grep("sitet2",names(x))],x[grep("sint",names(x))])
          return(data.frame(p=tmp_test$p.value,
                            sitet2vssint=mean(x[grep("sitet2",names(x))])/mean(x[grep("sint",names(x))])))
          }))
dim(wilcoxon_res)
wilcoxon_res$gene_id <- gsub("\\..*","",row.names(wilcoxon_res))

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
wilcoxon_res_SYMBOL <- bitr(wilcoxon_res$gene_id,
                            fromType = "ENSEMBL",toType = "SYMBOL",
                            OrgDb = org.Hs.eg.db) %>%
  inner_join(.,wilcoxon_res,by=c("ENSEMBL" = "gene_id"))
wilcoxon_res_SYMBOL <-
  merge(wilcoxon_res_SYMBOL,
        sitet2_RNAseq_FPKM_df[,c("gene_id",colnames(sitet2_RNAseq_FPKM_df)[grep("sitet2|sint",colnames(sitet2_RNAseq_FPKM_df))])],
      by.x = "ENSEMBL", by.y="gene_id")
head(wilcoxon_res_SYMBOL)
dim(wilcoxon_res_SYMBOL)
write.table(wilcoxon_res_SYMBOL,
            file = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\mcf7_sitet2vssint_RNAseq_wilcox.txt",
            quote = F,
            sep = "\t",
            row.names = F)
  