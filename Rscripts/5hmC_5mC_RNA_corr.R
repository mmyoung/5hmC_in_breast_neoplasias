## correlation between 5hmC and gene expression
## correlation in different genomic regions in each sample
library(dplyr)
## test in one sample
data_dir <- "E:\\20220921-WSL-5hmc\\data\\RNAseq_data"
normalized_count_df <- 
  readRDS(paste0(data_dir,"//normalized_count_df.rds"))
normalized_count_df <- 
  normalized_count_df %>%
  mutate(gene_id=row.names(.))

head(normalized_count_df)
gene_expr <-
  read.table("E:\\20220921-WSL-5hmc\\data\\RNAseq_data\\UDH_ADH_counts.txt",
             header = T,row.names = 1)
colnames(gene_expr) <- c("J","U3","U5","AW","A2","A10")
gene_anno <- 
  read.table("D:\\biologic_data\\annotation_files\\Hs.GRCh37.87.txt",
             stringsAsFactors = F,header = F, sep = "\t") %>%
  mutate(gen_len=V3-V2)
peak_dir <- "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\"

## A2


normalized_count_df <- 
  normalized_count_df %>%
  left_join(gene_anno[,c("gen_len","V4","V5")], by = c("gene_id"="V4"))
head(normalized_count_df)

normalized_count_df[,1:9] <- 
  normalized_count_df[,1:9] * 1000/normalized_count_df$gen_len

region_5hmC_RNA_corr <- 
  function(sample){
    region_corr_df <- data.frame()
    peak_ano <- 
      read.table(paste0(peak_dir,sample,"H_peaks.narrowPeak_anno"),
                 header = T,sep = "\t",quote = "\"",stringsAsFactors = F)
    peak_ano <- 
    peak_ano %>%
      left_join(normalized_count_df[,c(sample,"gene_id")],
                by = c("ENSEMBL"="gene_id")) %>%
      mutate(quantile_group=ntile(!!sample,3)) %>%
      mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                          ifelse(grepl("^Exon",annotation),"Exon",annotation)))
    
    for(i in unique(peak_ano$group)){
      tmp <- 
        peak_ano %>%
        filter(group==!!i)
      tmp_df <- 
          data.frame(corr = cor.test(tmp$V7,tmp[,sample],method = "spearman",exact = F)$estimate,
       Pvalue  = cor.test(tmp$V7,tmp[,sample], method = "spearman",exact = F)$p.value,
       annotation = i,
        sample = sample
      )
      region_corr_df <- rbind(region_corr_df,tmp_df)
    }
    return(region_corr_df)
  }

res <- 
do.call(rbind,lapply(c("A2","A10","AW","D3","D10","D20","J",
         "U3","U5"),
       region_5hmC_RNA_corr))
write.table(res,
            file = "E:\\20220921-WSL-5hmc\\analysis\\corr_RNA_5hmC_5mC\\Corr_5hmC_RNA_diff_region.txt",
            row.names = F,
            quote = F,
            sep = "\t")


ggplot(res,aes(y=-log10(Pvalue),x=corr,col=annotation))+
  geom_point()


## 2023-4-10
## calculate correlation using only the highest site in each region

region_5hmC_RNA_corr <- 
  function(sample){
    region_corr_df <- data.frame()
    peak_ano <- 
      read.table(paste0(peak_dir,sample,"H_peaks.narrowPeak_anno"),
                 header = T,sep = "\t",quote = "\"",stringsAsFactors = F)
    peak_ano <- 
      peak_ano %>%
      left_join(normalized_count_df[,c(sample,"gene_id")],
                by = c("ENSEMBL"="gene_id")) %>%
      mutate(quantile_group=ntile(!!sample,3)) %>%
      mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                          ifelse(grepl("^Exon",annotation),"Exon",annotation))) %>%
      group_by(group,SYMBOL) %>%
      slice_max(V7) %>%
      ungroup() %>%
      as.data.frame()
    
    for(i in unique(peak_ano$group)){
      tmp <- 
        peak_ano %>%
        filter(group==!!i)
      tmp_df <- 
        data.frame(corr = cor.test(tmp$V7,tmp[,sample],method = "spearman",exact = F)$estimate,
                   Pvalue  = cor.test(tmp$V7,tmp[,sample],method = "spearman",exact = F)$p.value,
                   annotation = i,
                   sample = sample
        )
      region_corr_df <- rbind(region_corr_df,tmp_df)
    }
    return(region_corr_df)
  }
res <- 
  do.call(rbind,lapply(c("A2","A10","AW","D3","D10","D20","J",
                         "U3","U5"),
                       region_5hmC_RNA_corr))
write.table(res,
            file = "E:\\20220921-WSL-5hmc\\analysis\\corr_RNA_5hmC_5mC\\Corr_Highest5hmC_RNA_region.txt",
            row.names = F,
            quote = F,
            sep = "\t")
## calculate correlation using the number of modification sites in each region

region_5hmC_RNA_corr <- 
  function(sample){
    region_corr_df <- data.frame()
    peak_ano <- 
      read.table(paste0(peak_dir,sample,"H_peaks.narrowPeak_anno"),
                 header = T,sep = "\t",quote = "\"",stringsAsFactors = F)
    peak_ano <- 
      peak_ano %>%
      left_join(normalized_count_df[,c(sample,"gene_id")],
                by = c("ENSEMBL"="gene_id")) %>%
      mutate(quantile_group=ntile(!!sample,3)) %>%
      mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                          ifelse(grepl("^Exon",annotation),"Exon",annotation))) %>%
      group_by(group,SYMBOL) %>%
      mutate(peak_num=n()) %>%
      ungroup() %>%
      as.data.frame()
     
    
    for(i in unique(peak_ano$group)){
      tmp <- 
        peak_ano %>%
        filter(group==!!i)
      tmp_df <- 
        data.frame(corr = cor.test(tmp$peak_num,tmp[,sample],method = "spearman",exact = F)$estimate,
                   Pvalue  = cor.test(tmp$peak_num,tmp[,sample], method = "spearman",exact = F)$p.value,
                   annotation = i,
                   sample = sample
        )
      region_corr_df <- rbind(region_corr_df,tmp_df)
    }
    return(region_corr_df)
  }
res <- 
  do.call(rbind,lapply(c("A2","A10","AW","D3","D10","D20","J",
                         "U3","U5"),
                       region_5hmC_RNA_corr))
write.table(res,
            file = "E:\\20220921-WSL-5hmc\\analysis\\corr_RNA_5hmC_5mC\\Corr_NUM5hmC_RNA_region.txt",
            row.names = F,
            quote = F,
            sep = "\t")

#### correlation between 5mC and RNA expression

region_5mC_RNA_corr <- 
  function(sample){
    region_corr_df <- data.frame()
    peak_ano <- 
      read.table(paste0(peak_dir,sample,"M_peaks.narrowPeak_anno"),
                 header = T,sep = "\t",quote = "\"",stringsAsFactors = F)
    peak_ano <- 
      peak_ano %>%
      left_join(normalized_count_df[,c(sample,"gene_id")],
                by = c("ENSEMBL"="gene_id")) %>%
      mutate(quantile_group=ntile(!!sample,3)) %>%
      mutate(group=ifelse(grepl("^Intron",annotation),"Intron",
                          ifelse(grepl("^Exon",annotation),"Exon",annotation))) %>%
      group_by(group,SYMBOL) %>%
      slice_max(V7) %>%
      ungroup() %>%
      as.data.frame()
    
    for(i in unique(peak_ano$group)){
      tmp <- 
        peak_ano %>%
        filter(group==!!i)
      tmp_df <- 
        data.frame(corr = cor.test(tmp$V7,tmp[,sample],method = "spearman",exact = F)$estimate,
                   Pvalue  = cor.test(tmp$V7,tmp[,sample],method = "spearman",exact = F)$p.value,
                   annotation = i,
                   sample = sample
        )
      region_corr_df <- rbind(region_corr_df,tmp_df)
    }
    return(region_corr_df)
  }
res <- 
  do.call(rbind,lapply(c("A2","A10","AW","D3","D10","D20","J",
                         "U3","U5"),
                       region_5mC_RNA_corr))
write.table(res,
            file = "E:\\20220921-WSL-5hmc\\analysis\\corr_RNA_5hmC_5mC\\Corr_Highest5mC_RNA_region.txt",
            row.names = F,
            quote = F,
            sep = "\t")
