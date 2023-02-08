## correlation between RNA expression with 5hmC level
## split RNA into groups according to their expression quantile value

## load required libraries
library(dplyr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(tidyverse)

## load RNA expression data
data_dir <- "E:\\20220921-WSL-5hmc\\data\\RNAseq_data"
# ADH_DCIS_count <-
#   read.table(paste0(data_dir,"\\ADH_DCIS_counts.txt"),
#              row.names = 1,header = T) %>%
#   filter(!rowSums(.)==0)
# UDH_ADH_count <-
#   read.table(paste0(data_dir,"\\UDH_ADH_counts.txt"),
#              row.names = 1,header = T) %>%
#   filter(!rowSums(.)==0)
# 
# ## normalize counts
# ADH_DCIS_count_nor <- 
#   sweep(ADH_DCIS_count, 2, colSums(ADH_DCIS_count), "/") * 10e6
# colnames(ADH_DCIS_count_nor) <- c("AW","A2","A10","D3","D10","D20")
# 
# ## normalize counts
# UDH_ADH_count_nor <- 
#   sweep(UDH_ADH_count, 2, colSums(UDH_ADH_count), "/") * 10e6
# colnames(UDH_ADH_count_nor) <- c("J","U3","U5","AW","A2","A10")
# 
# UDH_ADH_count_nor <- UDH_ADH_count_nor[,c("J","U3","U5")]
# 
# normalized_count_df <- 
#   merge(ADH_DCIS_count_nor, UDH_ADH_count_nor, by=0) %>%
#   column_to_rownames(var="Row.names")
# 
# saveRDS(normalized_count_df,
#         file=paste0(data_dir,"//normalized_count_df.rds"))

## Get needed data for peak annotation
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

genebody <- getBioRegion(TxDb = txdb,
                         by = "gene",
                         type = "body")

## function to get the tagMatrix for each group of peaks
get_TagMat <- 
  function(group,peak_file){
    first_TagMat <- 
      getTagMatrix(GRanges(peak_file[which(peak_file[,"group"]==group),]),
                   windows = genebody, nbin = 800,
                   upstream = 1000,downstream = 1000)
    ss <- colSums(first_TagMat)
    ss <- ss/sum(ss)
    freq_df <- 
      data.frame(pos=seq.int(1,length(ss),by=1),
                 Peak_count_freq=ss,
                 group=group)
    return(freq_df)
  }

plot_RNA_group_Tagcov <- 
  function(RNA_expr_nor, annotated_peak, sample_name){
    ## RNA_expr_nor need to be a dataframe with:
    ## gene_id as rownames, expression value as column (named value)
    
    # get gene ids in each quantile group
    RNA_expr_nor <- 
      RNA_expr_nor %>%
        filter(!value==0)
    AW_quan <- quantile(RNA_expr_nor$value,probs=c(0,0.25,0.5,0.75,1),na.rm=T,names=T)
    first_group <- 
      row.names(RNA_expr_nor)[
        RNA_expr_nor$value>AW_quan['0%'] & RNA_expr_nor$value<=AW_quan['25%']]
    second_group <- 
      row.names(RNA_expr_nor)[
        RNA_expr_nor$value>AW_quan['25%'] & RNA_expr_nor$value<=AW_quan['50%']]
    third_group <- 
      row.names(RNA_expr_nor)[
        RNA_expr_nor$value>AW_quan['50%'] & RNA_expr_nor$value<=AW_quan['75%']]
    forth_group <- 
      row.names(RNA_expr_nor)[
        RNA_expr_nor$value>AW_quan['75%'] & RNA_expr_nor$value<=AW_quan['100%']]
    
    # load annotated peak file
    peak_anno<- 
      read.table(annotated_peak,
                 header = T, stringsAsFactors = T,sep = "\t",comment.char = "",
                 quote="\"") %>%
      filter(ENSEMBL %in% row.names(RNA_expr_nor))
    peak_anno <-
      peak_anno %>%
      mutate(group = case_when(ENSEMBL %in% first_group ~ "first",
                              ENSEMBL %in% second_group ~ "second",
                              ENSEMBL %in% third_group ~ "third",
                              ENSEMBL %in% forth_group ~ "forth"
      ))
    print(table(peak_anno$group))
    total_TagMat <- 
      do.call(rbind,lapply(c("first","second","third","forth"),
                           get_TagMat,peak_file=peak_anno))
    total_TagMat$group <- factor(total_TagMat$group,levels=c("first","second","third","forth")) 
    p <- 
      ggplot(total_TagMat,aes(pos,Peak_count_freq,color=group))+
      geom_line()+
      geom_vline(xintercept = c(80,880),linetype=2)+
      scale_color_manual(values = c("#ff595e","#ffca3a","#1982c4","#8ac926"))+
      scale_x_continuous(name = "",
                         breaks = c(0,80,280,480,680,880,960),
                         labels = c("-1Kb","TSS","25%","50%","75%","TES","+1Kb"))+
      scale_y_continuous(name = "peak count frequency",
                         labels = function(x) format(x, scientific = TRUE))+
      labs(title=sample_name)+
      theme_bw()+
      theme(panel.grid = element_blank(),
            axis.text = element_text(colour = "black",size = 10),
            plot.title = element_text(hjust = 0.5,vjust = 0.5,size=12))
    ggsave(p, filename=paste0("E:\\20220921-WSL-5hmc\\analysis\\corr_RNA_5hmC_5mC\\",sample_name,".pdf"),
           width=13, height=8) 
  }

normalized_count_df <- loadRDS(paste0(data_dir,"//normalized_count_df.rds"))
peak_dir <- "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\"
for(i in colnames(normalized_count_df)){
  for(j in c("H","M")){
    RNA_expr <-
      data.frame(value=normalized_count_df[,i],
                 row.names=row.names(normalized_count_df))
    if(j=="H"){
      plot_RNA_group_Tagcov(RNA_expr, paste0(peak_dir,i,j,"_peaks.narrowPeak_anno"),paste0(i," 5hmC"))
    }
    else{
      plot_RNA_group_Tagcov(RNA_expr, paste0(peak_dir,i,j,"_peaks.narrowPeak_anno"),paste0(i," 5mC"))
    }
  }
}
