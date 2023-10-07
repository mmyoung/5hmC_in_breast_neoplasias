BiocManager::install("rGREAT")
library(rGREAT)
library(dplyr)

pair_ls <- c("ADH_DCIS","UDH_ADH","DCIS_PT")
run_GREAT = function(pair,out_path="E:\\20220921-WSL-5hmc\\analysis\\GREAT_DHMRs_res\\"){
  diff_file = paste0("E:\\20220921-WSL-5hmc\\data\\diff_peak\\",pair,"_H_diff_peak_deseq2.txt")
  compare_fold1 <- 
    read.table(diff_file,
               header = T, row.names = 1) %>%
    filter(Fold>1 & p.value<0.05) %>%
    select(1:3) %>%
    makeGRangesFromDataFrame()
  
  job = submitGreatJob(compare_fold1,species = "hg19")
  tbl = getEnrichmentTables(job)
  write.table(as.data.frame(tbl$`GO Biological Process`),
              file = paste0(out_path,pair,"_fold1_BP.txt"),
              sep="\t",quote = F,row.names = F)
  write.table(as.data.frame(tbl$`GO Molecular Function`),
              file = paste0(out_path,pair,"_fold1_MF.txt"),
              sep="\t",quote = F,row.names = F
  )
  write.table(as.data.frame(tbl$`GO Cellular Component`),
              file = paste0(out_path,pair,"_fold1_CC.txt"),
              sep="\t",quote = F,row.names = F
  )
  
  
  compare_foldM1 <- 
    read.table(diff_file,
               header = T, row.names = 1) %>%
    filter(Fold<(-1) & p.value<0.05) %>%
    select(1:3) %>%
    makeGRangesFromDataFrame()
  
  job = submitGreatJob(compare_foldM1,species = "hg19")
  tbl = getEnrichmentTables(job)
  write.table(as.data.frame(tbl$`GO Biological Process`),
              file = paste0(out_path,pair,"_foldM1_BP.txt"),
              sep="\t",quote = F,row.names = F)
  write.table(as.data.frame(tbl$`GO Molecular Function`),
              file = paste0(out_path,pair,"_foldM1_MF.txt"),
              sep="\t",quote = F,row.names = F
  )
  write.table(as.data.frame(tbl$`GO Cellular Component`),
              file = paste0(out_path,pair,"_foldM1_CC.txt"),
              sep="\t",quote = F,row.names = F
  )
}

lapply(pair_ls,run_GREAT)




## GREAT for continuously changed regions
out_path <- "E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\"
high_peak <-
  read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuously_high_peak.txt",
           header = T, sep="\t", comment.char = "",
           quote = "\"") %>%
  makeGRangesFromDataFrame()
job = submitGreatJob(high_peak,species = "hg19")
tbl = getEnrichmentTables(job)
write.table(as.data.frame(tbl$`GO Biological Process`),
            file = paste0(out_path,"continuous_high_peak_GREAT.txt"),
            sep="\t",quote = F,row.names = F)

## GREAT for continuously low regions
low_peak <-
  read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuously_low_peak.txt",
             header = T, sep="\t", comment.char = "",
             quote = "\"") %>%
  makeGRangesFromDataFrame()
job = submitGreatJob(low_peak,species = "hg19")
tbl = getEnrichmentTables(job)
write.table(as.data.frame(tbl$`GO Biological Process`),
            file = paste0(out_path,"continuous_low_peak_GREAT.txt"),
            sep="\t",quote = F,row.names = F)

## GREAT for stage-specific high regions
peak_df <-
  read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\stage_spec_high_peak.txt",
             header = T, sep="\t", comment.char = "",
             quote = "\"") %>%
  makeGRangesFromDataFrame()
job = submitGreatJob(peak_df,species = "hg19")
tbl = getEnrichmentTables(job)
write.table(as.data.frame(tbl$`GO Biological Process`),
            file = paste0(out_path,"stage_spec_high_GREAT.txt"),
            sep="\t",quote = F,row.names = F)

## GREAT analysis for new continuously changed peaks
library(openxlsx)
library(GenomicRanges)
library(dplyr)
library(rGREAT)
out_path <- "E:\\20220921-WSL-5hmc\\analysis\\GREAT_DHMRs_res\\"

for (i in c(2,3,4)){
  tmp<- 
    read.xlsx("E:\\20220921-WSL-5hmc\\pre-info\\连续分开统计low表.xlsx",
              sheet=i) %>%
    makeGRangesFromDataFrame()
  job = submitGreatJob(tmp,species = "hg19")
  tbl = getEnrichmentTables(job,download_by='tsv')
  write.table(as.data.frame(tbl$`GO Biological Process`),
              file = paste0(out_path,"continuous_low_peak_GREAT_",i,".txt"),
              sep="\t",quote = F,row.names = F)
}

for (i in c(1,2,3,4)){
  tmp<- 
    read.xlsx("E:\\20220921-WSL-5hmc\\pre-info\\连续分开统计high表.xlsx",
              sheet=i) %>%
    makeGRangesFromDataFrame()
  job = submitGreatJob(tmp,species = "hg19")
  tbl = getEnrichmentTables(job,download_by='tsv')
  write.table(as.data.frame(tbl$`GO Biological Process`),
              file = paste0(out_path,"continuous_high_peak_GREAT_",i,".txt"),
              sep="\t",quote = F,row.names = F)
}

