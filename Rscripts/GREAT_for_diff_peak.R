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





