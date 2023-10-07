## 2023-01-17
## calculate the correlation between 5mC and 5hmC
library(GenomicRanges)
library(ChIPseeker)
## A2 test


## read 5hmC data
hmC_mC_corr_cal <-
function(sample_name){
  H_file <- 
    paste0("E:\\20220921-WSL-5hmc\\data\\sample_peak\\",sample_name,"H_peaks.narrowPeak")
  M_file <-
    paste0("E:\\20220921-WSL-5hmc\\data\\sample_peak\\",sample_name,"M_peaks.narrowPeak")
  if(file.exists(H_file) & file.exists(M_file)){
    H_raw <-
      readPeakFile(H_file,
                   header=F)
    M_raw<- 
      readPeakFile(M_file,
                   header = F)
    ## pearson test between enrichment scores from 5mC and 5hmC
    OL_HM <- as.data.frame(findOverlaps(H_raw,M_raw))
    cor_test_res <- 
      cor.test(as.data.frame(H_raw[OL_HM[,"queryHits"],"V7"])[,"V7"],
               as.data.frame(M_raw[OL_HM[,"subjectHits"],"V7"])[,"V7"])
    return(data.frame(P=cor_test_res$p.value,Corr=cor_test_res$estimate)) 
  }
  else{
    return(data.frame(P=NA,Corr=NA))
  }
}

library(stringr)
all_sample_name <-
  unique(str_replace_all(str_split(list.files(path = "E:\\20220921-WSL-5hmc\\data\\sample_peak",
           pattern = "[^P]*\\.narrowPeak"),"_",simplify = T)[,1],"H|M",""))
all_sample_name
corr_table <- do.call(rbind,lapply(all_sample_name,hmC_mC_corr_cal))
write.table(corr_table,
            "E:\\20220921-WSL-5hmc\\analysis\\")


