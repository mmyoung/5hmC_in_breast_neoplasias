
## Date: 2023-05-11
## Auther: Lin Yang

## This was ran on the server in ChIPseeker conda environment
out_path <-
  "/data/nym/20220921-WSL-5hmc/analysis/MCF7_TFs_correlation/DhMRs_tet2_TF_binding_anno"
DhMRs_TF_annotation <-
  function(stage_DhMR){
    
    ## load DhMR files with tet2 binding annotation
    tet2_DhMR_OL <- 
      paste0("/data/nym/20220921-WSL-5hmc/analysis/GSE153251_tet2_ChIPseq/DhMRs_tet2_OL/",stage_DhMR,"_tet2_ol.txt")
    tet2_binding_sites <-
      read.table(tet2_UA_DhMR, sep = "\t", header = T, quote="\"",comment.char="") %>%
      makeGRangesFromDataFrame(keep.extra.columns=T)
    tmp_peak_df <- as.data.frame(tet2_binding_sites)
    ## add the annotation of all TFs into the original file
    for(TF in c("ESR1","FOS","FOSL2","FOXA1","FOXM1","GATA3","JUNB")){
      TF_ChIPseq_bed <- 
        paste0("/data/nym/20220921-WSL-5hmc/data/ENCODE_MCF7_ChIPseq/",TF,"_ChIPseq_MCF7.bed")
      TF_binding_sites <-
        read.table(TF_ChIPseq_bed,sep="\t", header = F, quote="\"",comment.char="") %>%
        makeGRangesFromDataFrame(seqnames.field = "V1",
                                 start.field = "V2",
                                 end.field = "V3",
                                 keep.extra.columns=T)
      overlap_index <- 
        as.data.frame(GenomicRanges::findOverlaps(tet2_binding_sites,
                                                  TF_binding_sites))$queryHits
      tmp_peak_df[,TF] <- "no"
      tmp_peak_df[overlap_index,TF] <- "yes"
    }
    
    write.table(tmp_peak_df,
                file = paste0(out_path,"/",stage_DhMR,"_TF_anno.txt"),
                quote = F, row.names = F,sep = "\t")
}

lapply(c("UDH_ADH_H","UDH_ADH_M","ADH_DCIS_H","ADH_DCIS_M","DCIS_PT_H"),
       DhMRs_TF_annotation)