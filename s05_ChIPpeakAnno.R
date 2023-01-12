

library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
##diff_peak_file <- list.files("/data/nym/20220921-WSL-5hmc/data/diff_peak",full.names=TRUE)[grep("up|down",list.files("/data/nym/20220921-WSL-5hmc/data/diff_peak"))]
diff_peak_file <- list.files("/data/nym/20220921-WSL-5hmc/data/sample_peak",full.names=TRUE)[grep("narrowPeak",list.files("/data/nym/20220921-WSL-5hmc/data/sample_peak"))]
res_df  <- data.frame()
for (peak_path in diff_peak_file){
	peak_basename <- basename(peak_path)
	peak_basename <- str_split(peak_basename,"_",simplify=TRUE)[1]
	peak_file <- read.table(peak_path, sep="\t",stringsAsFactors=F,header=F)
	colnames(peak_file) <- c("seqnames","start","end","length","strand")
	gr1 <- toGRanges(peak_file[,1:5])
	aCR<-assignChromosomeRegion(gr1, nucleotideLevel=FALSE,
                           precedence=c("Promoters", "immediateDownstream",
                                         "fiveUTRs", "threeUTRs",
                                         "Exons", "Introns"),
                           TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
	print(paste0(peak_basename," annotated region: ", length(aCR)))
	res_df <- rbind(res_df,data.frame(sample=peak_basename, region=names(aCR$percentage),percentage=as.numeric(aCR$percentage)))
}
#write.table(res_df,file="/data/nym/20220921-WSL-5hmc/analysis/ChIPpeakAnno_res/diff_peak_genomicRegion_Freq.txt",quote=F,sep="\t",row.names=F)
write.table(res_df,file="/data/nym/20220921-WSL-5hmc/analysis/ChIPpeakAnno_res/sample_peak_genomicRegion_Freq.txt",quote=F,sep="\t",row.names=F)
