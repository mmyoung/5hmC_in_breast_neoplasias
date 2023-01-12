library(dplyr)
setwd("/data/nym/20220921-WSL-5hmc/analysis/CGI_region_distri")
res_df <- data.frame()
for (i in c("A10","A2","ADH_","AW","D10","D20","D3","DCIS_","J","U3","U5","UDH_")){
	raw_file <- read.table(paste0(i,"5hmC_mC_in_CGI_region.txt"),header=T,stringsAsFactors=F)
	for (region_H in c("shelf","shore","CGI","sea")){
		for (region_M in c("shelf","shore","CGI","sea")){
			H_region_level<- raw_file[which(raw_file$region == region_H),c("seqname","start","end","hmC_level")]
			M_region_level<- raw_file[which(raw_file$region == region_M),c("seqname","start","end","mC_level")]
			if(dim(H_region_level)[1]>3 & dim(M_region_level)[1]>3){
				merged_df <- merge(H_region_level, M_region_level, by=c("seqname","start","end"))
				print(paste0("processing correlation analysis between 5hmC in ",region_H," and 5mC in ",region_M," in sample ",i))
				res_corr <- cor.test(merged_df$hmC_level,merged_df$mC_level)$estimate
				res_p <- cor.test(merged_df$hmC_level,merged_df$mC_level)$p.value
				res_df <- rbind(res_df,data.frame(sample=i, hmC_region=region_H, mC_region=region_M, corr=res_corr, p=res_p))
			}
		}
	}
}
write.table(res_df,file="hmC_mC_corr_each_sample.txt",sep="\t",quote=F,row.names=F)



