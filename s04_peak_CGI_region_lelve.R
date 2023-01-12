
library(dplyr)
setwd("/data/nym/20220921-WSL-5hmc/analysis/CGI_region_distri")
for (i in c("A10","A2","ADH_","AW","D10","D20","D3","DCIS_","J","U3","U5","UDH_")){

	hm_file <- paste0(i,"H_peaks.narrowPeak.cgi_distance")
	m_file <- paste0(i,"M_peaks.narrowPeak.cgi_distance")
	sample_5hmC <- read.table(hm_file,header=F,stringsAsFactors=F,sep="\t")
	sample_5mC <- read.table(m_file,header=F,stringsAsFactors=F,sep="\t")
	tmp <- sample_5hmC %>% filter(grepl("^chr",V11)) %>% group_by(V11,V12,V13,V18) %>% summarise(mean_level=mean(V7))
	tmp_M <- sample_5mC %>% filter(grepl("^chr",V11)) %>% group_by(V11,V12,V13,V18) %>% summarise(mean_level=mean(V7))
	res <- merge(tmp,tmp_M,by=c("V11","V12","V13","V18"))
	colnames(res) <- c("seqname","start","end","region","hmC_level","mC_level")
	write.table(res,paste0(i,"5hmC_mC_in_CGI_region.txt"),row.names=F,quote=F,sep="\t")
}


