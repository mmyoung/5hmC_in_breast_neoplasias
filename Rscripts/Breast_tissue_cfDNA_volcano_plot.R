## make volcano plots for DhMRs in breast tumor and cfDNA

library(ggrastr)
library(ggplot2)
library(dplyr)
library(stringr)

make_plot <- 
  function(data,data_prefix){
    data <- read.table(data,
                       header = T,stringsAsFactors = F,
                       comment.char="",sep="\t",quote="\"")
    pair_ls <- str_split(data_prefix,pattern = "_",simplify = T)
    
    up_num <- data %>%
      filter(p.value <0.05 & Fold < -1) %>%
      nrow()
    down_num <- data %>%
      filter(p.value <0.05 & Fold > 1) %>%
      nrow()
    
    ( ggplot(data=data, aes(x=-Fold, y =-log2(p.value)))+
        geom_point_rast(data=subset(data,abs(data$Fold)<1 | data$p.value>0.05),color="#828C8C",alpha=0.7)+
        ## down
        geom_point_rast(data=subset(data,data$p.value <0.05 & data$Fold > 1),color="#1f70a9",alpha=0.7) +
        ## up
        geom_point_rast(data=subset(data,data$p.value<0.05 & data$Fold < -1),color="#c22f2f",alpha=0.7) +
        #geom_point_rast(x=-5,y=40,color="#c22f2f",size=1.5)+
        #geom_point_rast(x=-5,y=38,color="#1f70a9",size=1.5)+
        annotate("text",x=3.5,y=40,label=paste0("UP ",up_num),size=3)+
        annotate("text",x=-3.5,y=40,label=paste0("DOWN ",down_num),size=3)+
        coord_cartesian(xlim = c(-5,5),
                        clip = 'off') +
        #scale_color_manual(values = c("down"="#1f70a9",
        #                              "up"="#c22f2f",
        #                              "nosig"="#828C8C"))+
        geom_hline(yintercept = -log2(0.05),lty=4,lwd=0.6,alpha=0.8,color="#828C8C")+
        geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8,color="#828C8C")+
        theme_bw()+
        theme(legend.position='none',
              plot.title = element_text(hjust = 0.5,vjust = 0.5),
              #panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),   
              axis.line = element_line(colour = 'black'),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12)
        )+
        labs(x='log2 (Foldchange)',y='-log2 (p-value)',
             title=paste0(pair_ls[2]," vs ",pair_ls[1]))) %>%
      ggsave(filename = paste0("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\",data_prefix,"_diffPeak_volcano.pdf"),
             width=4,height = 5)
  }

make_plot(data="E:\\20220921-WSL-5hmc\\data\\diff_peak_new\\UDH_PT_H_diff_peak_deseq2.txt",
          data_prefix = "UDH_PT_H")
make_plot(data="E:\\20220921-WSL-5hmc\\data\\cfDNA\\breast_healthy_diff_peak_deseq2.txt",
          data_prefix = "N_cfDNA T_H")


## count the number of DMRs in each pair
peak_region_count <- data.frame()
diff_res <- "E:\\20220921-WSL-5hmc\\data\\cfDNA\\breast_healthy_diff_peak_deseq2.txt"
base_name <- paste(str_split(basename(diff_res),"_",simplify = T)[1:3],collapse = "-")
diff_res_df <- read.table(diff_res,header = T,sep = "\t",comment.char = "",quote = "\"")
tmp_df <-
  diff_res_df %>%
  mutate(diff_status=ifelse(p.value<0.05 & Fold>1,"P<0.05&Fold>1",
                            ifelse(p.value<0.05 & Fold<(-1), "P<0.05&Fold<-1","NS")))%>%
  filter(!diff_status=="NS") %>%
  group_by(diff_status) %>%
  summarise(n()) %>%
  as.data.frame() %>%
  mutate(comparison=!!base_name)
peak_region_count <- rbind(peak_region_count,tmp_df)

write.table(peak_region_count,
            file = "E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno\\cfDNA_DhMR_counts.txt",
            quote = F,row.names = F,sep = "\t")

peak_region_count <- data.frame()
diff_res <- "E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno\\cfDNA_diff_5hmC_ChIPseeker_anno.txt"
base_name <- paste(str_split(basename(diff_res),"_",simplify = T)[1:3],collapse = "-")
diff_res_df <- read.table(diff_res,header = T,sep = "\t",comment.char = "",quote = "\"")
tmp_df <-
  diff_res_df %>%
  mutate(diff_status=ifelse(p.value<0.05 & Fold>1,"P<0.05&Fold>1",
                            ifelse(p.value<0.05 & Fold<(-1), "P<0.05&Fold<-1","NS")))%>%
  filter(!diff_status=="NS") %>%
  group_by(diff_status) %>%
  filter(!ENSEMBL=="NA") %>% 
  filter(!duplicated(ENSEMBL)) %>%
  summarise(n()) %>%
  as.data.frame() %>%
  mutate(comparison=!!base_name)
peak_region_count <- rbind(peak_region_count,tmp_df)

write.table(peak_region_count,
            file = "E:\\20220921-WSL-5hmc\\analysis\\cfDNA_peak_anno\\cfDNA_DhMG_counts.txt",
            quote = F,row.names = F,sep = "\t")
