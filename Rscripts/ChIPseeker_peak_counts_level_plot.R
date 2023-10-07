peak_region_counts <- 
read.table("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\01_peak_region_counts.txt",
           header = T, stringsAsFactors = F,
           quote = "\"",sep = "\t")
head(peak_region_counts)
hmC_peak_region <- 
  peak_region_counts[grep("H",peak_region_counts$comparison),]
head(hmC_peak_region)
library(dplyr)
hmC_peak_region <-
  hmC_peak_region %>%
  mutate(stage=case_when(comparison %in% c("A10H","A2H","AWH") ~ "ADH",
                            comparison %in% c("D10H","D20H","D3H") ~ "DCIS",
                            comparison %in% c("JH","U3H","U5H") ~ "UDH",
                            comparison %in% c("PT3H","PT5H","PT6H","PT7H") ~ "PT"
  ))
head(hmC_peak_region)
hmC_peak_region$comparison <- as.character(hmC_peak_region$comparison)
hmC_peak_region$stage <- factor(hmC_peak_region$stage,
                                levels = c("UDH","ADH","DCIS","PT"))

library(ggplot2)
library(ggsci)
p <-
  ggplot(hmC_peak_region,aes(x=comparison,y=n..,fill=group))+
  geom_bar(stat = "identity",position = "fill")+
  facet_grid(.~stage,scales = "free")+
  scale_fill_nejm()+
  theme_bw()+
  labs(x="sample",y="peak proportion")
ggsave(p,filename = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\peak_region_counts.pdf",
       width = 10,height = 8)


peak_region_level <- 
  read.table("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\01_peak_region_avg.txt",
             header = T, stringsAsFactors = F,
             quote = "\"",sep = "\t")

hmC_peak_level <- 
  peak_region_level[grep("H",peak_region_level$sample),]
head(hmC_peak_level)
library(dplyr)
hmC_peak_level <-
  hmC_peak_level %>%
  mutate(stage=case_when(sample %in% c("A10H","A2H","AWH") ~ "ADH",
                         sample %in% c("D10H","D20H","D3H") ~ "DCIS",
                         sample %in% c("JH","U3H","U5H") ~ "UDH",
                         sample %in% c("PT3H","PT5H","PT6H","PT7H") ~ "PT"
  ))
hmC_peak_level$stage <- factor(hmC_peak_level$stage,
                                levels = c("UDH","ADH","DCIS","PT"))
p<-
  ggplot(hmC_peak_level,aes(x=sample,y=name,color=group))+
  geom_point(size=2)+
  #geom_bar(stat = "identity",position = "dodge")+
  facet_grid(.~stage,scales = "free")+
  scale_color_nejm()+
  theme_bw()+
  labs(x="sample",y="Fold Enrichment")

ggsave(p,filename = "E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_peak_anno\\peak_region_level.pdf",
       width = 10,height = 8)
