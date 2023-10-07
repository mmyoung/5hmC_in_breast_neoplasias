library(openxlsx)
library(ggplot2)
TF_tf <-
  read.xlsx("E:\\20220921-WSL-5hmc\\data\\TET2bind DhMRs motif.xlsx",sheet = 2)
TF_tf$Pair <- factor(TF_tf$Pair,
                     levels = c("UDHvsADH","DCISvsADH","IDCvsDCIS"),
                     labels = c("Stage I","Stage II","Stage III"))
TF_tf$TFs <- factor(TF_tf$TFs,
                     levels = c("Fra1","Fosl2","JunB","Fos","FOXA1","ERE"))

(ggplot(TF_tf,aes(x=Pair,y=TFs,fill=-log10(Pvalue)))+
  geom_tile()+
  labs(title = "tet2-binding DhMRs")+
  scale_x_discrete(expand = c(0,0),name="")+
  scale_y_discrete(expand = c(0,0),name="")+
  scale_fill_gradientn(colours = c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))+
  theme(panel.background = element_rect(fill = "grey",colour = "grey"),
        panel.grid = element_line(color = "grey"),
        axis.ticks = element_blank(),
        axis.text = element_text(colour = "black",size = 11),
        plot.title = element_text(hjust = 0.5,vjust = 0.5)
        )) |>
  ggsave(filename = "E:\\20220921-WSL-5hmc\\analysis\\GSE153251_tet2_ChIPSeq\\motif_heatmap.pdf",
         width = 5,height = 4)
