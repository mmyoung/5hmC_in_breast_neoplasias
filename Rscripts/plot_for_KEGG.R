library(ggplot2)
library(stringr)

library(RColorBrewer)
library(scales)
data_dir <- "E:\\20220921-WSL-5hmc\\analysis\\KEGG\\"
data_dir <- "/Users/linyang/.mounty/Seagate Expansion Drive/电脑备份/E盘/20220921-WSL-5hmc/analysis/KEGG/"
setwd(data_dir)
#data <- read.table(paste0(data_dir,"continuous_up.txt"),
data <- read.table(paste0(data_dir,"continuous_down.txt"),
           header = T,stringsAsFactors = F,sep = "\t")

data$Term <- str_split(data$Term,":",simplify = T)[,2]
data <- data[order(data$`X.`),]
data$Term <- factor(data$Term,levels = data$Term)
colnames(data)

(ggplot(data,aes(x=Term,y=`X.`,color=-log10(PValue),size=`X.`))+
  #geom_bar(stat = "identity",width = 0.6,position = position_dodge(0))+
  geom_point()+
  labs(y = "Gene Ratio (%)", x = "KEGG Term",fill="-log10(PValue)")+
  theme_bw()+
  scale_y_continuous(expand = c(0,0.5))+
  coord_flip()+
    scale_color_gradient2(
      #high = muted("#c22f2f"),
      high = "#c22f2f",
      mid = "grey",
      low = "#1f70a9",
      #low = muted("#1f70a9"),
      midpoint = 6)+
      #midpoint = 4.5)+
    labs(size="Gene Ratio (%)")+
  theme(panel.grid =element_blank())
  ) %>%
  ggsave(filename = "continuous_down.pdf",width = 7,height = 5)
  #ggsave(filename = "continuous_up.pdf",width = 7,height = 5)


## 

up_data <- read.table(paste0(data_dir,"cfDNA_up.txt"),
#up_data <- read.table(paste0(data_dir,"tissue_up.txt"),
                  header = T,stringsAsFactors = F,sep = "\t")
down_data <- read.table(paste0(data_dir,"cfDNA_down.txt"),
#down_data <- read.table(paste0(data_dir,"tissue_down.txt"),
                   header = T,stringsAsFactors = F,sep = "\t")

data_for_plot <- rbind(up_data,down_data)
data_for_plot$status <- 
  c(rep("up",nrow(up_data)),rep("down",nrow(down_data)))

data_for_plot <- data_for_plot[order(data_for_plot$status,data_for_plot$Count),]
data_for_plot$Term <- str_split(data_for_plot$Term,":",simplify = T)[,2]

data$Term <- factor(data$Term,levels = data$Term)
data_for_plot$status <- factor(data_for_plot$status,levels = c("up","down"))
(ggplot(data_for_plot,aes(x=status,y=Term,color=-log10(PValue),size=Count))+
  geom_point()+
  scale_color_gradient2(
    high = "#c22f2f",
    mid = "grey",
    low = "#1f70a9",
    midpoint = 6
    #midpoint = 9
    )+
  labs(title = "cfDNA DhMG",x="",y="KEGG Term")+
  #labs(title = "tissue DhMG",x="",y="KEGG Term")+
  guides(shape="none")+
  theme_bw()) %>%
  ggsave(filename = "cfDNA_KEGG.pdf",width = 5,height = 6)
  #ggsave(filename = "tissue_KEGG.pdf",width = 5,height = 6)


## Plot for DEG & DhMGs

# I_data <- read.table(paste0(data_dir,"DEG-DhMGs phase1 DAVID-KEGG.txt"),
#                       header = T,stringsAsFactors = F,sep = "\t")
# II_data <- read.table(paste0(data_dir,"DEG-DhMGs phase2 DAVID-KEGG.txt"),
#                         header = T,stringsAsFactors = F,sep = "\t")
# III_data <- read.table(paste0(data_dir,"DEG-DhMGs phase3 DAVID-KEGG.txt"),
#                       header = T,stringsAsFactors = F,sep = "\t")
# data_for_plot <- rbind(I_data,II_data,III_data)
# 
# data_for_plot$status <-
#   c(rep("I",nrow(I_data)),rep("II",nrow(II_data)),rep("III",nrow(III_data)))

I_data <- read.table(paste0(data_dir,"DEG-DMGs phase1 DAVID-KEGG.txt"),
                     header = T,stringsAsFactors = F,sep = "\t")
II_data <- read.table(paste0(data_dir,"DEG-DMGs phase2 DAVID-KEGG.txt"),
                      header = T,stringsAsFactors = F,sep = "\t")

data_for_plot <- rbind(I_data,II_data)
data_for_plot$status <-
  c(rep("I",nrow(I_data)),rep("II",nrow(II_data)))

data_for_plot <- data_for_plot[order(data_for_plot$status,data_for_plot$Count),]
data_for_plot$Term <- str_split(data_for_plot$Term,":",simplify = T)[,2]

(ggplot(data_for_plot,aes(x=status,y=Term,color=-log10(PValue),size=Count))+
    geom_point()+
    scale_color_gradient2(
      high = "#c22f2f",
      mid = "grey",
      low = "#1f70a9",
      #midpoint = 3.5
      midpoint = 4
      )+
    guides(shape="none")+
    theme_bw()) %>%
  ggsave(filename = "DEG_DMGs_KEGG.pdf",width = 7,height = 6)
  #ggsave(filename = "DEG_DhMGs_KEGG.pdf",width = 7,height = 6)

## make plot for DhMGs and DMGs

## organizing the data
data_dir <- "E:\\20220921-WSL-5hmc\\analysis\\KEGG\\S3D\\"
files_ls <- list.files(path = data_dir,pattern = "*.txt")

test <-
do.call(rbind,lapply(files_ls,
       function(x){
         data.frame(read.table(paste0(data_dir,x),
                    header = T,stringsAsFactors = F,sep = "\t"),
                    x_axis=str_replace_all(x,".txt",""),
                    phase=str_split(str_replace_all(x,".txt",""),pattern = "_",simplify = T)[1],
                    change=str_split(str_replace_all(x,".txt",""),pattern = "_",simplify = T)[2])
       }
       ))
       
test$x_axis <- 
  factor(test$x_axis,levels = c("Phase2_Hu","Phase3_Hu","Phase1_Hd","Phase2_Hd","Phase3_Hd"))
test$phase <- factor(test$phase,
                     levels = c("Phase1","Phase2","Phase3"),
                     labels = c("I","II","III"))
test$change <- factor(test$change,
                      levels = c("Hu","Hd"),
                      labels = c("UP","DOWN"))
library(ggplot2)  
(ggplot(test,aes(x=phase,y=Description,
                 color=-log10(pvalue),size=Count))+
    geom_point()+
    facet_grid(.~change,scales = "free_x")+
    scale_color_gradient2(
      high = "#c22f2f",
      mid = "grey",
      low = "#1f70a9",
      #midpoint = 3.5
      midpoint = 12.5
    )+
    #guides(shape="none")+
    theme_bw()) %>%
  ggsave(filename = paste0(data_dir,"DhMG_KEGG.pdf"),width = 7,height = 6)

data_dir <- "E:\\20220921-WSL-5hmc\\analysis\\KEGG\\S4D\\"
files_ls <- list.files(path = data_dir,pattern = "*.txt")

test <-
  do.call(rbind,lapply(files_ls,
                       function(x){
                         data.frame(read.table(paste0(data_dir,x),
                                               header = T,stringsAsFactors = F,sep = "\t"),
                                    x_axis=str_replace_all(x,".txt",""),
                                    phase=str_split(str_replace_all(x,".txt",""),pattern = "_",simplify = T)[1],
                                    change=str_split(str_replace_all(x,".txt",""),pattern = "_",simplify = T)[2])
                       }
  ))

test$phase <- factor(test$phase,
                     levels = c("Phase1","Phase2"),
                     labels = c("I","II"))
test$change <- factor(test$change,
                      levels = c("Mu","Md"),
                      labels = c("UP","DOWN"))
library(ggplot2)  
(ggplot(test,aes(x=phase,y=Description,
                 color=-log10(pvalue),size=Count))+
    geom_point()+
    facet_grid(.~change,scales = "free_x")+
    scale_color_gradient2(
      high = "#c22f2f",
      mid = "grey",
      low = "#1f70a9",
      #midpoint = 3.5
      midpoint = 6
    )+
    #guides(shape="none")+
    theme_bw()) %>%
  ggsave(filename = paste0(data_dir,"DMG_KEGG.pdf"),width = 7,height = 6)
