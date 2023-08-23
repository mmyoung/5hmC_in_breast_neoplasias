
setwd("E:\\20220921-WSL-5hmc\\analysis\\enhancer_analysis")

DCIS_PT_foldM1 <-
  read.table("DCIS_PT_H_foldM1_enhancer.txt",
           sep = "\t", stringsAsFactors = F,header = F)
DCIS_PT_foldM1$V19 <- factor(DCIS_PT_foldM1$V19,
                             levels = c(0,1),
                             labels = c("non-SE","SE"))
table(DCIS_PT_foldM1$V19)
head(DCIS_PT_foldM1)
library(ggplot2)
ggplot(DCIS_PT_foldM1,aes(x=V19,y=V9))+
  geom_boxplot()

UDH_ADH_fold1 <-
  read.table("UDH_ADH_H_fold1_enhancer.txt",
             sep = "\t", stringsAsFactors = F,header = F)
UDH_ADH_fold1$V19 <- factor(UDH_ADH_fold1$V19,
                             levels = c(0,1),
                             labels = c("non-SE","SE"))
head(UDH_ADH_fold1)
library(ggplot2)
ggplot(UDH_ADH_fold1,aes(x=V19,y=V9))+
  geom_boxplot()


## GO annotation for the enhancer targets with altered 5hmC 
UDH_ADH_fold1_enhancer <-
  read.table("UDH_ADH_H_fold1_enhancer.txt",
             header = F, stringsAsFactors = F,
             comment.char = "",
             sep = "\t",
             quote = "\"")
DCIS_PT_foldM1_enhancer <-
  read.table("DCIS_PT_H_foldM1_enhancer.txt",
             header = F, stringsAsFactors = F,
             comment.char = "",
             sep = "\t",
             quote = "\"")
library(clusterProfiler)
library(org.Hs.eg.db)
ego_1<- 
enrichGO(UDH_ADH_fold1_enhancer$V25,
         keyType = "ENSEMBL",
         OrgDb = org.Hs.eg.db,
         ont = "BP",
         readable = T)

ego_2<- 
  enrichGO(DCIS_PT_foldM1_enhancer$V25,
           keyType = "ENSEMBL",
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = T)

ekegg_1 <- enrichKEGG(UDH_ADH_fold1_enhancer$V25)
ekegg_2 <- enrichKEGG(DCIS_PT_foldM1_enhancer$V25)

library(DOSE)
dotplot(ego_2)
