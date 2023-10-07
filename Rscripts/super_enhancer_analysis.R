
## enhancer related analyses
setwd("E:\\20220921-WSL-5hmc\\analysis\\enhancer_analysis")

## reproduce the enhancer ranking plot
enhancer_df <-
  read.table("MCF7_enhancer_AllEnhancers.table.txt",header = T)
enhancer_df$rank <- seq.int(from=nrow(enhancer_df),
                            to=1)

library(ggplot2)

enhancer_df$rank[enhancer_df$MCF7_H3K27ac_align_rep1_sorted.bam==19377.435]

p <-
ggplot(enhancer_df,aes(rank,MCF7_H3K27ac_align_rep1_sorted.bam))+
  geom_point()+
  geom_vline(xintercept = 39936,linetype="dotted")+
  geom_text(x=10000,y=400000,
            label="Super-enhancer identified: 1135")+
  labs(y="H3K27ac signal",x="Ranked Enhancers")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"))
ggsave(p,
       filename = "E:\\20220921-WSL-5hmc\\analysis\\enhancer_analysis\\SE_identification.pdf",
       width = 5,height = 4)


## compare the fold change between peaks located in SE and non-SE region
DCIS_PT_foldM1_df  <- read.table("DCIS_PT_H_foldM1_SE_anno.txt",
           header = F,stringsAsFactors = F)
UDH_ADH_fold1_df  <- read.table("UDH_ADH_H_fold1_SE_anno.txt",
                                 header = F,stringsAsFactors = F)
table(DCIS_PT_foldM1_df$V19)
#0    1 
#8383 2265 
table(UDH_ADH_fold1_df$V19)
#0    1 
#5786 1589

t.test(abs(DCIS_PT_foldM1_df$V9[which(DCIS_PT_foldM1_df$V19==0)]),
       abs(DCIS_PT_foldM1_df$V9[which(DCIS_PT_foldM1_df$V19==1)]))
#t = 3.3475, df = 3645.6, p-value = 0.0008236

p <-
  ggplot(DCIS_PT_foldM1_df,aes(as.factor(V19),abs(V9),fill=as.factor(V19)))+
    geom_boxplot()+
    geom_text(x=1.5,y=3.7,label="p-value = 0.0008236",size=4)+
    labs(x="",y="abs(log2FoldChange)",title = "hyper-5hmC (DCIS->PT)")+
    scale_x_discrete(labels=c("non-SE","SE"))+
    scale_fill_manual(values=c("#a8dadc","#457b9d"))+
    scale_y_continuous(limits = c(1,4))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=10,colour = "black"),
          plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 12),
          legend.position = "NA")
ggsave(p,
       filename = "E:\\20220921-WSL-5hmc\\analysis\\enhancer_analysis\\DCIS_PT_foldM1_SEvsnonSE_fold_boxplot.pdf",
       width = 5,height = 4)


ggplot(UDH_ADH_fold1_df,aes(as.factor(V19),abs(V9)))+
  geom_boxplot()

t.test(abs(UDH_ADH_fold1_df$V9[which(UDH_ADH_fold1_df$V19==0)]),
       abs(UDH_ADH_fold1_df$V9[which(UDH_ADH_fold1_df$V19==1)]))
#t = 2.7959, df = 2620.2, p-value = 0.005213

p <-
  ggplot(UDH_ADH_fold1_df,aes(as.factor(V19),abs(V9),fill=as.factor(V19)))+
  geom_boxplot()+
  geom_text(x=1.5,y=3.7,label="p-value = 0.005213",size=4)+
  labs(x="",y="abs(log2FoldChange)",title = "hypo-5hmC (UDH->ADH)")+
  scale_x_discrete(labels=c("non-SE","SE"))+
  scale_fill_manual(values=c("#a8dadc","#457b9d"))+
  scale_y_continuous(limits = c(1,4))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10,colour = "black"),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 12),
        legend.position = "NA")
ggsave(p,
       filename = "E:\\20220921-WSL-5hmc\\analysis\\enhancer_analysis\\UDH_ADH_fold1_SEvsnonSE_fold_boxplot.pdf",
       width = 5,height = 4)

## 
DCIS_PT_foldM1_SE <-
  gsub(" ","",apply(DCIS_PT_foldM1_df[which(DCIS_PT_foldM1_df$V19==1),c("V12","V13","V14")],1,paste,collapse = "_"))
UDH_ADH_fold1_SE <-
  gsub(" ","",apply(UDH_ADH_fold1_df[which(UDH_ADH_fold1_df$V19==1),c("V12","V13","V14")],1,paste,collapse = "_"))
common_change_SE <- intersect(DCIS_PT_foldM1_SE,UDH_ADH_fold1_SE)

MCF7_SE_anno <- 
  read.table("MCF7_SE_anno.txt",sep = "\t",
           comment.char = "",
           quote = "\"",
           header = T)
row.names(MCF7_SE_anno) <- 
  paste(MCF7_SE_anno$Chr,MCF7_SE_anno$Start-1,MCF7_SE_anno$End,sep = "_")
common_SE_gene <- MCF7_SE_anno[common_change_SE,"Gene.Name"]
library(clusterProfiler)
library(org.Hs.eg.db)
options(clusterProfiler.download.method = "wget")
ego <- 
  enrichGO(common_SE_gene,
           keyType = "SYMBOL",
           OrgDb = org.Hs.eg.db,
           ont = "BP")
ekegg <- 
enrichKEGG(common_SE_gene,organism = "hsa") ## cannot download KEGG, implemented in DAVID
write.table(common_SE_gene,file = "common_SE_gene.txt",
            quote=F,row.names = F)


