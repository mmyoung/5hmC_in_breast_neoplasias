## Correlation between DEGs and DhMRs
## stage UDH->ADH
library(dplyr)
## load RNA-seq data
make_pointplot <- 
  function(RNAseq_prefix,DhMR_prefix,comparison,region){
    RNAseq<- 
      read.table(paste0("E:\\20220921-WSL-5hmc\\analysis\\RNAseq_2014paper_DEGs\\",RNAseq_prefix,"_RNAseq_diff_gene.txt"),
      header = T)

    RNAseq$diff_sign <- 
      ifelse(abs(RNAseq$logFC)>1 & RNAseq$P.Value<0.05,"diff","no")
    
    DhMRs<-
    #read.table(paste0("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\",DhMR_prefix,"_H_diff_peak_deseq2.txt_anno"),
      read.table(paste0("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\",DhMR_prefix,"_M_diff_peak_deseq2.txt_anno"),
                stringsAsFactors = F, header = T,
               sep = "\t",comment.char = "",quote = "\"") 
    DhMRs$diff_sign <- 
      ifelse(abs(DhMRs$Fold)>1 & DhMRs$p.value<0.05,"diff","no")
    
    DhMRs_promoter <- 
      DhMRs[grep(paste(region,collapse = "|"),DhMRs$annotation),c("Fold","SYMBOL","diff_sign")]
    table(DhMRs_promoter$diff_sign)
    
    data_for_plot <-
    merge(RNAseq[,c("Row.names","logFC","diff_sign")],DhMRs_promoter,
          by.x = "Row.names",by.y = "SYMBOL")
    library(ggplot2)
    p <- 
    ggplot(data_for_plot,aes(-Fold,logFC))+
      geom_point(color="grey",alpha=0.7)+
      geom_point(data=subset(data_for_plot,abs(Fold)>1),aes(-Fold,logFC),color="blue",alpha=0.7)+
      geom_point(data=subset(data_for_plot,diff_sign.x=="diff" & diff_sign.y=="diff"),aes(-Fold,logFC),color="red")+
      scale_y_continuous(limits = c(-20,20),name = "log2FC (RNA expression)")+
      #scale_x_continuous(limits = c(-4,4),name = "log2FC (5hmC)")+
      scale_x_continuous(limits = c(-4,4),name = "log2FC (5mC)")+
      geom_hline(yintercept = 0,linetype="dashed")+
      labs(title=comparison)+
      theme_bw()+
      theme(panel.grid = element_blank(),
            axis.text = element_text(size = 10,colour = "black"),
            axis.title = element_text(size = 12,colour = "black"),
            plot.title = element_text(hjust = 0.5,vjust = 0.5))
    ggsave(p, 
           filename = paste0("E:\\20220921-WSL-5hmc\\analysis\\RNAseq_2014paper_DEGs\\",comparison,"_",paste(region,collapse=""),"_DMRs.pdf"),
           height = 4,width = 5)
    }

make_pointplot(RNAseq_prefix = "EN_normal",
               DhMR_prefix = "UDH_ADH",
               comparison = "ADHvsUDH",
               region = "Promoter")
make_pointplot(RNAseq_prefix = "DCIS_EN",
               DhMR_prefix = "ADH_DCIS",
               comparison = "DCISvsADH",
               region = "Promoter")
make_pointplot(RNAseq_prefix = "IDC_DCIS",
               DhMR_prefix = "DCIS_PT",
               comparison = "IDCvsDCIS",
               region = "Promoter")

make_pointplot(RNAseq_prefix = "EN_normal",
               DhMR_prefix = "UDH_ADH",
               comparison = "ADHvsUDH",
               region = "Intergenic")
make_pointplot(RNAseq_prefix = "DCIS_EN",
               DhMR_prefix = "ADH_DCIS",
               comparison = "DCISvsADH",
               region = "Intergenic")
make_pointplot(RNAseq_prefix = "IDC_DCIS",
               DhMR_prefix = "DCIS_PT",
               comparison = "IDCvsDCIS",
               region = "Intergenic")

make_pointplot(RNAseq_prefix = "EN_normal",
               DhMR_prefix = "UDH_ADH",
               comparison = "ADHvsUDH",
               region = c("Exon","Intron"))
make_pointplot(RNAseq_prefix = "DCIS_EN",
               DhMR_prefix = "ADH_DCIS",
               comparison = "DCISvsADH",
               region = c("Exon","Intron"))
make_pointplot(RNAseq_prefix = "IDC_DCIS",
               DhMR_prefix = "DCIS_PT",
               comparison = "IDCvsDCIS",
               region = c("Exon","Intron"))


## make venn plots for DEG and DMGs
library(VennDiagram)
make_pointplot <- 
  function(RNAseq_prefix,DhMR_prefix,title){
    RNAseq<- 
      read.table(paste0("E:\\20220921-WSL-5hmc\\analysis\\RNAseq_2014paper_DEGs\\",RNAseq_prefix,"_RNAseq_diff_gene.txt"),
                 header = T) %>%
      filter(abs(logFC)>1 & P.Value<0.05) %>%
      filter(!is.na(Row.names)) %>%
      pull(Row.names)
    
    DhMRs<-
      read.table(paste0("E:\\20220921-WSL-5hmc\\analysis\\ChIPseeker_diffPeak_anno\\",DhMR_prefix,"_diff_peak_deseq2.txt_anno"),
                 stringsAsFactors = F, header = T,
                 sep = "\t",comment.char = "",quote = "\"") %>%
      filter(abs(Fold)>1 & p.value<0.05) %>%
      filter(!is.na(SYMBOL)) %>%
      pull(SYMBOL)
    
    venn.diagram(
      x=list(DEG=RNAseq,
             DhMRs=DhMRs),
      disable.logging = T,
      filename = paste0("E:\\20220921-WSL-5hmc\\analysis\\corr_RNA_5hmC_5mC\\",title,"_venn.svg"),
      output=TRUE,
      cex = 2,
      main.cex = 3,
      fontfamily = "sans",
      cat.cex = 2,
      imagetype = "svg",
      main=DhMR_prefix,
      height = 10,
      width = 10)
  }

make_pointplot(RNAseq_prefix = "EN_normal",
               DhMR_prefix = "UDH_ADH_H",
               title = "ADH_UDH_DEG_DhMR")
make_pointplot(RNAseq_prefix = "DCIS_EN",
               DhMR_prefix = "ADH_DCIS_H",
               title = "DCIS_ADH_DEG_DhMR")
make_pointplot(RNAseq_prefix = "IDC_DCIS",
               DhMR_prefix = "DCIS_PT_H",
               title = "IDC_DCIS_DEG_DhMR")

make_pointplot(RNAseq_prefix = "EN_normal",
               DhMR_prefix = "UDH_ADH_M",
               title = "ADH_UDH_DEG_DMR")
make_pointplot(RNAseq_prefix = "DCIS_EN",
               DhMR_prefix = "ADH_DCIS_M",
               title = "DCIS_ADH_DEG_DMR")

## count the up- and down-regulated gene nubmer in each stage

do.call(rbind,lapply(c("EN_normal","DCIS_EN","IDC_DCIS"),
        function(x){
          RNAseq<- 
            read.table(paste0("E:\\20220921-WSL-5hmc\\analysis\\RNAseq_2014paper_DEGs\\",
                              x,
                              "_RNAseq_diff_gene.txt"),
                       header = T) %>%
            mutate(change=case_when(.default="no",logFC>1 & P.Value<0.05 ~ "up",
                                    logFC<(-1) & P.Value<0.05 ~ "down")) %>%
            group_by(change) %>%
            summarise(n()) %>%
            mutate(group=!!x)
        })) %>%
  as.data.frame() %>%
  write.table(file = "E:\\20220921-WSL-5hmc\\analysis\\RNAseq_2014paper_DEGs\\GEO_DEG_number_3phases.txt",
              quote = F,sep = "\t",
              row.names = F)
