library(dplyr)

RNAseq_dir <- "E:\\20220921-WSL-5hmc\\analysis\\RNAseq_2014paper_DEGs\\"

continuous_high <- 
  read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuously_high_peak.txt",
           stringsAsFactors = F, header = T,
           sep = "\t",comment.char = "",quote = "\"")
continuous_low <- 
  read.table("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\continuously_low_peak.txt",
             stringsAsFactors = F, header = T,
             sep = "\t",comment.char = "",quote = "\"")

make_sankey_plot <-
function(continuous_high, file_name){
    
  Stage_I <- 
    read.table(paste0(RNAseq_dir,"EN_normal_RNAseq_diff_gene.txt"),
               header = T)
  Stage_II <- 
    read.table(paste0(RNAseq_dir,"DCIS_EN_RNAseq_diff_gene.txt"),
               header = T)
  Stage_III <- 
    read.table(paste0(RNAseq_dir,"IDC_DCIS_RNAseq_diff_gene.txt"),
               header = T)
  
  data_for_plot <- 
  data.frame(SYMBOL = Stage_I$Row.names[match(continuous_high$SYMBOL,Stage_I$Row.names,nomatch = 0)],
             Stage_I[match(continuous_high$SYMBOL,Stage_I$Row.names,nomatch = 0),c("logFC","P.Value")],
             Stage_II[match(continuous_high$SYMBOL,Stage_II$Row.names,nomatch = 0),c("logFC","P.Value")],
             Stage_III[match(continuous_high$SYMBOL,Stage_III$Row.names,nomatch = 0),c("logFC","P.Value")])
  head(data_for_plot)
  colnames(data_for_plot) <-
    c("SYMBOL",paste0(rep(c("StageI_","StageII_","StageIII_"),each=2),c("logFC","P.value")))
  dim(data_for_plot)
  # 1887  7
  data_for_plot <- 
    data_for_plot %>%
      mutate(StageI_sign = case_when(
        StageI_logFC<(-1) & StageI_P.value<0.05 ~ "down",
        StageI_logFC>1 & StageI_P.value<0.05 ~ "up",
        .default = "no-sig"
      )) %>%
      mutate(StageII_sign = case_when(
        StageII_logFC<(-1) & StageII_P.value<0.05 ~ "down",
        StageII_logFC>1 & StageII_P.value<0.05 ~ "up",
        .default = "no-sig"
      )) %>%
      mutate(StageIII_sign = case_when(
        StageIII_logFC<(-1) & StageIII_P.value<0.05 ~ "down",
        StageIII_logFC>1 & StageIII_P.value<0.05 ~ "up",
        .default = "no-sig"
      ))
  head(data_for_plot)  
  
  data_for_plot %>%
    select(SYMBOL,StageI_sign,StageII_sign,StageIII_sign) %>%
    filter((StageI_sign=="no-sig")+(StageII_sign=="no-sig")+(StageIII_sign=="no-sig") <3) %>%
    reshape2::melt(id.vars="SYMBOL") %>%
    as.data.frame() %>%
  ggplot(aes(x=variable,y=SYMBOL,fill=value))+
    geom_tile()+
    scale_fill_manual(values=c("blue","grey","red"))
  
  ## Sankey plot
  library(ggsankey)
  library(ggplot2)
  
  df_skey <- data_for_plot %>%
    select(StageI_sign,StageII_sign,StageIII_sign) %>%
    filter(rowSums(.=="no-sig")<3) %>%
    make_long(StageI_sign,StageII_sign,StageIII_sign)
  
  head(df_skey)
  (ggplot(df_skey, aes(x = x
                      , next_x = next_x
                      , node = node
                      , next_node = next_node
                      , fill = factor(node)
                      , label = node)) +
    geom_sankey(flow.alpha = 0.5
                ,node.color = "grey"
                ,show.legend = TRUE)+
    scale_fill_manual(values = c("#2a9d8f","grey","#e76f51"),name="RNA level\nchange")+
    scale_x_discrete(name="",labels=c("Stage I","Stage II","Stage III"))+
    theme(panel.background = element_rect(fill = "white"),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 11,colour = "black"))) |>
    ggsave(filename = paste0("E:\\20220921-WSL-5hmc\\analysis\\make_diff_heatmap\\",file_name),
           width = 6,height = 4)
}


make_sankey_plot(continuous_high = continuous_high,
                 file_name = "continuous_high_peak_RNA_change.pdf")
make_sankey_plot(continuous_high = continuous_low,
                 file_name = "continuous_low_peak_RNA_change.pdf")
