
## plot the expression of tet2 in each 
RNAseq_df <- 
read.table("E:\\20220921-WSL-5hmc\\pre-info\\13059_2014_3265_MOESM4_ESM.txt",
           skip = 2,sep = "\t",comment.char = "",
           row.names = 1)
meta_infor <- read.table("E:\\20220921-WSL-5hmc\\pre-info\\13059_2014_3265_MOESM4_ESM.txt",
                         nrows = 2,
                         sep = "\t",comment.char = "")

tet2_expr <- 
  data.frame(value=unlist(RNAseq_df["TET2",]),
             group=unlist(meta_infor[2,-1])) %>%
  filter(!group=="nl")
tet2_expr$group <- factor(tet2_expr$group,
                          levels = c("normal","en","dcis","idc"))
ggplot(tet2_expr,aes(group,value,fill=group))+
  geom_violin()+
  geom_point()+
  guides()