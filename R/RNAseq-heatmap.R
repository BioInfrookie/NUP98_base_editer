library(ggplot2)
library(edgeR)
library(dplyr)

rm(list = ls())
options(stringsAsFactors = F)

ID_trans <- read.table("/home/qzh/ref/Saccharomyces/Ensemble_id_trans.txt",
                       header = T,sep = "\\",quote = "",comment.char = "")
ID_trans <- ID_trans[,c(1:3)]
names(ID_trans) <- c("genes","GENENAME","SYMBOL")
setwd("/home/qzh/Analysis/Shengwang/YY539_RNAseq_Nup/02-results/")
Gene_count = read.table("./YY539_merge_count.txt",
                        header = T,sep = "\t")
Gene_count <- Gene_count[-c(1,2),]

Gene_merge_count <- Gene_count[,c(1,2,4,3)]
rownames(Gene_merge_count) <- Gene_merge_count[,1]
Gene_merge_count <- Gene_merge_count[,-1]
Gene_merge_count <- Gene_merge_count[which(rowSums(Gene_merge_count) > 0),]
group <- c("WT","Mut","Mut")
DEG_list <- DGEList(counts = Gene_merge_count[,c(1:3)],
                    genes = rownames(Gene_merge_count))
Gene_merge_cpm <- data.frame(log10(cpm(DEG_list)))

Gene_merge_cpm <- Gene_merge_cpm %>%
  mutate_all(~ replace(., is.infinite(.), -3))
Gene_merge_cpm[Gene_merge_cpm>=3] <- 3
Gene_merge_cpm[Gene_merge_cpm<=-3] <- -3

names(Gene_merge_cpm) <- c("A105_background","A107_DEN","A107_DEN_EST")
heat_plot <- pheatmap::pheatmap(Gene_merge_cpm,
                                cluster_rows = T,
                                cluster_cols = F,
                                fontsize_col=10,
                                fontsize_row = 10,
                                show_rownames = F,
                                show_colnames = T,
                                border_color="white",
                                legend_breaks=seq(-3,3,1.5),
                                # legend_labels = "log10(CPM)",
                                color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                angle_col = "90")

ggsave(heat_plot,filename="./Nup98_gene_cpm_heatmap.pdf",width = 3,height = 7)
