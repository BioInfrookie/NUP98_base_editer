library(ggplot2)
library(edgeR)
library(org.Sc.sgd.db)
library(clusterProfiler)
library(AnnotationHub)

rm(list = ls())
options(stringsAsFactors = F)

ID_trans <- read.table("/home/qzh/ref/Saccharomyces/Ensemble_id_trans.txt",
                       header = T,sep = "\\",quote = "",comment.char = "")
ID_trans <- ID_trans[,c(1:3)]
names(ID_trans) <- c("genes","GENENAME","SYMBOL")
Gene_count = read.table("./A04_count.txt",
                        header = T,sep = "\t")

# Dox 100
A04_100_table <- Gene_count[,c(1,5,2)]
rownames(A04_100_table) <- A04_100_table[,1]
A04_100_table <- A04_100_table[,-1]
A04_100_table <- A04_100_table[which(rowSums(A04_100_table) > 0),]
names(A04_100_table) <- c("CK","Dox_100")
group <- c("WT","Mut")
DEG_list <- DGEList(counts = A04_100_table[,c(1,2)],
                    genes = rownames(A04_100_table),group = group)
keep <- rowSums(cpm(DEG_list)>1) >= 1
DEG_list <- DEG_list[keep, , keep.lib.sizes=FALSE]

DEG_list <- calcNormFactors(DEG_list)
DEG_list_bcv <- DEG_list
DEG_cpm <- data.frame(cpm(DEG_list))
DEG_cpm$genes <- rownames(DEG_cpm)
bcv <- 0.4
et <- exactTest(DEG_list_bcv, dispersion = bcv ^ 2)
DEG <- cbind(et$genes,et$table)
DEG$logFC <- DEG$logFC*-1
DEG$change = as.factor(ifelse(abs(DEG$logFC) >=1,
                              ifelse(DEG$logFC >= 1,'UP','DOWN'),'NOT'))

DEG_order <- DEG[order(DEG$logFC,decreasing = T),]
DEG_order$Rank <- round(rank(-DEG_order$logFC,ties.method = "first"))
DEG_table <- left_join(DEG_order,DEG_cpm,by="genes")
DEG_table <- left_join(DEG_table,ID_trans,by="genes")
# write.table(DEG_table,file = "A104_vs_A04_DEGlist.txt",sep = "\t",
#             col.names = T,row.names = F,quote = F)
this_tile <- paste0('\nThe number of up gene is ',nrow(DEG_table[DEG_table$change == 'UP',]),
                    '\nThe number of down gene is ',nrow(DEG_table[DEG_table$change == 'DOWN',]))

Scatter2 <- ggplot(data=DEG_table, aes(x=log10(CK), y=log10(Dox_100),color=change)) +
  geom_point(data=DEG_table[which(DEG_table$change=="NOT"),],
             aes(x=log10(CK), y=log10(Dox_100)),size=1,fill="grey") +
  geom_point(data=DEG_table[which(DEG_table$change=="UP"),], 
             aes(x=log10(CK), y=log10(Dox_100)),size=1.5,fill="red") +
  geom_point(data=DEG_table[which(DEG_table$change=="DOWN"),],
             aes(x=log10(CK), y=log10(Dox_100)),size=1.5,fill="darkgreen") +
  ylab(paste0("log10(CPM) Dox_100")) +
  xlab(paste0("log10(CPM) CK")) +
  # ggtitle(this_tile) +
  theme_classic() +
  theme(legend.position = "none") + 
  scale_x_continuous(expand = c(0,0),limits = c(0,5)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,5)) +
  geom_abline(intercept = c(log10(2),log10(0.5)), slope = 1, linetype="dotted", colour="black") +
  geom_abline(intercept = 0, slope = 1 , colour="black") +
  scale_color_manual(values=c(NOT = "grey", UP = "red",
                              DOWN="darkgreen")) + 
  theme(
    axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10),
    axis.line = element_line(linewidth =0),
    panel.border = element_rect(color = "black", size = 0.8, fill = NA))
Scatter2
ggsave(Scatter2,filename="./A04_Dox_100_vs_CK.pdf",width = 3,height = 3)


# Dox_50
A04_50_table <- Gene_count[,c(1,5,4)]
rownames(A04_50_table) <- A04_50_table[,1]
A04_50_table <- A04_50_table[,-1]
A04_50_table <- A04_50_table[which(rowSums(A04_50_table) > 0),]
names(A04_50_table) <- c("CK","Dox_50")
group <- c("WT","Mut")
DEG_list <- DGEList(counts = A04_50_table[,c(1,2)],genes = rownames(A04_50_table),group = group)
keep <- rowSums(cpm(DEG_list)>1) >= 1
DEG_list <- DEG_list[keep, , keep.lib.sizes=FALSE]

DEG_list <- calcNormFactors(DEG_list)
DEG_list_bcv <- DEG_list
DEG_cpm <- data.frame(cpm(DEG_list))
DEG_cpm$genes <- rownames(DEG_cpm)
bcv <- 0.4
et <- exactTest(DEG_list_bcv, dispersion = bcv ^ 2)
DEG <- cbind(et$genes,et$table)
DEG$logFC <- DEG$logFC*-1
DEG$change = as.factor(ifelse(abs(DEG$logFC) >=1 ,
                              ifelse(DEG$logFC >= 1,'UP','DOWN'),'NOT'))

DEG_order <- DEG[order(DEG$logFC,decreasing = T),]
DEG_order$Rank <- round(rank(-DEG_order$logFC,ties.method = "first"))
DEG_table <- left_join(DEG_order,DEG_cpm,by="genes")
DEG_table <- left_join(DEG_table,ID_trans,by="genes")
# write.table(DEG_table,file = "A104_vs_A04_DEGlist.txt",sep = "\t",
#             col.names = T,row.names = F,quote = F)
this_tile <- paste0('\nThe number of up gene is ',nrow(DEG_table[DEG_table$change == 'UP',]),
                    '\nThe number of down gene is ',nrow(DEG_table[DEG_table$change == 'DOWN',]))

Scatter1 <- ggplot(data=DEG_table, aes(x=log10(CK), y=log10(Dox_50),color=change)) +
  geom_point(data=DEG_table[which(DEG_table$change=="NOT"),],
             aes(log10(CK), y=log10(Dox_50)),size=1,fill="grey") +
  geom_point(data=DEG_table[which(DEG_table$change=="UP"),], 
             aes(log10(CK), y=log10(Dox_50)),size=1.5,fill="red") +
  geom_point(data=DEG_table[which(DEG_table$change=="DOWN"),],
             aes(log10(CK), y=log10(Dox_50)),size=1.5,fill="darkgreen") +
  ylab(paste0("log10(CPM) Dox 50")) +
  xlab(paste0("log10(CPM) CK")) +
  # ggtitle(this_tile) +
  theme_classic() +
  theme(legend.position = "none") + 
  scale_x_continuous(expand = c(0,0),limits = c(0,5)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,5)) +
  geom_abline(intercept = c(log10(2),log10(0.5)), slope = 1, linetype="dotted", colour="black") +
  geom_abline(intercept = 0, slope = 1 , colour="black") +
  scale_color_manual(values=c(NOT = "grey", UP = "red",
                              DOWN="darkgreen")) + 
  theme(
    axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10),
    axis.line = element_line(linewidth =0),
    panel.border = element_rect(color = "black", size = 0.8, fill = NA))
Scatter1
ggsave(Scatter1,filename="./A04_Dox50_vs_CK_MA.pdf",width = 3,height = 3)

# Dox_200
A04_200_table <- Gene_count[,c(1,5,3)]
rownames(A04_200_table) <- A04_200_table[,1]
A04_200_table <- A04_200_table[,-1]
A04_200_table <- A04_200_table[which(rowSums(A04_200_table) > 0),]
names(A04_200_table) <- c("CK","Dox_200")
group <- c("WT","Mut")
DEG_list <- DGEList(counts = A04_200_table[,c(1,2)],genes = rownames(A04_200_table),group = group)
keep <- rowSums(cpm(DEG_list)>1) >= 1
DEG_list <- DEG_list[keep, , keep.lib.sizes=FALSE]

DEG_list <- calcNormFactors(DEG_list)
DEG_list_bcv <- DEG_list
DEG_cpm <- data.frame(cpm(DEG_list))
DEG_cpm$genes <- rownames(DEG_cpm)
bcv <- 0.4
et <- exactTest(DEG_list_bcv, dispersion = bcv ^ 2)
DEG <- cbind(et$genes,et$table)
DEG$logFC <- DEG$logFC*-1
DEG$change = as.factor(ifelse(abs(DEG$logFC) >=1 ,
                              ifelse(DEG$logFC >= 1,'UP','DOWN'),'NOT'))

DEG_order <- DEG[order(DEG$logFC,decreasing = T),]
DEG_order$Rank <- round(rank(-DEG_order$logFC,ties.method = "first"))
DEG_table <- left_join(DEG_order,DEG_cpm,by="genes")
DEG_table <- left_join(DEG_table,ID_trans,by="genes")
# write.table(DEG_table,file = "A104_vs_A04_DEGlist.txt",sep = "\t",
#             col.names = T,row.names = F,quote = F)
this_tile <- paste0('\nThe number of up gene is ',nrow(DEG_table[DEG_table$change == 'UP',]),
                    '\nThe number of down gene is ',nrow(DEG_table[DEG_table$change == 'DOWN',]))

Scatter3 <- ggplot(data=DEG_table, aes(x=log10(CK), y=log10(Dox_200),color=change)) +
  geom_point(data=DEG_table[which(DEG_table$change=="NOT"),],
             aes(log10(CK), y=log10(Dox_200)),size=1,fill="grey") +
  geom_point(data=DEG_table[which(DEG_table$change=="UP"),], 
             aes(log10(CK), y=log10(Dox_200)),size=1.5,fill="red") +
  geom_point(data=DEG_table[which(DEG_table$change=="DOWN"),],
             aes(log10(CK), y=log10(Dox_200)),size=1.5,fill="darkgreen") +
  ylab(paste0("log10(CPM) Dox 200")) +
  xlab(paste0("log10(CPM) CK")) +
  # ggtitle(this_tile) +
  theme_classic() +
  theme(legend.position = "none") + 
  scale_x_continuous(expand = c(0,0),limits = c(0,5)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,5)) +
  geom_abline(intercept = c(log10(2),log10(0.5)), slope = 1, linetype="dotted", colour="black") +
  geom_abline(intercept = 0, slope = 1 , colour="black") +
  scale_color_manual(values=c(NOT = "grey", UP = "red",
                              DOWN="darkgreen")) + 
  theme(
    axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10),
    axis.line = element_line(linewidth =0),
    panel.border = element_rect(color = "black", size = 0.8, fill = NA))
Scatter3
ggsave(Scatter3,filename="./A04_Dox200_vs_CK_MA.pdf",width = 3,height = 3)

up_gene <- DEG_table[which(DEG_table$change %in% "UP"),]
down_gene <- DEG_table[which(DEG_table$change %in% "DOWN"),]
diff_gene <- DEG_table[which(!(DEG_table$change %in% "NOT")),]
diff_gene <- diff_gene[order(abs(diff_gene$logFC),decreasing = T),]
diff_gene <- diff_gene[c(1:10),]

Scatter <- Scatter1 + Scatter2 + Scatter3
ggsave(Scatter,filename="./A04_Dox_vs_CK_MA.pdf",width = 9,height = 3)


