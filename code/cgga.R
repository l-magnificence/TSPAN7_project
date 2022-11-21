setwd("~/bioinfo_mill/tspan7_proj")
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)

gene<-"TSPAN7"
batch<-"mRNA_325" 
type<-"Primary"
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
cgga_expr <- readRDS("~/bioinfo_mill/dataset/CGGA_data/processed_data/cgga_expr.rds")
expr <-2^cgga_expr-1
expr <- apply(expr,2,fpkmToTpm)
expr <- log2(expr+1)

meta <- readRDS("~/bioinfo_mill/dataset/CGGA_data/processed_data/cgga_clinic.rds")
meta<-meta[-which(is.na(meta$Grade)),]
meta[meta$Grade=="WHO II","Grade"]<-"G2"
meta[meta$Grade=="WHO III","Grade"]<-"G3"
meta[meta$Grade=="WHO IV","Grade"]<-"G4"

meta<-meta[meta$batch==batch & meta$PRS_type==type,]
expr<-expr[,meta$CGGA_ID]
expr<-as.data.frame(expr)
gene_select<-as.data.frame(t(expr[gene,]))
meta$gene<-gene_select[,gene]
meta$group<-ifelse(meta$gene>median(meta$gene),"High","Low")
identical(meta$CGGA_ID,colnames(expr))

###--------------Immune infiltration landscape--------------
cibersort_result<-readRDS(file="~/bioinfo_mill/lxj_glioma_proj/processed_data/cgga/cgga_cibersort_result.rds")
identical(meta$CGGA_ID,rownames(cibersort_result))

cibersort_result<-cibersort_result[,1:22]
cibersort_result<-as.data.frame(cibersort_result)
cibersort_result$group<-meta$group

tmp<-reshape2::melt(cibersort_result)

pdf(file="./result/cgga_325/tspan7_cibersort_dist.pdf",bg="white",width=12,height=4.5,pointsize=12)
ggplot(tmp, aes(variable, value,fill =group))+
  scale_fill_manual(values=c( "#BB173A","#32769B"))+
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  # geom_point(shape = 21, size=1, # 点的性状和大小
  #            position = position_jitterdodge(), # 让点散开
  #            color="black", alpha = 0.2) +
  theme_classic()+
  stat_compare_means(aes(group = group),method = "wilcox.test",hide.ns = F,label = "p.signif") +
  labs(fill = "")+
  ylab("Relative immune cell infiltration") +
  xlab("") +
  #ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))
dev.off()

pdf(file="./result/cgga_325/tspan7_cibersort_dist2.pdf",bg="white",width=16,height=5,pointsize=12)
ggplot(tmp, aes(x=value,..scaled.., fill=group)) +
  geom_density(alpha=0.3)+
  facet_wrap(~variable,nrow =  3)+
  scale_fill_manual(values=c( "#BB173A","#32769B"))+
  theme_classic()+
  labs(title = "CGGA glioma datasets")+
  ylab("Density") +
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

###-----------------correlation between tspan7 group and inhibitory immune checkpoints-----------------
immune_gene<-c("HAVCR2", #TIM-3
               "CD274", #PD-L1
               "CD86","LAG3","LAIR1","PVR","IDO1","CD80","CTLA4",
               "PDCD1", #PD-1
               "TIGIT","CD200R1","CEACAM1","CD276","CD200",#"KIR3DL1",
               "BTLA","ADORA2A","LGALS3","VTCN1"
) # got from paper PMID: 33537076

immune_expr<-as.data.frame(t(expr[immune_gene,]))
identical(rownames(immune_expr),meta$CGGA_ID)
immune_expr$group<-meta$group
immune_expr<-reshape2::melt(immune_expr)

pdf(file="./result/cgga_325/cgga_325_tspan7_ICB_expr.pdf",bg="white",width=12,height=4.5,pointsize=12)
ggplot(immune_expr, aes(variable, value,fill =group))+
  scale_fill_manual(values=c( "#BB173A","#32769B"))+
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  # geom_point(shape = 21, size=1, # 点的性状和大小
  #            position = position_jitterdodge(), # 让点散开
  #            color="black", alpha = 0.2) +
  theme_classic()+
  stat_compare_means(aes(group = group),method = "wilcox.test",hide.ns = F,label = "p.signif") +
  labs(fill = "")+
  ylab("Relative gene expression level") +
  xlab("") +
  #ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))
dev.off()

immune_expr<-as.data.frame(t(expr[immune_gene,]))
identical(rownames(immune_expr),meta$CGGA_ID)
immune_expr$gene<-as.data.frame(t(expr[gene,]))[,1]
immune_expr2<-reshape2::melt(immune_expr,"gene")

pdf(file="./result/cgga_325/cgga_325_tspan7_ICB_cor.pdf",bg="white",width=14,height=7,pointsize=12)
ggplot(immune_expr2,aes(gene,value))+
  geom_point(size=2,alpha =0.3)+
  geom_smooth(method="lm",color="#fdae61")+
  stat_cor(method = "spearman",color = "#585658")+
  facet_wrap(~variable,nrow =  3)+
  theme_classic()+
  labs(title = "CGGA glioma datasets")+
  ylab("Relative expression of inhibitory immune checkpoints") +
  xlab("TSPAN7 expression")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
