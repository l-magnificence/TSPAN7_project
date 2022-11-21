setwd("~/bioinfo_mill/tspan7_proj")
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GSVA)
library(cowplot)
library(aplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(maftools)
library(corrplot)
library(survival)
library(survminer)
gmtPathways <- function(gmt.file) {
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  names(pathways) <- sapply(pathwayLines, head, 1)
  pathways
}
min.max.norm <- function(x){
  ((x-min(x))/(max(x)-min(x)))
}
gene<-"TSPAN7"
expr <- readRDS("~/bioinfo_mill/dataset/tcga_glioma/processed/tcga_exprSet_log2tpm_glioma_UCSCmatch2016cell.rds.rds")
meta <- readRDS("~/bioinfo_mill/dataset/tcga_glioma/processed/tcga_meta_glioma_match2016cell.rds")

gene_select<-as.data.frame(t(expr[gene,]))
meta$gene<-gene_select[,gene]
meta$group<-ifelse(meta$gene>median(meta$gene),"High","Low")
identical(rownames(meta),colnames(expr))

###-----------------DEG-----------------
### limma
group <- meta$group
group <- factor(group,levels = c("Low","High"))

design <- model.matrix(~group)
fit <- lmFit(expr,design)
fit2 <- eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 

allDiff$significant<-ifelse(allDiff$logFC>2&allDiff$adj.P.Val<0.05,"up",
                            ifelse(allDiff$logFC< -2&allDiff$adj.P.Val<0.05,"down","unsignificant"))
table(allDiff$significant)
saveRDS(allDiff,file = "./processed/tcga_glioma_allDiff.rds")

pdf(file="./result/volcano.pdf",bg="white",width=7.3,height=5.2,pointsize=8)
ggplot(data =allDiff,aes(x=logFC,y= -1*log10(adj.P.Val)))+geom_point(aes(color=significant),alpha=0.4, size=3.5)+
  scale_color_manual(values = c("blue","grey","red")) + xlim(c(-8.5, 8.5))+
  labs(title="Differentially expressed genes",x=expression((log[2](FC)),y=expression(-log[10](padj)) ))+
  #geom_hline(yintercept=1.5,linetype=4)+
  geom_vline(xintercept=c(-2,2),linetype=4)+
  theme_classic()+
  theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        legend.key.size = unit(0.8, "cm"),legend.text = element_text(size=16),legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5), legend.position="right")
dev.off()


DEG<-allDiff[allDiff$significant %in%c("up","down"),]
ids <- bitr(row.names(DEG), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db",drop = TRUE)
DEG$gene<-rownames(DEG)
DEG_list=merge(DEG,ids,by.x='gene',by.y='SYMBOL')
DEG_list<-DEG_list[order(DEG_list$logFC,decreasing = T),]
table(DEG_list$significant)

##gsea
geneList<-DEG_list$logFC
names(geneList)<-DEG_list$ENTREZID

gseGO <- gseGO(geneList = geneList, 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENTREZID", 
                ont="BP")
gseGO <-setReadable(gseGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
saveRDS(gseGO,file = "./processed/tcga_glioma_allDiff_gseGO.rds")

tmp<-as.data.frame(gseGO)
enrichment<-tmp[c(grep("MAPK cascade",tmp$Description),grep("regulation of cell migration",tmp$Description),grep("cell population proliferation",tmp$Description),
                  grep("stress",tmp$Description),grep("negative regulation of programmed cell death",tmp$Description),
                  grep("tumor necrosis factor",tmp$Description),grep("angiogenesis",tmp$Description),grep("transcription, DNA-templated",tmp$Description)),]
enrichment<-enrichment[-c(14,17),]
tmp<-rbind(tmp[order(tmp$NES,decreasing = T),][1:15,],enrichment[order(enrichment$NES,decreasing = F),][1:15,])
tmp$Description<-factor(tmp$Description,levels = c(tmp$Description))
tmp$group<-ifelse(tmp$NES>0,"UP","DOWN")

pdf(file="./result/gsea_gobp.pdf",bg="white",width=6.5,height=5.4,pointsize=12)
ggplot(data = tmp, aes(x =NES , y = Description, fill = group)) +
  scale_fill_manual(values=c("#32769B", "#BB173A"))+
  geom_bar(stat = "identity",position = "identity",color="NA",size=0.25) +
  theme_classic()+
  labs(title = "GSEA GO BP")+
  ylab("") +
  xlab("NES") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")
dev.off()

gseKEGG <- gseKEGG(geneList = geneList, 
                   organism = "hsa",
                   keyType = "kegg", 
                   pvalueCutoff =0.05)
gseKEGG <-setReadable(gseKEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
saveRDS(gseKEGG,file = "./processed/tcga_glioma_allDiff_gseKEGG.rds")

tmp2<-as.data.frame(gseKEGG)
pdf(file="./result/gseKEGG_heatplot.pdf",bg="white",width=13,height=2.3,pointsize=12)
heatplot(gseKEGG,showCategory =10,foldChange=geneList)+
  scale_fill_gradient2(low = "#336392",mid = "#BEC8CC",high = "#A94322",midpoint = 0)+
  labs(title = "GSEA KEGG")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

###--------------hallmark gsea--------------
hallmark <- gmtPathways("~/software/msigdb_v7.2/h.all.v7.2.symbols.gmt")
all_pathways2 <- gmtPathways('~/software/msigdb_v7.2/c2.cp.biocarta.v7.2.symbols.gmt')
all_pathways3 <- gmtPathways('~/software/msigdb_v7.2/c2.cp.kegg.v7.2.symbols.gmt')

for (i in 1:length(names(hallmark))) { 
  names(hallmark)[i]<-substring(names(hallmark)[i],10,nchar(names(hallmark)[i]))
}
Marker<-list("G2M CHECKPOINT"=hallmark$G2M_CHECKPOINT,"EMT"=hallmark$EPITHELIAL_MESENCHYMAL_TRANSITION,
             "ANGIOGENESIS"=hallmark$ANGIOGENESIS,"HYPOXIA"=hallmark$HYPOXIA,"DNA REPAIR"=hallmark$DNA_REPAIR,
             "MAPK"=all_pathways2[["BIOCARTA_MAPK_PATHWAY"]])
saveRDS(Marker,file = "./processed/gene_list.rds")

hallmark_gsea<-gsva(as.matrix(expr), Marker ,method="ssgsea")
hallmark_gsea<-as.data.frame(t(hallmark_gsea))
tmp<-hallmark_gsea
identical(rownames(meta),rownames(tmp))
hallmark_gsea<-cbind(tmp,meta)
saveRDS(hallmark_gsea,file = "./processed/tcga_hallmark_gsea.rds")

tmp<-hallmark_gsea
tmp[,1:length(Marker)]<-apply(tmp[,1:length(Marker)], 2, min.max.norm)
tmp$group<-factor(tmp$group,levels = c("High","Low"))

pdf(file="./result/tcga_pathway_compare.pdf",bg="white",width=5.4,height=4.5,pointsize=12)
ggviolin(tmp, x = "group",
         y = c(names(Marker)),
         combine = TRUE,ncol=3, 
         color = "group", palette = "nejm",
         ylab = "Enrichment score", xlab =F,
         add = "boxplot",    
         add.params = list(binwidth = 0.1, dotsize = 0.3))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 0, hjust = 0.5))+
  stat_compare_means(method="wilcox.test",label.x.npc= 'center',label = "p.signif")
dev.off()

## corelation
cor_data<-data.frame("hallmark"=1,"cor"=1,"p"=1)
for (i in 1:length(Marker)) {
  m<-cor.test(hallmark_gsea[,i],hallmark_gsea$gene)
  cor_data[i,1]<-colnames(hallmark_gsea)[i]
  cor_data[i,2]<-m$estimate
  cor_data[i,3]<-m$p.value
}
cor_data$logp<-(-log10(cor_data$p))

##heatmap
meta<-meta[order(meta$gene,decreasing = T),]

is_sig = cor_data$p < 0.01
pch = rep("*", 10)
pch[!is_sig] = NA
row_ha = rowAnnotation( Cor = anno_barplot(cor_data$cor,gp = gpar(fill =paletteer::paletteer_d("ggsci::default_nejm"))),
                        #`-log10(Pvalue)` = cor_data$logp,
                        Sig=anno_simple(cor_data$logp, pch=pch,col = circlize::colorRamp2(c(seq(15, 69, 1.8)[-30]), paletteer::paletteer_c("grDevices::Red-Purple", 30,direction =-1))),
                        col = list(`-log10(Pvalue)`=circlize::colorRamp2(c(seq(15, 69, 1.8)[-30]), paletteer::paletteer_c("grDevices::Red-Purple", 30,direction =-1))),
                        na_col = "white",border = F,
                        gap = unit(3, "points"),width= unit(2, "cm"),
                        annotation_name_side = "top",annotation_name_rot=0)
ha1 = HeatmapAnnotation(TSPAN7=meta$gene,
                        Age=ifelse(meta$`Age (years at diagnosis)`>60,">60","<60"),
                        Gender=meta$Gender,
                        Grade=meta$Grade,
                        Histology=meta$Histology,
                        Transcriptome_subtype =meta$`Transcriptome Subtype`,
                        IDH=meta$`IDH status`,
                        Chr1p19q=meta$`1p/19q codeletion`,
                        ATRX=meta$`TERT promoter status`,
                        MGMT=meta$`MGMT promoter status`,
                        #purity=anno_points(as.numeric(meta$`ABSOLUTE purity`)),
                        col = list(TSPAN7=circlize::colorRamp2(c(seq(2.5, 10.6, 0.27)[-30]), paletteer::paletteer_c("grDevices::Inferno", 30,direction = -1)),
                                   Age =structure(names = c(">60", "<60"), c("#FC8D62","#E5C494")),
                                   Grade=structure(names = c("G2", "G3","G4"), c("#FEE0D2","#FB6A4A","#A50F15")),
                                   Histology=structure(names = c("astrocytoma", "oligodendroglioma","oligoastrocytoma","glioblastoma"), brewer.pal(8, "Set1")[1:4]),
                                   Gender =structure(names = c("female", "male"),c("#6A3D9A","#FB9A99")),
                                   Transcriptome_subtype =structure(names = c("ME", "NE","PN","CL","NA"), c(brewer.pal(4, "Set2"),"white")),
                                   IDH =structure(names = c("Mutant", "WT","NA"), c("#C63F41", "#2574AF","white")),
                                   Chr1p19q=structure(names = c("codel", "non-codel","NA"), c("#C63F41", "#2574AF","white")),
                                   ATRX =structure(names = c("Mutant", "WT","NA"), c("#C63F41", "#2574AF","white")),
                                   MGMT =structure(names = c("Methylated", "Unmethylated","NA"), c("#C63F41", "#2574AF","white"))),
                        na_col = "white",border = F,
                        gap = unit(3, "points"),
                        height = unit(8, "cm")
)

dat<-hallmark_gsea[,1:length(Marker)]
dat<-dat[rownames(meta),]
dat<-scale(dat)

pdf(file="./result/tcga_pathway_heatmap.pdf",bg="white",width=12,height=6,pointsize=12)
Heatmap(t(dat),name="Enrichment",cluster_columns = F,cluster_rows =F,
        col=colorRampPalette(c("royalblue","white","firebrick3"))(30),
        show_column_names =F,row_names_gp = gpar(fontsize = 11),width  = unit(10, "cm"),
        height  = unit(5, "cm"),
        show_row_names =T,top_annotation =ha1,left_annotation = row_ha)
dev.off()

p<-Heatmap(t(dat),name="Enrichment",cluster_columns = F,cluster_rows =F,
           col=colorRampPalette(c("royalblue","white","firebrick3"))(30),
           show_column_names =F,row_names_gp = gpar(fontsize = 11),width  = unit(10, "cm"),
           height  = unit(5, "cm"),
           show_row_names =T,top_annotation =ha1,left_annotation = row_ha)

lgd_pvalue = Legend(title = "-log10(p-value)", col_fun = circlize::colorRamp2(c(seq(15, 69, 1.8)[-30]), paletteer::paletteer_c("grDevices::Red-Purple", 30,direction =-1)) , at = c(0, 20, 40, 60,80), 
                    labels = c("0", "20", "40", "60","80"))
# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.01")
# these two self-defined legends are added to the plot by `annotation_legend_list`
draw(p, annotation_legend_list = list(lgd_pvalue, lgd_sig))


pdf(file="./result/tcga_pathway_heatmap2.pdf",bg="white",width=12,height=6,pointsize=12)
draw(p, annotation_legend_list = list(lgd_pvalue, lgd_sig))
dev.off()

###--------------Immune infiltration landscape--------------
cibersort_result<-readRDS("~/bioinfo_mill/lxj_glioma_proj/processed_data/tcga_cibersort_result.rds")
cibersort_result<-as.data.frame(cibersort_result)
cibersort_result<-cibersort_result[rownames(meta),]
cibersort_result<-cibersort_result[,1:22]
cibersort_result$group<-meta$group

tmp<-reshape2::melt(cibersort_result)

pdf(file="./result/tspan7_cibersort_dist.pdf",bg="white",width=12,height=4.5,pointsize=12)
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

pdf(file="./result/tspan7_cibersort_dist2.pdf",bg="white",width=16,height=5,pointsize=12)
ggplot(tmp, aes(x=value,..scaled.., fill=group)) +
  geom_density(alpha=0.3)+
  facet_wrap(~variable,nrow =  3)+
  scale_fill_manual(values=c( "#BB173A","#32769B"))+
  theme_classic()+
  labs(title = "TCGA glioma datasets")+
  ylab("Density") +
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

###--------------m1 m2 --------------
source("~/scRNA_GBM_integration/tacan_project/code/function.R")
tpm_expr<-2^expr-0.001
tmp<-t(log2(tpm_expr+1))
tmp<-as.data.frame(tmp)

pdf(file="./result/macrophage/M2_CD163.pdf",bg="white",width=3.5,height=3.5,pointsize=12)
cor_plot(tmp,"TSPAN7","CD163","spearman")
dev.off()

pdf(file="./result/macrophage/M2_TGFBI.pdf",bg="white",width=3.5,height=3.5,pointsize=12)
cor_plot(tmp,"TSPAN7","TGFBI","spearman")
dev.off()

pdf(file="./result/macrophage/M2_IL10.pdf",bg="white",width=3.5,height=3.5,pointsize=12)
cor_plot(tmp,"TSPAN7","IL10","spearman")
dev.off()

pdf(file="./result/macrophage/M1_TNF.pdf",bg="white",width=3.5,height=3.5,pointsize=12)
cor_plot(tmp,"TSPAN7","TNF","spearman")
dev.off()

pdf(file="./result/macrophage/M1_NOS2.pdf",bg="white",width=3.5,height=3.5,pointsize=12)
cor_plot(tmp,"TSPAN7","NOS2","spearman")
dev.off()

pdf(file="./result/macrophage/M1_IL1B.pdf",bg="white",width=3.5,height=3.5,pointsize=12)
cor_plot(tmp,"TSPAN7","IL1B","spearman")
dev.off()

###-----------------correlation between tspan7 group and inhibitory immune checkpoints-----------------
expr <- readRDS("~/bioinfo_mill/dataset/tcga_glioma/processed/tcga_exprSet_log2tpm1_glioma_UCSCmatch2016cell.rds")
meta <- readRDS("~/bioinfo_mill/dataset/tcga_glioma/processed/tcga_meta_glioma_match2016cell.rds")

gene<-"TSPAN7"
gene_select<-as.data.frame(t(expr[gene,]))
meta$gene<-gene_select[,gene]
meta$group<-ifelse(meta$gene>median(meta$gene),"High","Low")
identical(rownames(meta),colnames(expr))

immune_gene<-c("HAVCR2", #TIM-3
               "CD274", #PD-L1
               "CD86","LAG3","LAIR1","PVR","IDO1","CD80","CTLA4",
               "PDCD1", #PD-1
               "TIGIT","CD200R1","CEACAM1","CD276","CD200","KIR3DL1","BTLA","ADORA2A","LGALS3","VTCN1"
) # got from paper PMID: 33537076

immune_expr<-as.data.frame(t(expr[immune_gene,]))
identical(rownames(immune_expr),rownames(meta))
immune_expr$group<-meta$group
immune_expr<-reshape2::melt(immune_expr)

pdf(file="./result/immuotherapy/tspan7_ICB_expr.pdf",bg="white",width=12,height=4.5,pointsize=12)
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
identical(rownames(immune_expr),rownames(meta))
immune_expr$gene<-as.data.frame(t(expr[gene,]))[,1]
immune_expr2<-reshape2::melt(immune_expr,"gene")

pdf(file="./result/immuotherapy/tspan7_ICB_cor.pdf",bg="white",width=14,height=7,pointsize=12)
ggplot(immune_expr2,aes(gene,value))+
  geom_point(size=2,alpha =0.3)+
  geom_smooth(method="lm",color="#fdae61")+
  stat_cor(method = "spearman",color = "#585658")+
  facet_wrap(~variable,nrow =  3)+
  theme_classic()+
  labs(title = "TCGA glioma datasets")+
  ylab("Relative expression of inhibitory immune checkpoints") +
  xlab("TSPAN7 expression")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



