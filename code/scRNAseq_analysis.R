setwd("~/bioinfo_mill/tspan7_proj")

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(CellChat)
library(Scissor)
library(ggplot2)
library(ggpubr)
source("~/scRNA_GBM_integration/tacan_project/code/function.R")
min.max.norm <- function(x){
  ((x-min(x))/(max(x)-min(x)))
}
gene<-"TSPAN7"

###------------GSE84465------------
GSE84465_seurat <- readRDS("~/scRNA_GBM_integration/result/GSE84465_seurat_reassign.rds")
GSE84465_seurat <-subset(GSE84465_seurat,cell_type2 !="Unassigned immune cell")

GSE84465_seurat@meta.data$cell_type3 <- GSE84465_seurat@meta.data$cell_type2
GSE84465_seurat@meta.data$cell_type3 <- dplyr::recode(GSE84465_seurat@meta.data$cell_type3,
                                                       "Macrophage"="TAM",
                                                       "Microglia"="TAM")

pdf(file="./result/scRNA/GSE84465_celltype_tsne.pdf",bg="white",width=7,height=5.5,pointsize=8)
DimPlot(GSE84465_seurat, reduction = "tsne",group.by = "cell_type3",label = F)+labs(title="")+
  labs(title="GSE84465")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file="./result/scRNA/GSE84465_gene_expr.pdf",bg="white",width=6,height=3,pointsize=8)
StackedVlnPlot(GSE84465_seurat, gene, group.by="cell_type3",pt.size=0, cols=my36colors,ncol = 1)+coord_flip()
dev.off()

pdf(file="./result/scRNA/GSE84465_gene_expr_tsne.pdf",bg="white",width=4.7,height=4.1,pointsize=8)
FeaturePlot(GSE84465_seurat,features = gene, reduction = "tsne",cols =paletteer::paletteer_c("grDevices::OrRd", 5 ,direction = -1))
dev.off()

## subset tumor 
unique(GSE84465_seurat$cell_type3)
GSE84465_tumor<-subset(GSE84465_seurat,cell_type3=="Tumor cell")
DimPlot(GSE84465_tumor,reduction = "tsne",group.by = "cell_type3")

##define gene high express and low tumor
expr<-GSE84465_tumor@assays[["SCT"]]@data
expr<-as.data.frame(t(as.matrix(expr)))
expr<-dplyr::select(expr,"TSPAN7")
median(expr$TSPAN7)
mean(expr$TSPAN7)
TSPAN7_high_cell <- WhichCells(object = GSE84465_tumor, expression = TSPAN7 >=median(expr$TSPAN7), slot = "data")
TSPAN7_low_cell <- WhichCells(object = GSE84465_tumor, expression = TSPAN7 < median(expr$TSPAN7), slot = "data")

meta<-GSE84465_tumor@meta.data
meta$TSPAN7_group<-"a"
meta[TSPAN7_high_cell,"TSPAN7_group"]<-"TSPAN7_high"
meta[TSPAN7_low_cell,"TSPAN7_group"]<-"TSPAN7_low"
GSE84465_tumor$TSPAN7_group<-meta$TSPAN7_group
saveRDS(GSE84465_tumor,file = "~/bioinfo_mill/tspan7_proj/processed/GSE84465_tumor.rds")

meta<-GSE84465_seurat@meta.data
meta$cell_type4<-meta$cell_type3
meta$cell_type4<-as.character(meta$cell_type4)
meta[TSPAN7_high_cell,"cell_type4"]<-paste(gene,"_high",sep = "")
meta[TSPAN7_low_cell,"cell_type4"]<-paste(gene,"_low",sep = "")
table(meta$cell_type4)

GSE84465_seurat$cell_type4<-meta$cell_type4
saveRDS(GSE84465_seurat,file = "./processed/GSE84465_seurat.rds")

Idents(GSE84465_tumor)<-"TSPAN7_group"
GSE84465_tumor_DEG <- FindMarkers(GSE84465_tumor, ident.1="TSPAN7_high", ident.2 ="TSPAN7_low",
                                  only.pos = F, min.pct = 0.25, logfc.threshold = 0.5,test.use ="wilcox")
saveRDS(GSE84465_tumor_DEG,file = "~/bioinfo_mill/tspan7_proj/processed/GSE84465_tumor_DEG.rds")

## gsea
# DEG<-GSE84465_tumor_DEG[abs(GSE84465_tumor_DEG$avg_log2FC)>0 & GSE84465_tumor_DEG$p_val_adj<0.05 ,]
# ids <- bitr(row.names(DEG), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db",drop = TRUE)
# DEG$gene<-rownames(DEG)
# DEG_list=merge(DEG,ids,by.x='gene',by.y='SYMBOL')
# DEG_list<-DEG_list[order(DEG_list$avg_log2FC,decreasing = T),]
# DEG_list$significant<-ifelse(DEG_list$avg_log2FC>0.5 & DEG_list$p_val_adj<0.05,"up",
#                             ifelse(DEG_list$avg_log2FC< -0.5 & DEG_list$p_val_adj<0.05,"down","unsignificant"))
# table(DEG_list$significant)
# 
# enrichKEGG <- compareCluster(ENTREZID~significant, data=DEG_list, fun="enrichKEGG",
#                              organism="hsa", pvalueCutoff=0.05)
# enrichKEGG <- pairwise_termsim(enrichKEGG)   
# emapplot(enrichKEGG,showCategory =20,layout="kk", legend_n=3,repel=T,min_edge=0.2,cex_line=0.6)+
#   labs(title="KEGG enrichment", x="", y="")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# enrichGO <- compareCluster(ENTREZID~significant, data=DEG_list, fun ="enrichGO",ont="BP",OrgDb='org.Hs.eg.db') 
# enrichGO2 <- simplify(enrichGO, cutoff=0.7, by="p.adjust", select_fun=min)
# enrichGO2 <- pairwise_termsim(enrichGO)
# emapplot(enrichGO2,showCategory =20,layout="kk", legend_n=3,repel=T,min_edge=0.2,cex_line=0.6)+
#   labs(title="GO enrichment", x="", y="")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# dotplot(enrichGO, showCategory=15)

###------------GSE84465 Scissor------------
GSE84465_seurat<-readRDS("./processed/GSE84465_seurat.rds")
tumor<-subset(GSE84465_seurat,cell_type3=="Tumor cell")

expr <- readRDS("~/bioinfo_mill/dataset/tcga_glioma/processed/tcga_exprSet_log2tpm_glioma_UCSCmatch2016cell.rds.rds")
meta <- readRDS("~/bioinfo_mill/dataset/tcga_glioma/processed/tcga_meta_glioma_match2016cell.rds")

phenotype <- dplyr::select(meta,c("Survival (months)","Vital status (1=dead)"))
phenotype$`Survival (months)`<-as.numeric(phenotype$`Survival (months)`)
phenotype$`Vital status (1=dead)`<-as.double(phenotype$`Vital status (1=dead)`)
colnames(phenotype) <- c("time", "status")
head(phenotype)
identical(rownames(phenotype),colnames(expr))
phenotype<-as.matrix(phenotype)

GSE84465_tumor<-subset(GSE84465_seurat,cell_type3=="Tumor cell")
GSE84465_tumor<- Seurat_preprocessing(GSE84465_tumor@assays$SCT@counts, verbose = F)
DimPlot(GSE84465_tumor, reduction = 'tsne', label = T, label.size = 10)

infos1 <- Scissor(as.matrix(expr), GSE84465_tumor, phenotype, alpha = 0.05, 
                  family = "cox", Save_file = './processed/Scissor_survival.RData')

# 1 stand for Scissor+ cells associated with worse survival and 2 stand for Scissor- cells associated with good survival.
Scissor_select <- rep(0, ncol(GSE84465_tumor))
names(Scissor_select) <- colnames(GSE84465_tumor)
Scissor_select[infos1$Scissor_pos] <- 1
Scissor_select[infos1$Scissor_neg] <- 2
GSE84465_tumor <- AddMetaData(GSE84465_tumor, metadata = Scissor_select, col.name = "scissor")

identical(colnames(GSE84465_tumor),colnames(tumor))
tumor<- AddMetaData(tumor, metadata = Scissor_select, col.name = "scissor")

pdf(file="./result/scRNA/GSE84465_Scissor_tsne.pdf",bg="white",width=6,height=5.5,pointsize=8)
DimPlot(tumor, reduction = 'tsne', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
dev.off()

expr<-tumor@assays[["SCT"]]@data
expr<-as.data.frame(t(as.matrix(expr)))
expr<-dplyr::select(expr,gene)

meta<-tumor@meta.data
identical(rownames(meta),rownames(expr))
meta$gene<-expr[,1]
meta<-meta[meta$scissor !=0,]
meta$scissor<-as.character(meta$scissor)
meta[meta$scissor=="1","scissor"]<-"scissor+ cells"
meta[meta$scissor=="2","scissor"]<-"scissor- cells"

pdf(file="./result/scRNA/GSE84465_Scissor_expr.pdf",bg="white",width=3,height=3.8,pointsize=8)
ggstripchart(meta,x="scissor",y="gene", color="scissor",fill="scissor",
             palette="nejm",shape=21,size=2,
             add=c("boxplot"),# one of “none”, “reg.line”, “loess”
             add.params=list(color="black",fill="NA"),
             error.plot="errorbar",
             position=position_jitter(0.2),ggtheme=theme_classic())+
  labs(title="",x="", y = gene)+
  stat_compare_means(aes(group=scissor),method = "wilcox.test",hide.ns = F,label = "p.signif",label.x = 1.5)+
  theme(legend.position = "none")
dev.off()

###------------pathway enrich ------------
# GSE84465_tumor <- readRDS("~/bioinfo_mill/tspan7_proj/processed/GSE84465_tumor.rds")
# Marker <- readRDS("~/bioinfo_mill/tspan7_proj/processed/gene_list.rds")
# 
# expr<-GSE84465_tumor@assays[["SCT"]]@data
# expr<-as.data.frame(t(as.matrix(expr)))
# expr<-dplyr::select(expr,"TSPAN7")
# 
# GSE84465_tumor<-AddModuleScore(GSE84465_tumor,Marker,name=names(Marker))
# tmp<-GSE84465_tumor@meta.data
# colnames(tmp)[20:25]<-names(Marker)
# tmp$TSPAN7_group<-ifelse(tmp$TSPAN7_group=="TSPAN7_high","High","Low")
# tmp$TSPAN7_group<-factor(tmp$TSPAN7_group,levels = c("High","Low"))
# tmp[,20:25]<-apply(tmp[,20:25], 2, min.max.norm)
# identical(rownames(tmp),rownames(expr))
# tmp$TSPAN7<-expr$TSPAN7
# 
# library(AUCell)
# cells_rankings <- AUCell_buildRankings(GSE84465_tumor@assays$SCT@data) 
# cells_AUC <- AUCell_calcAUC(Marker, cells_rankings)
# res<-GSE84465_tumor@meta.data
# for (i in cells_AUC@NAMES) {
#   aucs <- as.numeric(getAUC(cells_AUC)[i,])
#   res$AUC<-aucs
#   colnames(res)[colnames(res)=="AUC"]<-i
# }
# res[,20:25]<-apply(res[,20:25], 2, min.max.norm)
# res$TSPAN7_group<-ifelse(res$TSPAN7_group=="TSPAN7_high","High","Low")
# res$TSPAN7_group<-factor(res$TSPAN7_group,levels = c("High","Low"))
# ggviolin(res, x = "TSPAN7_group",
#          y = c(names(Marker)),
#          combine = TRUE,ncol=3, 
#          color = "TSPAN7_group", palette = "nejm",
#          ylab = "Enrichment score", xlab =F,
#          add = "boxplot",    
#          add.params = list(binwidth = 0.1, dotsize = 0.3))+
#   theme(legend.position = "none",axis.text.x = element_text(angle = 0, hjust = 0.5))+
#   stat_compare_means(method="wilcox.test",label.x.npc= 'center',label = "p.signif")

###------------four subtype------------
source("./code/foursubtypes_gbm_function.R")
exprSet <- GSE84465_tumor@assays$SCT@data
exprSet <-as.data.frame(exprSet)
exprSet$bins <- sample(1:30, nrow(exprSet), replace=T)

geneset <- fread("~/scRNA_GBM_integration/data/IDHwt.GBM.MetaModules.tsv",data.table = F)
geneset = as.data.frame(geneset)

results <- calcSC(exprSet, geneset)
res<-results
res<-as.data.table(res)
res<-res[, c("act","random") := NULL]
res<-res %>% tidyr::pivot_wider(names_from = "geneset", values_from = "score")
res<-as.data.frame(res)
res$MESlike<-apply(res[,2:3],1,max)
res$NPClike<-apply(res[,6:7],1,max)
saveRDS(res,file = "./processed/foursubtypes_gbm_score.rds")

res2<-select(res,"sample","NPClike","OPClike","AClike","MESlike")
row.names(res2)<-res2$sample
res2<-res2[,-1]
res2$subtype<-apply(res2, 1, function(t) {colnames(res2)[which.max(t)]})
table(res2$subtype)

gene_expr<-as.data.frame(t(exprSet[gene,1:1021]))
identical(rownames(gene_expr),rownames(res2))
gene_expr<-cbind(gene_expr,res2)

# my_comparisons <- list(c("AClike","MESlike"),c("OPClike","MESlike"),c("MESlike","NPClike"))
pdf(file="./result/scRNA/foursubtypes_gbm_expr.pdf",bg="white",width=4,height=3.5,pointsize=8)
ggplot(gene_expr , aes(x=subtype, TSPAN7)) + 
  geom_jitter(aes(x=subtype, TSPAN7,colour = subtype), size = 2.5, width = 0.1, height = 0,alpha=0.5)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +     # 误差棒
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, size = 0.2)+  # 指代均数的水平横线
  theme_classic()+
  ggsci::scale_color_jama()+
  stat_compare_means(aes(group = subtype),method = "kruskal.test",hide.ns = F,label = "p.signif",label.x = 2.5) +
  labs(title = "")+
  ylab("TSPAN7 expression") +
  xlab("") +
  #ylim(35,37.5)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #legend.position = "none",
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),legend.position = "none")
dev.off()




