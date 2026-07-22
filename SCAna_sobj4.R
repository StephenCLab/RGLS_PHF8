use_colors <- c(
  "#E64B35",  
  "#4DBBD5", 
  "#3C5488", 
  "#F39B7F", 
  "#8491B4", 
  "#91D1C2", 
  "#DC0000", 
  "#7E6148", 
  "#B09C85",  
  "#F7B6B2",  
  "#B09CD6",  
  "#B3B3B3",  
  "#FFDC91", 
  "#A6D854", 
  "#D95F02",  
  "#7570B3", 
  "#E7298A", 
  "#E6AB02", 
  "#00A087", 
  "#66A61E" 
)


library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidydr)

####sobj4####
sobj=readRDS('sobj_v4.0427.rds')
sobj4 <- subset(sobj,subset=sample %in% c('WT+Model','WT+Model+Treatment','KO+Model+Treatment','KO+Model+Treatment+OverExpress'))
sobj4$sample <- factor(sobj4$sample,levels = c('WT+Model','WT+Model+Treatment','KO+Model+Treatment','KO+Model+Treatment+OverExpress'))
saveRDS(sobj4,'sobj4.rds')

sobj4 <- NormalizeData(sobj4, verbose = T)
sobj4 <- ScaleData(sobj4, verbose = T)
sobj4<- FindVariableFeatures(sobj4,nfeatures = 2000)
sobj4 <- RunPCA(sobj4, npcs = 30, verbose = T)
sobj4 <- RunUMAP(sobj4, reduction = "pca", dims = 1:30)
sobj4 <- FindNeighbors(sobj4, reduction = "pca", dims = 1:30)
sobj4 <- FindClusters(sobj4,resolution = 0.8)
saveRDS(sobj4,'sobj4_cluster.rds')

DimPlot(sobj4,cols = use_colors,label = TRUE) + theme_dr()+ theme(panel.grid=element_blank())+ NoLegend()
DimPlot(sobj4,group.by = 'sample',cols = use_colors,label = FALSE) + theme_dr()+ theme(panel.grid=element_blank())
DimPlot(sobj4,group.by = 'celltype',cols = use_colors,label = TRUE) + theme_dr()+ theme(panel.grid=element_blank())

markers=c('Cd3e','Cd3d','Cd4','Cd8a','Cd8b1',
          'Cd79a','Cd19','Ms4a1',
          'Krt18','Epcam','Foxp1',
          'Cd68','C1qa','C1qb',
          'Dcn','Pdgfra','Col1a1')
DotPlot(sobj4, features = markers,group.by = 'seurat_clusters',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')

DotPlot(sobj4, features = markers,,group.by = 'celltype',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')


for (gene in markers) {
  #  FeaturePlot
  p <- FeaturePlot(sobj4, features = gene, reduction = "umap") +
    ggtitle(gene) + theme_dr()+ theme(panel.grid=element_blank())  
  #  PDF
  pdf(file = paste0(gene, "_FeaturePlot_sobj4.pdf"), width = 6, height = 4)
  print(p)
  dev.off()
}

#sobj4
library(ggplot2)
library(dplyr)

count_table <- as.data.frame(table(sobj4@meta.data$sample, sobj4@meta.data$celltype))
names(count_table) <- c("sample", "celltype", "count")

count_table$Freq <- ave(count_table$count, count_table$sample, FUN = function(x) x / sum(x))
library(openxlsx)
write.xlsx(count_table,'Ratio.xlsx')

# NO label
p1=ggplot(count_table, aes(x = sample, y = Freq, fill = celltype))+
  geom_col()+
  theme_classic()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+
  scale_fill_manual(values = use_colors)+labs(x='',y='Percent(%)')
ggsave(filename = "p_Ratio-type.pdf", plot = p1, height = 4, width = 6)

####T_sobj4####
T_sobj=readRDS('T_sobj_v4.0427.rds')
T_sobj4 <- subset(T_sobj,subset=sample %in% c('WT+Model','WT+Model+Treatment','KO+Model+Treatment','KO+Model+Treatment+OverExpress'))
T_sobj4$sample <- factor(T_sobj4$sample,levels = c('WT+Model','WT+Model+Treatment','KO+Model+Treatment','KO+Model+Treatment+OverExpress'))
saveRDS(T_sobj4,'T_sobj4.rds')

T_sobj4 <- NormalizeData(T_sobj4, verbose = T)
T_sobj4 <- ScaleData(T_sobj4, verbose = T)
T_sobj4<- FindVariableFeatures(T_sobj4,nfeatures = 2000)
T_sobj4 <- RunPCA(T_sobj4, npcs = 30, verbose = T)
# t-SNE and Clustering
T_sobj4 <- RunUMAP(T_sobj4, reduction = "pca", dims = 1:30)
T_sobj4 <- FindNeighbors(T_sobj4, reduction = "pca", dims = 1:30)
T_sobj4 <- FindClusters(T_sobj4,resolution = 0.7)
saveRDS(T_sobj4,'T_sobj4_cluster.rds')

T_sobj4=readRDS('T_sobj4_cluster.rds')
T_markers <- c('Cd3e','Cd3d',
              'Cd4',
              'Cd8a','Cd8b1',
              'Cd69','Rag1',
              'Foxp3','Il2ra','Ctla4',
              'Klrk1','Klrc2')

for (gene in T_markers) {
  #  FeaturePlot
  p <- FeaturePlot(T_sobj4, features = gene, reduction = "umap") +
    ggtitle(gene) + theme_dr()+ theme(panel.grid=element_blank())  
  #  PDF
  pdf(file = paste0(gene, "_FeaturePlot_T_sobj4.pdf"), width = 6, height = 4)
  print(p)
  dev.off()
}
# T_obj4$subtype <- factor(T_obj4$subtype,levels=c('DN','DP','SP','ISP','Treg','NKT'))
# FP_NKT_T_sobj4.pdf
FeaturePlot(T_sobj4, features = c("Klrk1")) + theme_dr()+ theme(panel.grid=element_blank())
FeaturePlot(T_sobj4, features = c("Foxp3") )+ theme_dr()+ theme(panel.grid=element_blank())
FeaturePlot(T_sobj4, features = c("Cd69", "Rag1") , blend = TRUE,cols=c("lightgrey","#E64B35",'blue'))
FeaturePlot(T_sobj4, features = c("Cd4", "Cd8a") , blend = TRUE,cols=c("lightgrey","#E64B35",'blue'))

# DP_sc_T_sobj4.pdf
DimPlot(T_sobj4,cols = use_colors,label = TRUE) + theme_dr()+ theme(panel.grid=element_blank())+ NoLegend()
DimPlot(T_sobj4,group.by = 'sample',cols = use_colors,label = FALSE) + theme_dr()+ theme(panel.grid=element_blank())
DimPlot(T_sobj4,group.by = 'subtype',cols = use_colors[c(6,5,4,3,2,1)],label = TRUE) + theme_dr()+ theme(panel.grid=element_blank())

DotPlot(T_sobj4, features = T_markers,group.by = 'seurat_clusters',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')

DotPlot(T_sobj4, features = T_markers,,group.by = 'subtype',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')


count_table <- as.data.frame(table(T_sobj4@meta.data$sample, T_sobj4@meta.data$subtype))
names(count_table) <- c("sample", "subtype", "count")
# 
count_table$Freq <- ave(count_table$count, count_table$sample, FUN = function(x) x / sum(x))
library(openxlsx)
write.xlsx(count_table,'T_sobj4_Ratio.xlsx')

count_table$sample2=count_table$sample
count_table$sample <- 'WM'
count_table$sample[count_table$sample2== "WT+Model+Treatment"] <- "WMT"
count_table$sample[count_table$sample2== "KO+Model+Treatment"] <- "KMT"
count_table$sample[count_table$sample2== "KO+Model+Treatment+OverExpress"] <- "KMTO"
count_table$sample=factor(count_table$sample,levels=c('WM','WMT','KMT','KMTO'))

# NO label
p1=ggplot(count_table, aes(x = sample, y = Freq, fill = subtype))+
  geom_col()+
  theme_classic()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+
  scale_fill_manual(values = use_colors[c(6,5,4,3,2,1)])+labs(x='',y='Percent(%)')
# ggsave(filename = "p_Ratio-T_sobj4type.pdf", plot = p1, height = 4, width = 6)



nkt_data <- subset(count_table, subtype == 'NKT')
p_inset <- ggplot(nkt_data, aes(x = sample, y = Freq)) +
   #  geom_col OR geom_line + geom_point
  geom_col(fill = "#E64B35", width = 0.6) + 
  coord_cartesian(ylim = c(0, 0.1)) +
  labs(x = NULL, y = "NKT(%)", title = "") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7),axis.title.y = element_text(size=8))

p1+ inset_element(p_inset, 
                       left = 0.7,   
                       bottom = 0.5,
                       right = 1,  
                       top = 0.85,     
                       align_to = "panel",
                       on_top = TRUE)      



DimPlot(sobj4,group.by = "seurat_clusters",cols = use_colors,label=T)+ theme_dr()+ theme(panel.grid=element_blank())+ NoLegend()
DimPlot(sobj4,group.by = "celltype",label=T,split.by='sample2',ncol=2,cols = use_colors)+ theme_dr()+ theme(panel.grid=element_blank())+ NoLegend()

DimPlot(T_sobj4,group.by = "seurat_clusters",cols = use_colors,label=T)+ theme_dr()+ theme(panel.grid=element_blank())+ NoLegend()
DimPlot(T_sobj4,group.by = "sample",cols = use_colors)+ theme_dr()+ theme(panel.grid=element_blank())
DimPlot(T_sobj4,group.by = "subtype2",label=T,split.by='sample',ncol=2,cols = use_colors[c(5,4,3,2,1)])+ theme_dr()+ theme(panel.grid=element_blank())+ NoLegend()


T4_sub <- subset(T_sobj4,subset=subtype %in% c('DP','NKT'))
T4_sub <- NormalizeData(T4_sub, verbose = T)
T4_sub <- ScaleData(T4_sub, verbose = T)
T4_sub<- FindVariableFeatures(T4_sub,nfeatures = 2000)
T4_sub <- RunPCA(T4_sub, npcs = 30, verbose = T)
# t-SNE and Clustering
T4_sub <- RunUMAP(T4_sub, reduction = "pca", dims = 1:30)
T4_sub <- FindNeighbors(T4_sub, reduction = "pca", dims = 1:30)
T4_sub <- FindClusters(T4_sub,resolution = 0.7)
saveRDS(T4_sub,'T4_sub_cluster.rds')

DimPlot(T4_sub,label = TRUE) + theme_dr()+ theme(panel.grid=element_blank())+ NoLegend()

DotPlot(T4_sub, features = c('Klrk1','Klrc2','Mki67','Rag1','Itm2a'),group.by = 'seurat_clusters',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')

NKT 'Klrk1','Klrc2' 5
DP（DPblas，DP Mki67 + ）7,9
DP（DPres，DP Rag1high）0,1,2,3,6,8,10
DP（DPsels，DP Itm2a + ）4

DotPlot(T4_sub, features = c('Klrk1','Klrc2','Mki67','Rag1','Itm2a','Cd4','Cd24a','Cd44'),group.by = 'seurat_clusters',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')

DPsel/stage 0 iNKT-like（CD24a high CD44 low）4


T4_sub$anno1<- as.character(T4_sub$seurat_clusters)
table(T4_sub$anno1)
T4_sub@meta.data[T4_sub$seurat_clusters %in% c(5),]$anno1<- 'NKT'
T4_sub@meta.data[T4_sub$seurat_clusters %in% c(7,9),]$anno1<- 'DPblas'
T4_sub@meta.data[T4_sub$seurat_clusters %in% c(0,1,2,3,6,8,10),]$anno1<- 'DPres'
T4_sub@meta.data[T4_sub$seurat_clusters %in% c(4),]$anno1<- 'DPsels'
saveRDS(T4_sub,'T4_sub.rds')

DotPlot(T4_sub, features = c('Mki67','Rag1','Itm2a','Cd4','Cd24a','Cd44','Klrk1','Klrc2'),group.by = 'anno1',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')


T_markers <- c('Cd3e','Cd3d',
              'Cd4',
              'Cd8a','Cd8b1',
              'Cd69','Rag1',
              'Foxp3','Il2ra','Ctla4',
              'Klrk1','Klrc2')
DotPlot(T4_sub, features = T_markers,group.by = 'subtype2',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')

count_table <- as.data.frame(table(T4_sub@meta.data$sample2, T4_sub@meta.data$anno1))
names(count_table) <- c("sample", "subtype", "count")

count_table$Freq <- ave(count_table$count, count_table$sample, FUN = function(x) x / sum(x)*100)
library(openxlsx)
write.xlsx(count_table,'T4_sub_Ratio.xlsx')

ggplot(count_table, aes(x = sample, y = Freq, fill = subtype))+
    geom_col()+
    theme_classic()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+
    scale_fill_manual(values = use_colors[c(4,3,2,1,6,5,11,8)])+labs(x='',y='Percent(%)')

library(monocle)
sobj=readRDS('T4_sub.rds')
# monocle2 
expr_matrix <- as.data.frame(sobj@assays$RNA@counts)
#
sample_sheet <- sobj@meta.data
# 
gene_annotation <- data.frame(
  gene_short_name=rownames(sobj@assays$RNA),
  row.names = rownames(sobj@assays$RNA)
)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

#monocle
cd <- newCellDataSet(as(as.matrix(expr_matrix),"sparseMatrix"), 
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())
# cd=cd[, sample(ncol(cd), 3000)]
# expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * 0.01))

expressed_genes <- rowSums(exprs(cd) > 0) > dim(cd)[2] * 0.01  
cd <- cd[expressed_genes, ]
###size factors和dispersions
cd <- estimateSizeFactors(cd)
cd <- estimateDispersions(cd)

cd <- detectGenes(cd, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * 0.01))

length(expressed_genes)

disp_table <- dispersionTable(cd)

ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)

diff_test_res <- differentialGeneTest(cd[expressed_genes,],fullModelFormulaStr = "~anno1")

ordering_genes <- row.names (subset(diff_test_res, qval < 1e-5))
ordering_genes <- intersect(ordering_genes, expressed_genes)
cd <- setOrderingFilter(cd, ordering_genes)

cd <- reduceDimension(cd, max_components = 2, method = 'DDRTree')

cd <- orderCells(cd, reverse = F) 
head(cd@phenoData@data)
saveRDS(cd,file='T4_sub_m2_cd_ALL.rds')

pData(cd)$anno1 <- factor(pData(cd)$anno1, levels = c('DPblas','DPres','DPsels','NKT'))


pdf(file="cd.Pseudotime.pdf",width=6,height=4)
plot_cell_trajectory(cd, color_by = "Pseudotime")
dev.off()

pdf(file="cd.trajectory.pdf",width=6,height=4)
plot_cell_trajectory(cd, color_by = "anno1", cell_size = 0.5,cell_link_size = 0.5,label_cell = TRUE)
dev.off()

pdf(file="cd.trajectory_label_cell.pdf",width=12,height=8)
plot_cell_trajectory(cd, color_by = "anno1", cell_size = 0.5,cell_link_size = 0.5,label_cell = TRUE) + facet_wrap(~sample2, nrow = 2, scales = "free")+ggtitle("celltype") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()



NKT_Sig=c('Klrk1','Klrc2','Zbtb16','Tbx21','Ifng','Klrd1')
gene_sets <- list(
  NKT_Sig = NKT_Sig
)

T_sobj4 <- AddModuleScore(
  object = T_sobj4,
  features = gene_sets,
  name = "NKT_Sig"  
)

DotPlot(T_sobj4, features = 'NKT_Sig1',group.by = 'subtype2',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')


DotPlot(T_sobj4, features = NKT_Sig,group.by = 'subtype2',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')





