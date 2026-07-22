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

####sobj2####
sobj=readRDS('sobj_v4.0427.rds')
sobj2 <- subset(sobj,subset=sample %in% c('WT+Model+OverExpress-empty','WT+Model+OverExpress'))
sobj2$sample <- factor(sobj2$sample,levels = c('WT+Model+OverExpress-empty','WT+Model+OverExpress'))
saveRDS(sobj2,'sobj2.rds')

sobj2 <- NormalizeData(sobj2, verbose = T)
sobj2 <- ScaleData(sobj2, verbose = T)
sobj2<- FindVariableFeatures(sobj2,nfeatures = 2000)
sobj2 <- RunPCA(sobj2, npcs = 30, verbose = T)
sobj2 <- RunUMAP(sobj2, reduction = "pca", dims = 1:30)
sobj2 <- FindNeighbors(sobj2, reduction = "pca", dims = 1:30)
sobj2 <- FindClusters(sobj2,resolution = 0.8)
saveRDS(sobj2,'sobj2_cluster.rds')


sobj2=readRDS('sobj2_cluster.rds')
DimPlot(sobj2,cols = use_colors,label = TRUE) + theme_dr()+ theme(panel.grid=element_blank())+ NoLegend()
DimPlot(sobj2,group.by = 'sample',cols = use_colors,label = FALSE) + theme_dr()+ theme(panel.grid=element_blank())
DimPlot(sobj2,group.by = 'celltype',cols = use_colors,label = TRUE) + theme_dr()+ theme(panel.grid=element_blank())

markers=c('Cd3e','Cd3d','Cd4','Cd8a','Cd8b1',
          'Cd79a','Cd19','Ms4a1',
          'Krt18','Epcam','Foxp1',
          'Cd68','C1qa','C1qb',
          'Dcn','Pdgfra','Col1a1')
DotPlot(sobj2, features = markers,group.by = 'seurat_clusters',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')

DotPlot(sobj2, features = markers,,group.by = 'celltype',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')

for (gene in markers) {
  #  FeaturePlot
  p <- FeaturePlot(sobj2, features = gene, reduction = "umap") +
    ggtitle(gene) + theme_dr()+ theme(panel.grid=element_blank())  
  #  PDF
  pdf(file = paste0(gene, "_FeaturePlot_sobj2.pdf"), width = 6, height = 4)
  print(p)
  dev.off()
}

#sobj2
library(ggplot2)
library(dplyr)

count_table <- as.data.frame(table(sobj2@meta.data$sample, sobj2@meta.data$celltype))
names(count_table) <- c("sample", "celltype", "count")
count_table$Freq <- ave(count_table$count, count_table$sample, FUN = function(x) x / sum(x))
library(openxlsx)
write.xlsx(count_table,'Ratio.xlsx')


p1=ggplot(count_table, aes(x = sample, y = Freq, fill = celltype))+
  geom_col()+
  theme_classic()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+
  scale_fill_manual(values = use_colors)+labs(x='',y='Percent(%)')
ggsave(filename = "p_Ratio-type.pdf", plot = p1, height = 4, width = 6)

####T_sobj2####
T_sobj=readRDS('T_sobj_v4.0427.rds')
T_sobj2 <- subset(T_sobj,subset=sample %in% c('WT+Model+OverExpress-empty','WT+Model+OverExpress'))
T_sobj2$sample <- factor(T_sobj2$sample,levels = c('WT+Model+OverExpress-empty','WT+Model+OverExpress'))
saveRDS(T_sobj2,'T_sobj2.rds')

T_sobj2 <- NormalizeData(T_sobj2, verbose = T)
T_sobj2 <- ScaleData(T_sobj2, verbose = T)
T_sobj2<- FindVariableFeatures(T_sobj2,nfeatures = 2000)
T_sobj2 <- RunPCA(T_sobj2, npcs = 30, verbose = T)
# t-SNE and Clustering
T_sobj2 <- RunUMAP(T_sobj2, reduction = "pca", dims = 1:30)
T_sobj2 <- FindNeighbors(T_sobj2, reduction = "pca", dims = 1:30)
T_sobj2 <- FindClusters(T_sobj2,resolution = 0.7)
saveRDS(T_sobj2,'T_sobj2_cluster.rds')

T_sobj2=readRDS('T_sobj2_cluster.rds')
T_markers <- c('Cd3e','Cd3d',
              'Cd4',
              'Cd8a','Cd8b1',
              'Cd69','Rag1',
              'Foxp3','Il2ra','Ctla4',
              'Klrk1','Klrc2')

for (gene in T_markers) {
  #  FeaturePlot
  p <- FeaturePlot(T_sobj2, features = gene, reduction = "umap") +
    ggtitle(gene) + theme_dr()+ theme(panel.grid=element_blank())  
  #  PDF
  pdf(file = paste0(gene, "_FeaturePlot_T_sobj2.pdf"), width = 6, height = 4)
  print(p)
  dev.off()
}
# T_obj4$subtype <- factor(T_obj4$subtype,levels=c('DN','DP','SP','ISP','Treg','NKT'))
# FP_NKT_T_sobj2.pdf
FeaturePlot(T_sobj2, features = c("Klrk1")) + theme_dr()+ theme(panel.grid=element_blank())
FeaturePlot(T_sobj2, features = c("Foxp3") )+ theme_dr()+ theme(panel.grid=element_blank())
FeaturePlot(T_sobj2, features = c("Cd69", "Rag1") , blend = TRUE,cols=c("lightgrey","#E64B35",'blue'))
FeaturePlot(T_sobj2, features = c("Cd4", "Cd8a") , blend = TRUE,cols=c("lightgrey","#E64B35",'blue'))

# DP_sc_T_sobj2.pdf
DimPlot(T_sobj2,cols = use_colors,label = TRUE) + theme_dr()+ theme(panel.grid=element_blank())+ NoLegend()
DimPlot(T_sobj2,group.by = 'sample',cols = use_colors,label = FALSE) + theme_dr()+ theme(panel.grid=element_blank())
DimPlot(T_sobj2,group.by = 'subtype',cols = use_colors[c(6,5,4,3,2,1)],label = TRUE) + theme_dr()+ theme(panel.grid=element_blank())

DotPlot(T_sobj2, features = T_markers,group.by = 'seurat_clusters',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')

DotPlot(T_sobj2, features = T_markers,,group.by = 'subtype',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')


count_table <- as.data.frame(table(T_sobj2@meta.data$sample, T_sobj2@meta.data$subtype))
names(count_table) <- c("sample", "subtype", "count")

count_table$Freq <- ave(count_table$count, count_table$sample, FUN = function(x) x / sum(x))
library(openxlsx)
write.xlsx(count_table,'T_sobj2_Ratio.xlsx')

count_table$sample2=count_table$sample
count_table$sample <- 'WMO'
count_table$sample[count_table$sample2== "WT+Model+OverExpress-empty"] <- "WMOe"
count_table$sample=factor(count_table$sample,levels=c('WMO','WMOe'))


p1=ggplot(count_table, aes(x = sample, y = Freq, fill = subtype))+
  geom_col()+
  theme_classic()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+
  scale_fill_manual(values = use_colors[c(6,5,4,3,2,1)])+labs(x='',y='Percent(%)')
# ggsave(filename = "p_Ratio-T_sobj2type.pdf", plot = p1, height = 4, width = 6)

nkt_data <- subset(count_table, subtype == 'NKT')
p_inset <- ggplot(nkt_data, aes(x = sample, y = Freq)) +
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

NKT_Sig=c('Klrk1','Klrc2','Zbtb16','Tbx21','Ifng','Klrd1','Cd1d1','Cd1d2')
gene_sets <- list(
  NKT_Sig = NKT_Sig
)

T_sobj2 <- AddModuleScore(
  object = T_sobj2,
  features = gene_sets,
  name = "NKT_Sig")

DotPlot(T_sobj2, features = 'NKT_Sig1',group.by = 'subtype',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')

T_sobj2$sample2 <- 'WMO'
T_sobj2$sample2[T_sobj2$sample== "WT+Model+OverExpress-empty"] <- "WMOe"
T_sobj2$sample2=factor(T_sobj2$sample2,levels=c('WMO','WMOe'))

DotPlot(T_sobj2, features = NKT_Sig,group.by = 'sample2',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = FALSE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')

