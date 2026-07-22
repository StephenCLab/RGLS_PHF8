####load packages####
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidydr)

# geo:GSE253901

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


load('TCELL.RData')
TCELL$group2 <- 'RGLSH'
TCELL$group2[TCELL$group== "model"] <- "EB"
TCELL$group2[TCELL$group== "contrast"] <- "Con"

# anno umap
DimPlot(TCELL,group.by = "sub_clusters",label=T)+ theme_dr()+ theme(panel.grid=element_blank())+ NoLegend()
# subset DPres
res <- subset(TCELL,subset = new2=='DPres')
# V(D)J gene expression
DotPlot(res, features = c('Rag1','Rag2'),group.by = 'group2',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +
        guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+
        labs(x='',y='')+ 
        scale_size_continuous(breaks = c(0, 25, 50, 75, 100))
# celltype percent
count_table <- as.data.frame(table(TCELL@meta.data$group2, TCELL@meta.data$new2))
names(count_table) <- c("sample", "subtype", "count")
count_table$Freq <- ave(count_table$count, count_table$sample, FUN = function(x) x / sum(x))
library(openxlsx)
write.xlsx(count_table,'TCELL_Ratio.xlsx')
count_table$sample=factor(count_table$sample,levels=c('Con','EB','RGLSH'))
p1=ggplot(count_table, aes(x = sample, y = Freq, fill = subtype))+
  geom_col()+
  theme_classic()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+
  scale_fill_manual(values = use_colors[c(4,3,2,1,6,5,11,8)])+labs(x='',y='Percent(%)')

# DPres percent
DPres_data <- subset(count_table, subtype == 'DPres')
p_inset <- ggplot(DPres_data, aes(x = sample, y = Freq)) +
  geom_col(fill = "#E64B35", width = 0.6) + 
  lim = c(0, 0.8)) +
  labs(x = NULL, y = "DPres(%)", title = "") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7),axis.title.y = element_text(size=8))
p1+ inset_element(p_inset, 
                       left = 0.7,   
                       bottom = 0.5, 
                       right = 1,   
                       top = 0.85,    
                       align_to = "panel",
                       on_top = TRUE) 

# split data
resEC=subset(res,subset=group2 %in% c('EB','Con'))
resER=subset(res,subset=group2 %in% c('EB','RGLSH'))

# edgeR  Con vs EB 
library(edgeR)
data1=resEC@assays$RNA@data
dd_CM<- as.data.frame(as.matrix(data1))
group2<-  factor(resEC@meta.data$group2,levels = c( "EB","Con"))#Con vs EB
table(group2)
dgelist <- DGEList(counts =dd_CM, group = group2)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
design <- model.matrix(~group2)
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
sample1='Con'
sample2='EB'
# save deg results
write.table(lrt, sprintf('%svs%s.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)

library(openxlsx)
library(ggplot2)
# add change result
d1=read.table(sprintf('%svs%s.txt',sample1,sample2))
colnames(d1)
d1$gene<- rownames(d1)
cut_off_FDR = 0.05  
cut_off_logFC = 0.1        
d1$change<- 'Stable'
d1$change = ifelse(d1$FDR< cut_off_FDR & abs(d1$logFC) > cut_off_logFC, 
                        ifelse(d1$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(d1)
table(d1$change)
write.table(d1, sprintf('%svs%s_d1.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)

# add label col
d1<- read.xlsx(sprintf('%svs%s.xlsx',sample1,sample2))
p <- ggplot(
  d1, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  labs(x="log2(FC)",
       y="-log10(FDR)")+
  scale_x_continuous(limits = c(-0.3, 0.3))+
  scale_y_continuous(limits = c(0,20))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('Con v.s. EB') +
  
  geom_text_repel(
    data = d1,
    aes(label = label),
    box.padding = 0.5,
    segment.color = "grey50",
    max.overlaps = 20,color='black'
  )
# p
ggsave(filename = "diffgene_ConvsEB.pdf", height = 4, width = 6, plot = p)


#  RGLSH vs EB
library(edgeR)
data2=resER@assays$RNA@data
dd_CM<- as.data.frame(as.matrix(data2))
group2<-  factor(resER@meta.data$group2,levels = c( "EB","RGLSH"))
table(group2)
dgelist <- DGEList(counts =dd_CM, group = group2)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
design <- model.matrix(~group2)
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
sample1='RGLSH'
sample2='EB'
write.table(lrt, sprintf('%svs%s.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)

library(openxlsx)
library(ggplot2)
d1=read.table(sprintf('%svs%s.txt',sample1,sample2))
colnames(d1)
d1$gene<- rownames(d1)

cut_off_FDR = 0.05  
cut_off_logFC = 0.1       
d1$change<- 'Stable'
d1$change = ifelse(d1$FDR< cut_off_FDR & abs(d1$logFC) > cut_off_logFC, 
                        ifelse(d1$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(d1)
table(d1$change)
write.table(d1, sprintf('%svs%s_d1.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)

d1<- read.xlsx(sprintf('%svs%s.xlsx',sample1,sample2))
d1$FDR <- as.numeric(d1$FDR)

p <- ggplot(
  d1, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  labs(x="log2(FC)",
       y="-log10(FDR)")+
  scale_x_continuous(limits = c(-0.3, 0.3))+
  scale_y_continuous(limits = c(0,20))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('RGLSH v.s. EB') +
  
  geom_text_repel(
    data = d1,
    aes(label = label),
    box.padding = 0.5,
    segment.color = "grey50",
    max.overlaps = 20,color='black'
  )
# p
ggsave(filename = "diffgene_RGLSHvsEB.pdf", height = 4, width = 6, plot = p)


# sup1&F1
VlnPlot(TCELL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident") & 
    theme(plot.title = element_text(size=10))

DimPlot(TCELL,group.by = "new2",label=T,split.by='orig.ident',ncol=3,cols = use_colors[c(4,3,2,1,6,5,11,8)])+ theme_dr()+ theme(panel.grid=element_blank())+ NoLegend()



DP$new5 <- DP$new4
DP$new5[DP$new5 == "DPsel-7"] <- "NKT-like"
DimPlot(DP,group.by = "orig.ident",cols = use_colors)+ theme_dr()+ theme(panel.grid=element_blank())
DimPlot(DP,group.by = "new5",label=T,split.by='orig.ident',ncol=2,cols = use_colors)+ theme_dr()+ theme(panel.grid=element_blank())+ NoLegend()

count_table <- as.data.frame(table(DP@meta.data$group, DP@meta.data$new5))
names(count_table) <- c("sample", "subtype", "count")

count_table$Freq <- ave(count_table$count, count_table$sample, FUN = function(x) x / sum(x)*100)
library(openxlsx)
write.xlsx(count_table,'DP_Ratio.xlsx')

ggplot(count_table, aes(x = sample, y = Freq, fill = subtype))+
    geom_col()+
    theme_classic()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+
    scale_fill_manual(values = use_colors[c(4,3,2,1,6,5,11,8)])+labs(x='',y='Percent(%)')


count_table <- as.data.frame(table(TCELL@meta.data$orig.ident, TCELL@meta.data$new2))
names(count_table) <- c("sample", "subtype", "count")

count_table$Freq <- ave(count_table$count, count_table$sample, FUN = function(x) x / sum(x)*100)
library(openxlsx)
write.xlsx(count_table,'TCELL_Ratio.xlsx')

ggplot(count_table, aes(x = sample, y = Freq, fill = subtype))+
    geom_col()+
    theme_classic()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+
    scale_fill_manual(values = use_colors[c(4,3,2,1,6,5,11,8)])+labs(x='',y='Percent(%)')


DotPlot(DP, features = c('Mki67','Rag1','Itm2a','Cd24a','Cd44'),group.by = 'new4',
        cols = c("lightgrey", "#E64B35FF"),
        dot.scale = 6,
        scale = TRUE) +coord_flip()+
  guides(color = guide_colorbar(title = "Avg"), size = guide_legend(title = "Per"))+labs(x='',y='')

FeaturePlot(DP, features = c('Mki67','Rag1','Itm2a','Cd24a','Cd44'),ncol=3)

pseudo_counts <- AggregateExpression(
  DP,
  assays = "RNA",
  slot = "counts", 
  group.by = "pseudobulk_key",
  return.seurat = FALSE
)$RNA

pseudo_counts <- round(pseudo_counts)


dim(pseudo_counts) 
head(colnames(pseudo_counts))

library(dplyr)
library(tidyr)

coldata <- data.frame(
  sample_key = colnames(pseudo_counts),
  stringsAsFactors = FALSE
) %>%
  separate(sample_key, into = c("orig.ident1", "orig.ident2","new5"), sep = "_")

coldata$orig.ident=paste0(coldata$orig.ident1,'_',coldata$orig.ident2)
coldata$orig.ident <- factor(coldata$orig.ident)
coldata$new5 <- factor(coldata$new5)

coldata$group <- 'RGLSH'
coldata$group[coldata$orig.ident %in% c("model_1","model_2")] <- "EB"
coldata$group[coldata$orig.ident %in% c("contrast_1","contrast_2")] <- "Con"
coldata$group <- factor(coldata$group)

head(coldata)

library(DESeq2)

celltypes <- unique(coldata$new5)

all_results1 <- list()
wb1 <- createWorkbook()
for (ct in celltypes) {
  keep <- coldata$new5 == ct
  if (sum(keep) < 2) {
    next}
  counts_sub <- pseudo_counts[, keep, drop = FALSE]
  coldata_sub <- coldata[keep, ]
  rownames(coldata_sub)=paste0(coldata_sub$orig.ident,'_',coldata_sub$new5)
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = coldata_sub,
    design = ~ group  
  )

  keep_genes <- rowSums(counts(dds)) >= 10
  dds <- dds[keep_genes, ]

#   dds <- DESeq(dds)
  dds <- DESeq(dds, test = "LRT", reduced = ~ 1)

#   res <- results(dds, contrast = c("group", "Treat", "Ctrl"))
    res <- results(dds, contrast = c("group", "Con", "EB"))
    
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(new5 = ct) %>%
    arrange(padj)
  res_df=na.omit(res_df)
  cut_off_padj = 0.05 
  cut_off_avg_log2FC = 0.5 
  res_df$change<- 'Stable'
  res_df$change = ifelse(res_df$padj< cut_off_padj & abs(res_df$avg_log2FC) > cut_off_avg_log2FC, 
                            ifelse(res_df$avg_log2FC> 0 ,'Up','Down'),
                            'Stable')
    head(res_df)
    table(res_df$change)
    write.table(res_df, sprintf('%sConvsEB_res_df.txt',ct), sep = '\t', col.names = NA, quote = FALSE)
    p <- ggplot(
    res_df, aes(x = avg_log2FC, y = -log10(padj), colour=change)) +
    geom_point(alpha=0.4, size=2) +
    scale_color_manual(values=c("Up"="#ff4757", "Down"="#546de5", "Stable"="#d2dae2"))+
    geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(cut_off_padj),lty=4,col="black",lwd=0.8) +
    labs(x="log2(FC)",
        y="-log10(padj)")+
    theme_bw()+

    theme(plot.title = element_text(hjust = 0.5), 
            legend.position="right", 
            legend.title = element_blank())+
    ggtitle('Con v.s. EB') 
    
    ggsave(filename = sprintf("%s_diffgene_ConvsEB.pdf",ct), height = 4, width = 6, plot = p)

  all_results1[[ct]] <- res_df
  sheet_name <- substr(ct, 1, 31)
  addWorksheet(wb1, sheetName = sheet_name)
  writeData(wb1, sheet = sheet_name, x =res_df)
}

output_file <- "Pseudobulk_DESeq2_Results_ConvsEB.xlsx"
saveWorkbook(wb1, file = output_file, overwrite = TRUE)


for (ct in celltypes) {
    res_df <- all_results1[[ct]]

    upgene <- paste(ct, "upgene", 'ConvsEB',sep = "_")
    assign(upgene, res_df[res_df$change=='Up',]$gene)
    print(upgene)
    upgenes <- get(upgene)
    if (length(upgenes) > 0) {
        gs = bitr(upgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
        ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE) 
        ego.bp@result <- ego.bp@result %>%
        arrange(desc(Count))
        write.csv(ego.bp@result[ego.bp@result$pvalue<0.05,],
                file = sprintf('%s_GOBP_up.csv',upgene))

        pdf(file=sprintf("%s_GOBP_up.pdf",upgene),width=6,height=6)
        print(dotplot(ego.bp, showCategory=10,,title=sprintf("%s_GOBP_up",upgene),orderBy='Count'))
        dev.off()
    }
    
    downgene <- paste(ct, "downgene", 'ConvsEB',sep = "_")
    assign(downgene, res_df[res_df$change=='Down',]$gene)
    print(downgene)
    downgenes <- get(downgene)
    if (length(downgenes) > 0) {
        gs = bitr(downgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
        ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE) 
        ego.bp@result <- ego.bp@result %>%
        arrange(desc(Count))
        write.csv(ego.bp@result[ego.bp@result$pvalue<0.05,],
                file = sprintf('%s_GOBP_down.csv',downgene))

        pdf(file=sprintf("%s_GOBP_down.pdf",downgene),width=6,height=6)
        print(dotplot(ego.bp, showCategory=10,,title=sprintf("%s_GOBP_down",downgene),orderBy='Count'))
        dev.off()
    }
}


all_results2 <- list()
wb2 <- createWorkbook()
for (ct in celltypes) {

  keep <- coldata$new5 == ct
  if (sum(keep) < 2) {
    next
  }
  counts_sub <- pseudo_counts[, keep, drop = FALSE]
  coldata_sub <- coldata[keep, ]
  rownames(coldata_sub)=paste0(coldata_sub$orig.ident,'_',coldata_sub$new5)
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = coldata_sub,
    design = ~ group  
  )

  keep_genes <- rowSums(counts(dds)) >= 10
  dds <- dds[keep_genes, ]
#   dds <- DESeq(dds)
  dds <- DESeq(dds, test = "LRT", reduced = ~ 1)

#   res <- results(dds, contrast = c("group", "Treat", "Ctrl"))
    res <- results(dds, contrast = c("group", "RGLSH", "EB"))

  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(new5 = ct) %>%
    arrange(padj)
    res_df=na.omit(res_df)
      cut_off_padj = 0.05
    cut_off_avg_log2FC = 0.5   
    res_df$change<- 'Stable'
    res_df$change = ifelse(res_df$padj< cut_off_padj & abs(res_df$avg_log2FC) > cut_off_avg_log2FC, 
                            ifelse(res_df$avg_log2FC> 0 ,'Up','Down'),
                            'Stable')
    head(res_df)
    table(res_df$change)
    write.table(res_df, sprintf('%sRGLSHvsEB_res_df.txt',ct), sep = '\t', col.names = NA, quote = FALSE)

    p <- ggplot(
    res_df, aes(x = avg_log2FC, y = -log10(padj), colour=change)) +
    geom_point(alpha=0.4, size=2) +
    scale_color_manual(values=c("Up"="#ff4757", "Down"="#546de5", "Stable"="#d2dae2"))+
    geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(cut_off_padj),lty=4,col="black",lwd=0.8) +
    labs(x="log2(FC)",
        y="-log10(padj)")+
    # scale_x_continuous(limits = c(-0.3, 0.3))+
    # scale_y_continuous(limits = c(0,20))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5), 
            legend.position="right", 
            legend.title = element_blank())+
    ggtitle('RGLSH v.s. EB') 
    
    ggsave(filename = sprintf("%s_diffgene_RGLSHvsEB.pdf",ct), height = 4, width = 6, plot = p)


  all_results2[[ct]] <- res_df
  sheet_name <- substr(ct, 1, 31)
  addWorksheet(wb2, sheetName = sheet_name)
  writeData(wb2, sheet = sheet_name, x =res_df)
}

output_file <- "Pseudobulk_DESeq2_Results_RGLSHvsEB.xlsx"
saveWorkbook(wb2, file = output_file, overwrite = TRUE)


for (ct in celltypes) {
    res_df <- all_results2[[ct]]

    upgene <- paste(ct, "upgene", 'RGLSHvsEB',sep = "_")
    assign(upgene, res_df[res_df$change=='Up',]$gene)
    print(upgene)
    upgenes <- get(upgene)
    if (length(upgenes) > 0) {
        gs = bitr(upgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
        ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE) 
        ego.bp@result <- ego.bp@result %>%
        arrange(desc(Count))
        write.csv(ego.bp@result[ego.bp@result$pvalue<0.05,],
                file = sprintf('%s_GOBP_up.csv',upgene))

        pdf(file=sprintf("%s_GOBP_up.pdf",upgene),width=6,height=6)
        print(dotplot(ego.bp, showCategory=10,,title=sprintf("%s_GOBP_up",upgene),orderBy='Count'))
        dev.off()
    }
    
    downgene <- paste(ct, "downgene", 'RGLSHvsEB',sep = "_")
    assign(downgene, res_df[res_df$change=='Down',]$gene)
    print(downgene)
    downgenes <- get(downgene)
    if (length(downgenes) > 0) {
        gs = bitr(downgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
        ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Mm.eg.db,ont= "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE) 
        ego.bp@result <- ego.bp@result %>%
        arrange(desc(Count))
        write.csv(ego.bp@result[ego.bp@result$pvalue<0.05,],
                file = sprintf('%s_GOBP_down.csv',downgene))

        pdf(file=sprintf("%s_GOBP_down.pdf",downgene),width=6,height=6)
        print(dotplot(ego.bp, showCategory=10,,title=sprintf("%s_GOBP_down",downgene),orderBy='Count'))
        dev.off()
    }
}







