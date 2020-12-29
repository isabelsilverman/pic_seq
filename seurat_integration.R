library(Seurat)

kidney64 = readRDS("biopsy64/kidneyB64_final.rds")
kidney64[["celltype"]] <- kidney64@active.ident
kidney68 = readRDS("biopsy68/kidneyB68_final.rds")
kidney68[["celltype"]] <- kidney68@active.ident
kidney.list <- c(kidney64, kidney68)

kidney.anchors <- FindIntegrationAnchors(object.list = kidney.list, dims = 1:30)
kidney.integrated <- IntegrateData(anchorset = kidney.anchors, dims = 1:30)

library(ggplot2)
library(cowplot)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(kidney.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
kidney.integrated <- ScaleData(kidney.integrated, verbose = FALSE)
kidney.integrated <- RunPCA(kidney.integrated, npcs = 30, verbose = FALSE)
kidney.integrated <- FindNeighbors(kidney.integrated, dims = 1:9)
kidney.integrated <- FindClusters(kidney.integrated, resolution = 1.0)
kidney.integrated <- RunUMAP(kidney.integrated, reduction = "pca", dims = 1:9)
DimPlot(kidney.integrated, reduction = "umap")

#label clusters
kidney.markers <- FindAllMarkers(kidney.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

kidney.genes <- c('NPHS1','NPHS2','SLC5A2','SLC22A6','AGT','SLC34A1','LRP2', 'CUBN','SLC14A2',
                  'CLCNKA','SLC12A1','UMOD','OXTR','NOS1','SLC12A3','AQP2','SCNN1G',
                  'CALB1','SLC4A1','SLC26A4','PECAM1','VWF','VCAM1','ICAM1')
immune.genes <- c('IL7R','CCR7','S100A4','CD14','LYZ','MS4A1','CD8A','FCGR3A','MS4A7',
                  'GNLY','NKG7','FCER1A','CST3','CD58','CD3E','CCL5','CD68','C1QA','C1QB','C1QC')
for(gene in kidney.genes) {
  print(kidney.markers[kidney.markers$gene == gene,])
}
for(gene in immune.genes) {
  print(kidney.markers[kidney.markers$gene == gene,])
}
new.cluster.ids <- c('EC','1','EC','LOH','LOH','5','PC','EC','Macrophage','9',
                     'IC','Lymphocyte','LOH','PT','EC + T','PT','LOH + T','IC','18','PC + IC',
                     'PT + EC','PC','22','EC + Macrophage','LOH + Macrophage','PC','IC','LOH + PC','LOH + IC','IC + Macrophage')

names(new.cluster.ids) <- levels(kidney.integrated)
kidney.integrated <- RenameIdents(kidney.integrated, new.cluster.ids)
DimPlot(kidney.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(kidney.integrated, file = "kidney_final.rds")




