library(dplyr)
library(Seurat)
library(patchwork)

#load biopsy data
data_dir = "filtered_feature_bc_matrix"
kidney.data <- Read10X(data.dir = data_dir)
kidney <- CreateSeuratObject(counts = kidney.data, project = "kidney", min.cells = 3, min.features = 200)

#filter bad cells
kidney[["percent.mt"]] <- PercentageFeatureSet(kidney, pattern = "^MT-")
kidney <- subset(kidney, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

kidney <- NormalizeData(kidney, normalization.method = "LogNormalize", scale.factor = 10000)

kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(kidney)
kidney <- ScaleData(kidney, features = all.genes)

#overcluster with high res
kidney <- RunPCA(kidney, features = VariableFeatures(object = kidney))
kidney <- FindNeighbors(kidney, dims = 1:9)
kidney <- FindClusters(kidney, resolution = 1.0)

kidney <- RunUMAP(kidney, dims = 1:9)
DimPlot(kidney, reduction = "umap")

#label clusters
kidney.markers <- FindAllMarkers(kidney, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

kidney.genes <- c('NPHS1','NPHS2','SLC5A2','SLC22A6','AGT','SLC34A1','LRP2', 'CUBN','SLC14A2',
                  'CLCNKA','SLC12A1','UMOD','OXTR','NOS1','SLC12A3','AQP2','SCNN1G',
                  'CALB1','SLC4A1','SLC26A4','PECAM1','VWF')
immune.genes <- c('IL7R','CCR7','S100A4','CD14','LYZ','MS4A1','CD8A','FCGR3A','MS4A7',
                  'GNLY','NKG7','FCER1A','CST3','CD58','CD3E','CCL5','CD68','C1QA','C1QB','C1QC')
for(gene in kidney.genes) {
  print(kidney.markers[kidney.markers$gene == gene,])
}
for(gene in immune.genes) {
  print(kidney.markers[kidney.markers$gene == gene,])
}
#68
# new.cluster.ids <- c('EC','EC','PT','LOH','4','LOH','T+EC','PT','Macrophage','PC',
#                     'IC','LOH+T','IC','LOH','Lymphocytes','EC','16','DC','PT+LOH','PT+NK',
#                     'PC', 'T+EC','PT+T','IC+NK','PC+EC','EC+Macrophage','IC','PT+IC')

#64
#new.cluster.ids <- c('EC','EC','PC','3','LOH','5','6','EC','LOH','IC',
#                     'Macrophages','EC','LOH','LOH','Lymphocytes','PC + IC','16','LOH + Macrophage','PC','Macrophage',
#                     'EC + T','LOH + T','LOH + PC','LOH + PC','IC','EC + Macrophage','LOH + EC','LOH + IC', 'IC + Macrophage')

#65
#new.cluster.ids <- c('Macrophages','1','2','immune','EC','EC','6','LOH','immune','EC',
#                     'PT','Lymphocytes','PT','EC','immune','IC','16','DC','18','Lymphocytes(B?)')
  
names(new.cluster.ids) <- levels(kidney)
kidney <- RenameIdents(kidney, new.cluster.ids)
DimPlot(kidney, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(kidney, file = "kidney_final.rds")
