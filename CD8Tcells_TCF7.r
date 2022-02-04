#using R, only this version is ok
module load R/3.5.1-gccmkl
R
library(Seurat)
library(RColorBrewer) ;
mycolors3 <- brewer.pal(name="YlOrRd", n=9) ;
mycolors3 <- c("gray95", mycolors3) ;

kpkpl.data <- Read10X(data.dir = "/endosome/archive/bioinformatics/Li_lab/s186258/Huiyu_Project/Cellranger_output/AGGkpkpl/outs/filtered_feature_bc_matrix") ;
kpkpl.data <- Read10X(data.dir = "D:\\scRNAseq\\KPKPLtrial\\") ;

tmp <- colnames(kpkpl.data) ; 
tmp1 <- sapply(strsplit(tmp, "-"), function(x) x[1]) ;
tmp2 <- sapply(strsplit(tmp, "-"), function(x) x[2]) ;
library(plyr) ;
tmp3 <- revalue(tmp2, c("1"="KP9_3", "2"="KPL9_3_1")) ;

tmp4 <- paste(tmp1, tmp3, sep="-") ;
colnames(kpkpl.data) <- tmp4 ;

kpkpl <- CreateSeuratObject(counts = kpkpl.data, project = "kpkpl", min.cells = 3, min.features = 100) ;
kpkpl@meta.data$sample <- sapply(strsplit(colnames(kpkpl), "-"), function(x) x[2]);

#Human reference genome is ^MT-
kpkpl[["percent.mt"]] <- PercentageFeatureSet(kpkpl, pattern = "^mt-")
VlnPlot(kpkpl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kpkpl <- subset(kpkpl, subset = nFeature_RNA > 200 & percent.mt < 10) ; #9624 cells

#Normalizing the data
kpkpl <- NormalizeData(kpkpl, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features
kpkpl <- FindVariableFeatures(kpkpl, selection.method = "vst", nfeatures = 2000)

#Scaling the data
all.genes <- rownames(kpkpl)
kpkpl <- ScaleData(kpkpl, features = all.genes)

#Linear dimensional reduction
kpkpl <- RunPCA(kpkpl, features = VariableFeatures(object=kpkpl))

# t-SNE and Clustering
kpkpl <- RunUMAP(kpkpl, reduction = "pca", dims = 1:50)
kpkpl <- FindNeighbors(kpkpl, reduction = "pca", dims = 1:50)
kpkpl <- FindClusters(kpkpl, resolution = seq(0.1, 1.5, by=0.1))
ElbowPlot(kpkpl, ndims = 50) ;
DimPlot(kpkpl, group.by="sample", cols=longchao_cols) ;
DimPlot(kpkpl, group.by="RNA_snn_res.0.3", label=TRUE, cols=llc_cols) ;
## set cluster ID
cluster.ids.kpkpl <- sort(as.numeric(unique(as.character(kpkpl@meta.data$RNA_snn_res.0.3))))
annotation <- c('Monocytes/MDSCs', 'Macrophages', 'Cd4 T Cells', 'Regulatory T Cells','Malignant Cells', 'Cd8 T Cells','Double Negative T Cells', 
                'Natural Killer Cells', 'Dendritic Cells','Inhibitory T Cells')
kpkpl@meta.data$annotation <- kpkpl@meta.data$RNA_snn_res.0.3
kpkpl@meta.data$annotation <- plyr::mapvalues(x = kpkpl@meta.data$annotation, 
                                                         from = cluster.ids.kpkpl, to = annotation)
DimPlot(kpkpl, group.by="annotation", label = TRUE, cols= llc_cols, pt.size = 0.5) ;
DimPlot(treatment.update, group.by="sample", label = FALSE, cols= longchao_cols, pt.size = 0.5);
## heatmap
tmp_copy <- kpkpl
My_levels <- c('Monocytes/MDSCs', 'Macrophages', 'Cd4 T Cells', 'Regulatory T Cells',
               'Malignant Cells', 'Cd8 T Cells','Double Negative T Cells', 
               'Natural Killer Cells', 'Dendritic Cells','Inhibitory T Cells')
Idents(tmp_copy) <- factor(Idents(tmp_copy), levels = My_levels)
tmp_copy$annotation <- factor(x =tmp_copy$annotation, levels = c('Monocytes/MDSCs', 'Macrophages', 'Cd4 T Cells', 'Regulatory T Cells',
                                                                  'Cd8 T Cells','Double Negative T Cells', 
                                                                 'Natural Killer Cells', 'Inhibitory T Cells','Dendritic Cells','Malignant Cells'))
markers.to.plot <- c("S100a9", "S100a8", "Cd14","Adgre1", "Naaa","C1qc", "Cd3e", "Cd3g", "Cd3d","Cd4","Foxp3","Ctla4","Cd8a", "Cd8b1", "Gzmb", "Klrb1c", "Ncr1","Flt3", "Batf3", "Il33", "Krt8", "Krt18") #"S100a4", "Itgae", "Ccr7" 
DoHeatmap(tmp_copy, features = markers.to.plot, size = 3, group.by="annotation", group.colors = c('#b3d3e8','#7dcade','#368dab', '#259b39', '#a2ce89', '#e4361e', '#93646d', '#c0a26e', '#edb000','#65b72f','#eb6938'), disp.min = -5) + scale_fill_gradientn(colors = c("blue", "white", "red")) + guides(color=FALSE)

## Write marker tables
kpkpl@meta.data$adj_cluster <- as.character(kpkpl@meta.data$RNA_snn_res.0.2) ;
kpkpl.markers <- FindAllMarkers(kpkpl, group.by= "adj_cluster", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
setwd("C:/Users/lixx4/OneDrive/Desktop")
write.csv(kpkpl.markers, file="kpkpl.markers.csv", row.names=TRUE) ;

## select out T cells
t.cell_m <- kpkpl@meta.data ;
apply(t.cell_m[, 5:ncol(t.cell_m)], 2, function(x) range(as.numeric(x))) ;
DimPlot(kpkpl, group.by="RNA_snn_res.0.5", label=TRUE) ;
FeaturePlot(kpkpl, features=c("Cd3e", "Cd3g", "Cd3d"), col=mycolors3)

kpkpl_Tcells <- subset(kpkpl, subset=RNA_snn_res.0.5%in%c(3, 8, 2, 7, 11)) ;
FeaturePlot(kpkpl_Tcells, features=c("Cd3e", "Cd3g", "Cd3d"), col=mycolors3)

#
kpkpl_Tcells <- NormalizeData(kpkpl_Tcells, normalization.method = "LogNormalize", scale.factor = 10000)
kpkpl_Tcells <- FindVariableFeatures(kpkpl_Tcells, selection.method = "vst", nfeatures = 2000)
kpkpl_Tcells <- RunPCA(kpkpl_Tcells, features = VariableFeatures(object=kpkpl_Tcells)) ;

kpkpl_Tcells <- RunUMAP(kpkpl_Tcells, reduction = "pca", dims = 1:50)
kpkpl_Tcells <- FindNeighbors(kpkpl_Tcells, reduction = "pca", dims = 1:50)
kpkpl_Tcells <- FindClusters(kpkpl_Tcells, resolution = 0.5)

DimPlot(kpkpl_Tcells, group.by="RNA_snn_res.0.5", label=TRUE) ;
DimPlot(kpkpl_Tcells, group.by="sample", label=TRUE) ;
FeaturePlot(kpkpl_Tcells, feature=c("Cd8a", "Cd4", "Foxp3", "Tcf7", "Fcer1g"), col=mycolors3)
dimnames(kpkpl_Tcells)
###############
kpkpl.markers <- FindAllMarkers(kpkpl_Tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

kpkpl.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top5 <- kpkpl.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

DoHeatmap(kpkpl, features = top5$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))

setwd("C:/Users/lixx4/OneDrive/Desktop")

write.table(kpkpl.markers, file="kpkpl.tmarkers.txt", sep = ",")
########

t.exp_m <- kpkpl_Tcells@assays$RNA@data ; 
Cd8_exp <- apply(t.exp_m[c("Cd8a", "Cd8b1"), ], 2, mean) ;

kpkpl_Tcells[["Cd8_status"]] <- Cd8_exp > 0;
kpkpl_CD8Tcells <- subset(kpkpl_Tcells, subset=Cd8_status) ;
kpkpl_CD8Tcells <- subset(kpkpl_CD8Tcells, subset = Cd4 < 0.01)
#
VlnPlot(kpkpl_CD8Tcells, feature="Tcf7", group.by="sample", pt.size=0, cols=(values=c("#403ffc","#ff7e7c")))
display.brewer.all()
VlnPlot(kpkpl_CD8Tcells, feature="Tcf7", group.by="sample", col=brewer.pal(3, "RdBu") , pt.size=0) ;

kpkpl_CD8Tcells <- RunUMAP(kpkpl_CD8Tcells, reduction = "pca", dims = 1:10, seed.use = 72) ;
kpkpl_CD8Tcells <- FindNeighbors(kpkpl_CD8Tcells, reduction = "pca", dims = 1:50)
kpkpl_CD8Tcells <- FindClusters(kpkpl_CD8Tcells, resolution = 0.2)
FeaturePlot(kpkpl_CD8Tcells, feature= "Tcf7", col=mycolors3)
FeaturePlot(kpkpl_CD8Tcells, feature= "Foxp3", col=mycolors3)
my_cols <- c('KP9_3'='#403ffc','KPL9_3_1'='#ff7e7c')
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
DimPlot(kpkpl_CD8Tcells, cols=my_cols2, group.by="sample", label=FALSE) ;
DimPlot(kpkpl_CD8Tcells, group.by="RNA_snn_res.0.2", label=TRUE, cols=col.cluster) ;

## set cluster ID
cluster.ids <- sort(as.numeric(unique(as.character(kpkpl_CD8Tcells@meta.data$RNA_snn_res.0.2))))
annotation <- c('Terminal Exhausted T', 'Central Memory T', 'Proliferating T', 'Inhibitory T')
kpkpl_CD8Tcells@meta.data$annotation <- kpkpl_CD8Tcells@meta.data$RNA_snn_res.0.2
kpkpl_CD8Tcells@meta.data$annotation <- plyr::mapvalues(x = kpkpl_CD8Tcells@meta.data$annotation, 
                                            from = cluster.ids, to = annotation)
DimPlot(kpkpl_CD8Tcells, group.by="annotation", label = TRUE, cols= longchao_cols, pt.size = 1) 

## Heatmap based on markers of your own choice
tmp_copy <- kpkpl_CD8Tcells
My_levels <- c('Terminal Exhausted T', 'Central Memory T', 'Proliferating T', 'Inhibitory T')
Idents(tmp_copy) <- factor(Idents(tmp_copy), levels = My_levels)
tmp_copy$RNA_snn_res.0.2 <- factor(x =tmp_copy$RNA_snn_res.0.2, levels = c('Terminal Exhausted T', 'Central Memory T', 'Proliferating T', 'Inhibitory T'))
markers.to.plot <- c("Lag3", "Tgfb1", "Pdcd1","Tcf7", "Klf2","Sell", "Mki67", "Cdca8", "Rrm1", "Foxp3","Klrg1", "Ikzf2") 
DoHeatmap(tmp_copy, features = markers.to.plot, size = 3, group.by="annotation", group.colors = c('#E7A427', '#3D9F43', '#3E8BA9', '#E96946', '#3E8BA9', '#B59671', '#E4CD83'), disp.min = -1) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + guides(color=FALSE)

## Visualize each group individually
t.cell_m1 <-kpkpl_CD8Tcells@meta.data ;
Hcells1 <- rownames(t.cell_m1)[t.cell_m1$sample=="KP9_3"] ;
DimPlot(kpkpl_CD8Tcells, label=FALSE, pt.size=1, cells.highlight=Hcells1, cols.highlight = "#0020AA") ;
Hcells2 <- rownames(t.cell_m1)[t.cell_m1$sample=="KPL9_3_1"] ;
DimPlot(kpkpl_CD8Tcells, label=FALSE, pt.size=1, cells.highlight=Hcells2, cols.highlight = "#0020AA") ;
###
features <- c("Tcf7","Pdcd1","Foxp3")
VlnPlot(kpkpl_CD8Tcells, features=features)

###############
kpkpl.markers <- FindAllMarkers(kpkpl_CD8Tcells, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

kpkpl.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top5 <- kpkpl.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

DoHeatmap(kpkpl, features = top5$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))

setwd("C:/Users/lixx4/OneDrive/Desktop")

write.table(kpkpl.markers, file="kpkpl.cd8tmarkers.txt", sep = ",")
########
kpkpl_CD8Tcells
pbmc.markers <- FindAllMarkers(kpkpl_CD8Tcells, only.pos = TRUE, logfc.threshold = 0)
cluster1.marker <- pbmc.markers[pbmc.markers$cluster==1, ] ;

t.index <- grep("^Rps|^Rpl", rownames(cluster1.marker))
cluster1.marker.update <- cluster1.marker[-t.index, ] ;


##
#setwd("/endosome/archive/bioinformatics/Li_lab/s186258/Huiyu_Project/RData") ;
#save.image(file="Tempo.RData") ;

###
kpkpl_CD8Tcells
pbmc.markers2 <- FindAllMarkers(kpkpl_CD8Tcells, logfc.threshold = 0, only.pos = FALSE)
cluster1.markerb <- pbmc.markers2[pbmc.markers2$cluster==1, ] ;

t.index <- grep("^Rps|^Rpl", rownames(cluster1.markerb))
cluster1.marker.updateb <- cluster1.markerb[-t.index, ] ;
diff.genes <- FindMarkers(cluster1.marker.updateb,min.pct = 0, logfc.threshold = 0)

cluster1.marker.updateb %>% top_n(n = 2, wt = avg_logFC)

write.table(cluster1.marker.updateb, file="kpkpl.cd8.violin4.txt", sep = ",")
table(kpkpl_CD8Tcells@meta.data$sample, kpkpl_CD8Tcells@meta.data$seurat_clusters)

## Violin plot
voldata <- read.delim(file = "kpkpl.cd8.violin4.txt", header = TRUE, row.names = NULL)

voldata$threshold <- as.factor(ifelse(voldata$p_val <0.01 & abs(voldata$avg_logFC) >=0.25,ifelse(voldata$avg_logFC > 0.25 ,'Up','Down'),'Not'))

g <- ggplot(data=voldata, aes(x=avg_logFC, y=logPvalue,  colour=threshold)) +
  scale_color_manual(values=c("blue", "dimgrey","red")) +
  geom_point(alpha=0.5, size=1.75) +
  xlab("log2 fold change") + ylab("-log10 p-value")+
  theme(legend.position="none")

g

