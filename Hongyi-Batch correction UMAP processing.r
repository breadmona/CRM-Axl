library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

col.cluster2.1 <- c('#1C79B7','#F38329','#2DA248','#DC403E','#976BA6','#8D574C','#D07DB0','#BDBF3B',
                    '#27BDD0','#B0C7E5','#F9BA79','#A0CC8A','#F09EA1','#C7B0D1','#D4BCB8','#F7CCDD',
                    '#DBDD8D','#A3D7E5','#B04E4F','#A38A58','#ED5351','#0C8945','#0F71A7','#A82764',
                    '#F8DAE4','#7A4B1F','#5E6BAE','#8AC997','#DAC9AC','#0F5148','#A0B5DC','#9F858E',
                    '#5C181D','#7B8380','#E8DDDA','#264220','#5AB747','#5169AE','#4B3D58','#CD428A',
                    '#62615A','#B82129','#66762E')
###Load 1st batch data
treatment.data <- Read10X(data.dir = "c:\\Users\\lixx4\\OneDrive\\Desktop\\KPLtreat1\\") ;
tmp <- colnames(treatment.data) ; 
tmp1 <- sapply(strsplit(tmp, "-"), function(x) x[1]) ;
tmp2 <- sapply(strsplit(tmp, "-"), function(x) x[2]) ;
library(plyr) ;
tmp3 <- revalue(tmp2, c("1"="Ctrl", "2"="PD_1", "3"="Combined")) ;

tmp4 <- paste(tmp1, tmp3, sep="-") ;
colnames(treatment.data) <- tmp4 ;

treatment <- CreateSeuratObject(counts = treatment.data, project = "treatment", min.cells = 3, min.features = 100) ;
treatment@meta.data$sample <- sapply(strsplit(colnames(treatment), "-"), function(x) x[2]);

treatment[["percent.mt"]] <- PercentageFeatureSet(treatment, pattern = "^mt-")
VlnPlot(treatment, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
treatment <- subset(treatment, subset = nFeature_RNA > 200 & percent.mt < 10) ; 

### Load 2nd batch data
bgb.data <- Read10X(data.dir = "c:\\Users\\lixx4\\OneDrive\\Desktop\\bgb\\") ;
tmp <- colnames(bgb.data) ; 
tmp1 <- sapply(strsplit(tmp, "-"), function(x) x[1]) ;
tmp2 <- sapply(strsplit(tmp, "-"), function(x) x[2]) ;
library(plyr) ;
tmp3 <- revalue(tmp2, c("1"="BGB324")) ;
tmp4 <- paste(tmp1, tmp3, sep="-") ;
colnames(bgb.data) <- tmp4 ;

bgb <- CreateSeuratObject(counts = bgb.data, project = "bgb", min.cells = 3, min.features = 100) ;
bgb@meta.data$sample <- sapply(strsplit(colnames(bgb), "-"), function(x) x[2]);

#Human reference genome is ^MT-
bgb[["percent.mt"]] <- PercentageFeatureSet(bgb, pattern = "^mt-")
VlnPlot(bgb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
bgb <- subset(bgb, subset = nFeature_RNA > 200 & percent.mt < 10) ; 

###Normalize seurat object
treatment <- NormalizeData(treatment, verbose = FALSE)
treatment <- FindVariableFeatures(treatment, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
bgb <- NormalizeData(bgb, verbose = FALSE)
bgb <- FindVariableFeatures(bgb, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

###Find anchors for integration
reference.list <- c("treatment","bgb")
kpl.anchors <- FindIntegrationAnchors(object.list = list(treatment, bgb), dims = 1:50)

kpl.integrated <- IntegrateData(anchorset=kpl.anchors, dims = 1:50)

#Scaling the data
all.genes <- rownames(kpl.integrated)
treatment <- ScaleData(kpl.integrated, features = all.genes)

#Linear dimensional reduction
treatment <- RunPCA(treatment, features = VariableFeatures(object=treatment))

# t-SNE and Clustering
treatment <- RunUMAP(treatment, reduction = "pca", dims = 1:50)
DimPlot(treatment, group.by="sample") ;

treatment <- FindNeighbors(treatment, reduction = "pca", dims = 1:50)
treatment <- FindClusters(treatment, resolution = seq(0.1, 1.5, by=0.1))
ElbowPlot(treatment, ndims = 50) ;

t.cell_n <- treatment@meta.data ;
apply(t.cell_n[, 5:ncol(t.cell_n)], 2, function(x) range(as.numeric(x))) ;
DimPlot(treatment, group.by="integrated_snn_res.0.2", label=TRUE, cols=llc_cols) ;

## Visualize general UMAP for supplement
treatment.update <- subset(treatment, subset=integrated_snn_res.0.2%in%c(0, 1,2, 3, 4, 5, 6, 7, 8, 9)) ;
DimPlot(treatment.update, group.by="integrated_snn_res.0.2", label=TRUE, cols=llc_cols) ;

## set cluster ID
cluster.ids.treatment <- sort(as.numeric(unique(as.character(treatment.update@meta.data$integrated_snn_res.0.2))))
annotation <- c('Monocytes/MDSCs', 'Macrophages', 'Regulatory T Cells', 'Cd4 T Cells', 'Cd8 T Cells', 
                'Malignant Cells', 'Inhibitory T Cells', 'Natural Killer Cells','Double Negative T Cells', 
                'Dendritic Cells')
treatment.update@meta.data$annotation <- treatment.update@meta.data$integrated_snn_res.0.2
treatment.update@meta.data$annotation <- plyr::mapvalues(x = treatment.update@meta.data$annotation, 
                                            from = cluster.ids.treatment, to = annotation)
DimPlot(treatment.update, group.by="annotation", label = TRUE, cols= llc_cols, pt.size = 0.5) ;
DimPlot(treatment.update, group.by="sample", label = FALSE, cols= longchao_cols, pt.size = 0.5);

treatment.update.markers <- FindAllMarkers(treatment.update, group.by= "integrated_snn_res.0.2", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- treatment.update.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(treatment.update, features = top5$gene) + NoLegend()

tmp_copy <- treatment.update
My_levels <- c('Monocytes/MDSCs', 'Macrophages', 'Regulatory T Cells', 'Cd4 T Cells', 'Cd8 T Cells', 
               'Malignant Cells', 'Inhibitory T Cells', 'Natural Killer Cells','Double Negative T Cells', 
               'Dendritic Cells')
Idents(tmp_copy) <- factor(Idents(tmp_copy), levels = My_levels)
tmp_copy$annotation <- factor(x =tmp_copy$annotation, levels = c('Monocytes/MDSCs', 'Macrophages', 'Regulatory T Cells', 'Cd4 T Cells', 'Cd8 T Cells', 
                                                                   'Inhibitory T Cells', 'Natural Killer Cells','Double Negative T Cells', 
                                                                   'Dendritic Cells', 'Malignant Cells'))
markers.to.plot <- c("S100a9", "S100a8", "Cd14","Adgre1", "Naaa","C1qc", "Foxp3","Ctla4","Cd3e", "Cd3g", "Cd3d","Cd4", "Cd8a", "Cd8b1", "Gzmb", "Klrb1c", "Ncr1","Flt3", "Batf3","Il33", "Krt8", "Krt18") #"S100a4", "Itgae", "Ccr7" 
DoHeatmap(tmp_copy, features = markers.to.plot, size = 3, group.by="annotation", group.colors = c('#b3d3e8','#7dcade','#368dab', '#259b39', '#65b72f', '#e4361e', '#93646d', '#c0a26e', '#edb000','#a2ce89', '#eb6938'), disp.min = -2.5) + scale_fill_gradientn(colors = c("blue", "white", "red")) + guides(color=FALSE)

## Select out T cells for downstream analysis
FeaturePlot(treatment, features=c("Cd3e", "Cd3g", "Cd3d"), col=mycolors3)
FeaturePlot(treatment, features=c("Klrb1c"), col=mycolors3)


treatment_Tcells <- subset(treatment, subset=integrated_snn_res.0.6%in%c(12, 4, 5, 3, 7, 9)) ;
FeaturePlot(treatment_Tcells, features=c("Cd3e", "Cd3g", "Cd3d"), col=mycolors3)

t.exp_n <- treatment_Tcells@assays$RNA@data ; 
Cd8_exp <- apply(t.exp_n[c("Cd8a", "Cd8b1"), ], 2, mean) ;

treatment_Tcells[["Cd8_status"]] <- Cd8_exp > 0;

treatment_CD8Tcells <- subset(treatment_Tcells, subset=Cd8_status) ;
treatment_CD8Tcells <- NormalizeData(treatment_CD8Tcells, normalization.method = "LogNormalize", scale.factor = 10000)
treatment_CD8Tcells <- FindVariableFeatures(treatment_CD8Tcells, selection.method = "vst", nfeatures = 2000)
treatment_CD8Tcells <- RunPCA(treatment_CD8Tcells, features = VariableFeatures(object=treatment_Tcells)) ;

treatment_CD8Tcells <- RunUMAP(treatment_CD8Tcells, reduction = "pca", dims = 1:10, seed.use=1) ;
DimPlot(treatment_CD8Tcells, group.by="sample") ;
FeaturePlot(treatment_CD8Tcells, features=c("Foxp3", "Cd4", "Cd8a", "Cd8b1"), col=mycolors3) ;

treatment_CD8Tcells <- FindNeighbors(treatment_CD8Tcells, reduction = "pca", dims = 1:10)
treatment_CD8Tcells <- FindClusters(treatment_CD8Tcells, resolution = seq(0.1, 1.5, by=0.1))

DimPlot(treatment_CD8Tcells, group.by="integrated_snn_res.0.4", label=TRUE) ; #0.2, 0.3, 0.4, 0.5

# filter out non-CD8 T cells
FeaturePlot(treatment_CD8Tcells, feature=c("Cd8a", "Cd4", "Foxp3", "Tcf7", "Fcer1g"), col=mycolors3)
treatment_filter_CD8Tcells <- subset(treatment_CD8Tcells, subset = Cd4 < 0.01)
#treatment_filter_CD8Tcells <- subset(treatment_CD8Tcells, subset=integrated_snn_res.0.4%in%c(0,1,2,5,6,7)) ;
treatment_filter_CD8Tcells <- RunUMAP(treatment_filter_CD8Tcells, reduction = "pca", dims = 1:10, seed.use=62) ;
DimPlot(treatment_filter_CD8Tcells, group.by="integrated_snn_res.0.4", label=TRUE) ; #0.2, 0.3, 0.4, 0.5
treatment_filter_CD8Tcells@meta.data$adj_cluster <- as.character(treatment_filter_CD8Tcells@meta.data$integrated_snn_res.0.4) ;

##plot use tran R
t.XY <- treatment_filter_CD8Tcells@reductions$umap@cell.embeddings ;
t.cell_m <- treatment_filter_CD8Tcells@meta.data ; 
t.cluster <- t.cell_m$adj_cluster ;
used.colors <- mycolors1[1:6] ;
names(used.colors) <- unique(t.cluster) ;
t.color <- used.colors[t.cluster] ;
plot(t.XY[,1:2], col=t.color, asp=1, cex=0.5, pch=21, xlab="UMAP(Dm1)", ylab="UMAP(Dm2)") ;
legend("topright", names(used.colors), col=used.colors, pch=21) ;

## Write marker tables
DimPlot(treatment_filter_CD8Tcells, group.by="sample", label=TRUE) ; #0.2, 0.3, 0.4, 0.5
FeaturePlot(treatment_filter_CD8Tcells, feature=c("Tcf7"), col=mycolors3)
FeaturePlot(treatment_filter_CD8Tcells, feature=c("S1pr1"), col=mycolors3)

VlnPlot(treatment_filter_CD8Tcells, feature="Tcf7", group.by="sample")
treatment_filter_CD8Tcells.markers <- FindAllMarkers(treatment_filter_CD8Tcells, group.by= "adj_cluster", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
setwd("C:/Users/lixx4/OneDrive/Desktop")
write.csv(treatment_filter_CD8Tcells.markers, file="treatment_filter_CD8Tcells.markers.csv", row.names=TRUE) ;

## ROE analysis
ROE.theory <- table(treatment_filter_CD8Tcells@meta.data$adj_cluster)
write.csv(ROE.theory, file="ROE.theory.csv", row.names=TRUE) ;
ROE.actual<- table(treatment_filter_CD8Tcells@meta.data$sample, treatment_filter_CD8Tcells@meta.data$adj_cluster)
write.csv(ROE.actual, file="ROE.actual.csv", row.names=TRUE) 
ROE <- chisq.test(ROE.actual)
names(ROE)
stats <- ROE$observed/ROE$expected
stats

## Adjust UMAP for better visualization
setwd("/project/bioinformatics/Li_lab/s186258/boli_lab/Temporary/Huiyu_Li/Data")
load("KPLtreat.Rdata") ;

DimPlot(treatment_filter_CD8Tcells, group.by="adj_cluster", label=TRUE) ; #0.2, 0.3, 0.4, 0.5
tmp <- RunUMAP(treatment_filter_CD8Tcells, reduction = "pca", dims = 1:15, seed.use = 70, group.by="adj_cluster") ;
DimPlot(tmp, group.by="adj_cluster", label = TRUE, cols= col.cluster, pt.size = 1) ;
DimPlot(tmp, group.by="sample", label=TRUE) ; #0.2, 0.3, 0.4, 0.5

## set cluster ID
cluster.ids <- sort(as.numeric(unique(as.character(tmp@meta.data$adj_cluster))))
annotation <- c('Bystander T', 'Helper T', 'Exhaustive Activated T', 'Proliferating T', 
                'Inhibitory T','Stem T', 'Proliferating Transitory T', 
                'Clonal Expanded T')
tmp@meta.data$annotation <- tmp@meta.data$adj_cluster
tmp@meta.data$annotation <- plyr::mapvalues(x = tmp@meta.data$annotation, 
                                            from = cluster.ids, to = annotation)
DimPlot(tmp, group.by="annotation", label = TRUE, cols= col.cluster, pt.size = 1) ;
FeaturePlot(tmp, feature=c("Tcf7", "S1pr1"), col=mycolors3)
FeaturePlot(tmp, feature=c("S1pr1"), col=mycolors3)
DotPlot(tmp, features = c("Tcf7","Pdcd1", "S1pr1"), group.by="annotation") + RotatedAxis()
DotPlot(tmp, features = c("Tcf7","Pdcd1", "S1pr1"), group.by="sample") + RotatedAxis()

## Visualize cell location in each treatment group
t.cell_m <- treatment_filter_CD8Tcells@meta.data ;
Hcells <- rownames(t.cell_m)[t.cell_m$sample=="Combined"] ;
DimPlot(tmp, label=FALSE, pt.size=1, cells.highlight=Hcells, cols.highlight = "#0020AA") ;


## Heatmap based on markers of your own choice
tmp_copy <- tmp
My_levels <- c('Bystander T', 'Helper T', 'Exhaustive Activated T', 'Proliferating T', 
               'Inhibitory T','Stem T', 'Proliferating Transitory T', 
               'Clonal Expanded T')
Idents(tmp_copy) <- factor(Idents(tmp_copy), levels = My_levels)
tmp_copy$adj_cluster <- factor(x =tmp_copy$adj_cluster, levels = c('Bystander T', 'Helper T', 'Exhaustive Activated T', 'Proliferating T', 
                                                                   'Inhibitory T','Stem T', 'Proliferating Transitory T', 
                                                                   'Clonal Expanded T'))
markers.to.plot <- c("Ccl5", "S100a6", "Ccl4","Ccr7", "Tnfrsf4","Tnfsf4", "Pdcd1", "Lag3", "Gzma", "Gzmb","Trbv14", "Trbv13-2", "Ccna2", "Mki67", "Cdca8","Tcf7", "Klf2", "Sell", "Mcm3", "Mcm5", "Mcm7","Trbv29", "Trav12-2") #"S100a4", "Itgae", "Ccr7" 
DoHeatmap(tmp_copy, features = markers.to.plot, size = 3, group.by="adj_cluster", group.colors = c('#E5D2DD', '#F3B1A0', '#F1BB72', '#57C3F3', '#D6E7A3','#E95C59','#476D87', '#53A85F'), disp.min = -1) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + guides(color=FALSE)

##Zhaoning code for heatmap
Idents(tmp_copy) <- "annotation"
table(tmp_copy@active.ident)
My_levels <- c('Bystander T', 'Helper T', 'Exhaustive Activated T', 'Proliferating T', 
               'Inhibitory T','Stem T', 'Proliferating Transitory T', 
               'Clonal Expanded T')
markers.to.plot <- c("Ccl5", "S100a6", "Ccl4","Ccr7", "Tnfrsf4","Tnfsf4", "Pdcd1", "Lag3", "Gzma", "Gzmb","Trbv14", "Trbv13-2", "Ccna2", "Mki67", "Cdca8","Tcf7", "Klf2", "Sell", "Mcm3", "Mcm5", "Mcm7","Trbv29", "Trav12-2") #"S100a4", "Itgae", "Ccr7" 
DoHeatmap(tmp_copy, features = markers.to.plot, size = 3, group.colors = c('#E5D2DD', '#F3B1A0', '#F1BB72', '#57C3F3', '#D6E7A3','#E95C59','#476D87', '#53A85F'), disp.min = -1) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + guides(color=FALSE)
## Add Module score                                                                                                  
marker <- list(c('Tcf7','S1pr1'))
tmp <- AddModuleScore(object = tmp, features = marker, ctrl = 5, name = 'X2gene_score')
FeaturePlot(object = tmp, features = "X2gene_score1", cols=mycolors3)

## Add TCR-TRB info to seurat object metadata
setwd("/endosome/archive/bioinformatics/Li_lab/s186258/Huiyu_Project/RData") ;
load("TCR_info.RData");

t.cell_m <- treatment_filter_CD8Tcells@meta.data ;
t.cell_m$TCRB <- TCR_info[match(rownames(t.cell_m), TCR_info$cell_id), "cdr3"] ;
t.cell_m$cloneSize <- TCR_info[match(rownames(t.cell_m), TCR_info$cell_id), "cloneSize"] ;

treatment_filter_CD8Tcells@meta.data$TCRB <- t.cell_m$TCRB ;
treatment_filter_CD8Tcells@meta.data$cloneSize <- t.cell_m$cloneSize ;

## Clone expansion TCR analysis
Clonectrl <- subset(treatment_filter_CD8Tcells, sample=="Ctrl") ; 
clonesize <- table(Clonectrl@meta.data$TCRB,Clonectrl@meta.data$adj_cluster, Clonectrl@meta.data$cloneSize)
Clonebgb <- subset(treatment_filter_CD8Tcells, sample=="BGB324") ; 
clonesizebgb <- table(Clonebgb@meta.data$TCRB,Clonebgb@meta.data$adj_cluster, Clonebgb@meta.data$cloneSize)
Clonepd1 <- subset(treatment_filter_CD8Tcells, sample=="PD_1") ; 
clonesizepd1 <- table(Clonepd1@meta.data$TCRB,Clonepd1@meta.data$adj_cluster, Clonepd1@meta.data$cloneSize)
Clonecom <- subset(treatment_filter_CD8Tcells, sample=="Combined") ; 
clonesizecom <- table(Clonecom@meta.data$TCRB,Clonecom@meta.data$adj_cluster, Clonecom@meta.data$cloneSize)

table(treatment_filter_CD8Tcells@meta.data$sample, treatment_filter_CD8Tcells$cloneSize)

setwd("C:/Users/lixx4/OneDrive/Desktop")
write.csv(clonesizecom, file="clone distributioncom.csv", row.names=TRUE) 

## TCR shared clone type analysis - result should be symmetrical
###cluster7 and cluster1
cluster7_TCRs <- t.cell_m$TCRB[ t.cell_m$adj_cluster==7] ;
cluster7_TCRs <- cluster7_TCRs[!is.na(cluster7_TCRs)] ;

cluster1_TCRs <- t.cell_m$TCRB[ t.cell_m$adj_cluster==1] ;
cluster1_TCRs <- cluster1_TCRs[!is.na(cluster1_TCRs)] ;

C7_clonetypes <- unique(cluster7_TCRs) ;
C1_clonetypes <- unique(cluster1_TCRs) ;

numerator <- length(intersect(C7_clonetypes, C1_clonetypes)) ;
denominator <- length(C7_clonetypes) * length(C1_clonetypes) ;
C=10000
(numerator/denominator) * C #14.8

save.image(file = "c:\\Users\\lixx4\\OneDrive\\Desktop\\plustTCR.RData")


###
IFNG <- rnorm(100, 1, 1) ; # top50:C1, tail50:C2
#Ifng
IFNG <- subset(treatment, features = "Ifng") ; 
IFNG <- table(treatment@meta.data$nCount_RNA,treatment@meta.data$adj_cluster)
#@metadata$adj_cluster
IFNGR1 <- rnorm(100, 1, 1) ;
#Ifngr1

Cluster_index <- c(rep("C1", 50), rep("C2", 50)) ; 

Cluster_index <- sort(as.numeric(unique(as.character(tmp@meta.data$adj_cluster))))
IFNG_avg <- tapply(IFNG,Cluster_index, mean) ;
IFNGR1_avg <- tapply(IFNGR1, Cluster_index, mean) ;

E_Lc1 <- IFNG_avg["C1"] ;
E_Rc1 <- IFNGR1_avg["C1"]

Score <- E_Lc1 * E_Rc1 ;

exp_m <- rbind(IFNG, IFNGR1) ;

Rscore_list <- c() ;
for(i in 1:10000)
{
  t.Rindex <-sample(1:100, 100) ;
  t.Rexp <- exp_m[, t.Rindex] ;
  
  RIFNG_avg <- tapply(t.Rexp[1,], Cluster_index, mean) ;
  RIFNGR1_avg <- tapply(t.Rexp[2,], Cluster_index, mean) ;
  
  RE_Lc1 <- RIFNG_avg["C1"] ;
  RE_Rc1 <- RIFNGR1_avg["C1"]
  
  RScore <- RE_Lc1 * RE_Rc1 ;
  Rscore_list <- c(Rscore_list, RScore)
  print(i)
}

sum(Rscore_list > Score) /100

#exp 19000 * 20000
a <- exp[c("INFG", "TCF7"), ]

C1_index <- which(meta.data$cluster=="C1")
C1_cells <- rownames(meta.data)[C1_index] 




