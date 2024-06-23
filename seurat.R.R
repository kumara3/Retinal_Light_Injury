#install.packages('Seurat')

##Standard pre-processing workflow
# The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat.
# These represent the creation of a Seurat object,
# the selection and filtration of cells based on QC metrics,
# data normalization and scaling, and
# the detection of highly variable genes.

# usage:

usage1 = "

Usage:
Rscript Basic.Seurat.R (followed by options below)
    1. </FullPath/TenXCountDir> Can be either of below, please try and give FullPath
                                outs/raw_gene_bc_matrices/mm10/
                                outs/filtered_gene_bc_matrices/mm10/
    2. </FullPath/OutDir>       E.g. SeuratAnalysis/
    3. <OutPrefix>              E.g. SampleName
    4. <Genes of interest text file>     If you want to define the clusters as per your genes of interest 
					 E.g. genes_of_interest.txt
"

args = commandArgs(trailingOnly = TRUE)
if(length(args)!=4) {
  cat(usage1)
  stop("\n Wrong parameters. See usage above.\n")
}

options(echo=TRUE)
system("ulimit -n 10000")
suppressPackageStartupMessages(library(Seurat))
library(dplyr)

TenX_countDir = args[1]
output_dir = args[2]
dir.create(file.path(output_dir), showWarnings = TRUE)
prefix = args[3]
genes_of_interest = args[4]

# Load the dataset
TenX.data <- Read10X(data.dir = TenX_countDir)

##Reading the genes of interset file ###
interested_genes <- character(0)
open_genes_of_interest <- file(genes_of_interest,'r')
linn <- readLines(open_genes_of_interest)
for(i in 1:length(linn)){
  interested_genes <- c(interested_genes,linn[i])
}

#interested_genes
#length_interested_genes=length(interested_genes)


# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
TenX <- CreateSeuratObject(counts = TenX.data, min.cells = 3, min.features = 200, project = paste0("10X_",prefix))
#colnames(x = TenX.data) <- paste('TenX.data', colnames(x = TenX.data), sep = '_')


mito.features <- grep(pattern = "^MT-|^mt-", x = rownames(x = TenX), value = TRUE)
percent.mt <- Matrix::colSums(x = GetAssayData(object = TenX, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = TenX, slot = 'counts'))

# The [[ operator can add columns to object metadata, and is a great place to stash QC stats
TenX[['percent.mt']] <- percent.mt

tmp <- paste(output_dir,"/", prefix, ".nRNA_mito_Seurat.ViolinPlot.pdf",sep="")
pdf(tmp, height=10 , width=10)
VlnPlot(object = TenX, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
dev.off()

#FeatureScatter is typically used to visualize feature-feature relationships

tmp <- paste(output_dir,"/", prefix, ".Seurat.FeatureScatterPlot.pdf",sep="")
pdf(tmp, height=10 , width=10)
par(mfrow = c(1, 2))
FeatureScatter(object = TenX, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = TenX, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
print(paste0("The number of cells to start with is  :",length(TenX@meta.data$orig.ident)))
## Save the RDS file ###
tmp <- paste(output_dir,"/", prefix, ".Seurat.raw.TenX.rds",sep="")
saveRDS(TenX, file = tmp)


#Normalizing the data
TenX <- NormalizeData(object = TenX, normalization.method = "LogNormalize", scale.factor = 10000)


##Detection of variable genes across the single cells
TenX <- FindVariableFeatures(object = TenX, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(TenX), 10)
tmp <- paste(output_dir,"/", prefix, ".Seurat.VariableFeaturesDispersionPlot.pdf",sep="")
pdf(tmp, height=10 , width=10)
plot1 <- VariableFeaturePlot(TenX)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
print(paste0("The number of variable features selected for the downstreaming analysis",length(VariableFeatures(TenX))))
dev.off()

#Scaling the data and removing unwanted sources of variation
TenX <- ScaleData(object = TenX, features = rownames(TenX), vars.to.regress = "percent.mt")

## Save the normalized and scaled data ###
tmp <- paste(output_dir,"/", prefix, ".SeuratNormalized.Scaled.TenX.rds",sep="")
saveRDS(TenX, file = tmp)
#TenX <- readRDS(tmp)

### To go further we need to look for the distribution of different markesr in cluster of cells ###
TenX <- RunPCA(object = TenX, features = VariableFeatures(object = TenX))
tmp <- paste(output_dir,"/", prefix, ".Seurat.First5PCs.10geneEach.txt",sep="")
sink(tmp)
print(x = TenX[['pca']], dims= 1:5, nfeatures = 10, projected = FALSE)
sink()

tmp <- paste(output_dir,"/", prefix, ".Seurat.dim_1to5.PCAplot.pdf",sep="")
pdf(tmp, height=15, width=15)
VizDimLoadings(object = TenX, dims = 1:5,reduction = "pca")
dev.off()

#TenX <- ProjectDim(object = TenX)


tmp <- paste(output_dir,"/", prefix, ".Seurat.PCHeatmaps.pdf",sep="")
pdf(tmp, height=15, width=15)
DimHeatmap(object = TenX, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = TenX, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


#####Determine statistically significant principal components####

TenX <- JackStraw(object = TenX, num.replicate = 100)

## Save the .rds file for jackstraw ####
tmpA <- paste(output_dir,"/", prefix, ".Seurat.jackstraw.TenX.rds",sep="")
saveRDS(TenX, file = tmpA)
TenX <- ScoreJackStraw(object = TenX, dims=1:20)

### Plot the jackstraw ####
tmp <- paste(output_dir,"/", prefix, ".Seurat.JackStrawPlot.pdf",sep="")
pdf(tmp, height=15, width=15)
JackStrawPlot(object = TenX, dims = 1:15)
ElbowPlot(object = TenX)
dev.off()


######Cluster cells##################################
TenX <- FindNeighbors(object = TenX, dims = 1:10)
TenX <- FindClusters(object = TenX, resolution = 0.5)
TenX <- RunTSNE(object = TenX, dims = 1:10)
tmp <- paste(output_dir,"/", prefix, ".Seurat.tSNE.pdf",sep="")
pdf(tmp, height=15, width=15)
DimPlot(object = TenX, reduction = 'tsne',pt.size=2, label=TRUE, label.size=6)
dev.off()

################### Finding the differentialy expressed features ###################
TenX.markers <- FindAllMarkers(object = TenX,min.pct = 0.25)
tmpB <- paste(output_dir,"/", prefix, ".Seurat.FindAllMarkers.rds",sep="")
saveRDS(TenX, file = tmpB)

TenX.markers.print<-TenX.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
tmp <- paste(output_dir,"/", prefix, ".Markers_DE.txt",sep="")
write.table(TenX.markers, file=tmp, sep="\t",row.names=FALSE)

##########################Violin Plot for the markers ##############################
markers.to.plot=unique(TenX.markers.print$gene)
tmp <- paste(output_dir,"/", prefix, ".Seurat.VlnPlot.pdf",sep="")
pdf(tmp, height=15, width=15)
VlnPlot(object = TenX, features = markers.to.plot)
dev.off()

tmp <- paste(output_dir,"/", prefix, ".Seurat.FeaturePlots.pdf",sep="")
pdf(tmp, height=15, width=15)
FeaturePlot(object = TenX , features = markers.to.plot)
dev.off()

tmp <- paste(output_dir,"/", prefix, ".Seurat.DoHeatMapPlots.pdf",sep="")
pdf(tmp, height=15, width=15)
DoHeatmap(object = TenX , features = markers.to.plot)
dev.off()


