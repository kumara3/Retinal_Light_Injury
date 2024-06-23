rm(list=ls())
options(echo=TRUE)

# usage:
usage = "Usage: Rscript script.R  <A text file of geneList with sample names to be drawn in heatmap> <outdir>\n"
args <- commandArgs(trailingOnly = TRUE)

if(length(args)!=2) {
  
  stop("\n    Wrong parameters. \n    ", usage)
}

filename <- args[1]
outdir <- args[2]
dir.create(outdir,showWarnings = TRUE)

heatmapWithhc <- function(A)
{
  library(dynamicTreeCut) 
  library(gplots)
  library(MASS)
  library(RColorBrewer)
  library(reshape2)
  
  #Import the data as Matrix
  mM <- data.matrix(na.omit(read.table(A, row.names = 1, header = TRUE, sep = "\t", na.strings = "NA")))
  
  
  MyMatrix <- mM
  
  #Perform hierarchical clustering, one can update the distance matrix to something else if required.
  hr <- hclust(as.dist(1-cor(t(MyMatrix), method="spearman")), method="complete")
  
  hcd <- as.dendrogram(hr)
  
	png(file = paste(A, ".hc_clust.dendo.png", sep = ""),width=680, height=680, units="px")
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                  #cex = 0.7, col = "blue")
  # Customized plot; remove labels
  plot(hcd, xlab = "Height", nodePar = nodePar, leaflab = "none",
       horiz = TRUE, edgePar = list(col = 2:3, lwd = 2:1))
	plot(hr)
  plot(rect.hclust(hr,h=1))
	dev.off()
  
  #cuttre cuts the dendrogram or the hclust object as per K or h
 
	mycl <- cutree(hr, h=0.1)
	
 ## Trying the dynamic cut tree from WGCNA ###
	mycl <- cutreeDynamicTree(hr,deepSplit = TRUE, minModuleSize = 30) ## This uses the tree method
  cor_d = as.matrix(as.dist(1-cor(t(MyMatrix))))
 	
  mycl <- cutreeDynamic(hr,deepSplit = 4, minClusterSize = 20,method="hybrid",distM=cor_d) ## gives samller cluster (importat parameters to be varied are deepsplit and method)
	#rect.hclust(hr,h=1)
  
  # or simply add the cluster ID to your data
  MyMatrixC <- cbind(MyMatrix, clusterID=mycl)
  
  MyMatrixC <- apply(MyMatrixC[hr$order,], 2,rev)
  
  col = colnames(MyMatrixC)
  row = rownames(MyMatrixC)
  tab = cbind(row,data.frame(MyMatrixC))
  colnames(tab) = c("GENE_NAME",col)
  write.table(tab, file=paste(outdir,"/",A,".hc_clust.txt", sep = ""), sep="\t",row.names=F)
  
  #These are for various color options
  clusterCols <- rainbow(length(unique(mycl)))
	#myClusterSideBar <- clusterCols[mycl]   ## does not work as the cluster size can be 0.
  
	myClusterSideBar <- colorRampPalette(clusterCols)(length(mycl))[rank(mycl)] ## works with cluster size as 0. Check if cluster number 0 has any other meaning.
  bk = unique(c(seq(-5.73,-0.5, length=100),seq(-0.5,0, length=100), seq(0,5,length=100)))
  hmcols<- colorRampPalette(c("darkblue","lightblue", "gold"))(length(bk)-1)
  
  
  #Output to a pdf
  pdf(file = paste(outdir,"/",A, ".hc_clust.pdf", sep = ""), width = 10, height = 10)
  
  ## Actual heatmap command
#heatmap plot : https://stackoverflow.com/questions/52599180/partial-row-labels-heatmap-r
  hm = heatmap.2(MyMatrix, 
                 main = "Heatmap (Average expression)",
                 Rowv=as.dendrogram(hr),
                 Colv=NA, 
                 dendrogram="row",
                 scale="row",
                 RowSideColors = myClusterSideBar,
                 col=hmcols,
                 density.info="none",
                 trace="none",
#labRow = c("Apoe","mt-Atp6","Glul","mt-Co3","mt-Co2","Arc","Areg","Flnc","Foxf1","Gpr3"),
                 #labRow = "",
                 lhei = c(1,5),
                 srtCol = 45,
								 cexRow = 2.5,
                 cexCol =2.5,
                 margins = c(10,10))
  
  centered_data = t(scale(t(MyMatrixC[,-ncol(MyMatrixC)]))) # center rows, mean substracted
  hc_genes = as.dendrogram(hr)
	wil.data = MyMatrixC[,-ncol(MyMatrixC)]  

  gene_partition_assignments <- MyMatrixC[,"clusterID"]
  #gnames_in_cluster_order = row.names(centered_data)[order.dendrogram(hc_genes) ]
  partition_in_cluster_order = unique(MyMatrixC[,"clusterID"])
  
  gene_names = rownames(centered_data)
  num_cols = ncol(centered_data)
  
	#pdf(file="my_cluster_plots.tc.pdf", width=8, height=14)
  par(mfrow=c(3, 2))
	par(cex.lab=1.5,cex.main=1.5,cex.axis=1.5,font.lab=2)
  par(mar=c(7,5,4,2))
  
  for (i in seq_along(partition_in_cluster_order)) {
		clust_name_id = partition_in_cluster_order[i]
    partition_i = (gene_partition_assignments == partition_in_cluster_order[i])
    partition_centered_data = centered_data[partition_i,] 
   # if the looping cluster id has only one gene and belong to that cluster ID, then it returns a vector instead of a table
    if (sum(partition_i) == 1) {
      dim(partition_centered_data) = c(1,num_cols)
      colnames(partition_centered_data) = colnames(centered_data)
      rownames(partition_centered_data) = gene_names[partition_i]
    }
    partition_centered_data_long <- melt(partition_centered_data, id.vars=c(""))
    colnames(partition_centered_data_long) <- c("Gene","condition","zscore")
			
    ###one way anova test . This is a parametric test for comparing the means of multiple groups ##
    ### see http://www.sthda.com/english/wiki/one-way-anova-test-in-r ###
    
    res.aov <- aov(as.numeric(partition_centered_data_long$`zscore`)~as.factor(partition_centered_data_long$condition),data=partition_centered_data_long)
    Tukey <- TukeyHSD(res.aov)

    res.wil <- pairwise.wilcox.test(as.numeric(partition_centered_data_long$`zscore`), as.factor(partition_centered_data_long$condition) ,p.adjust.method = "BH")
		PT <- res.wil$p.value
    #PT1 <- as.data.frame(fullPTable(PT))
    PT1 <- as.data.frame(PT)
		

	### wilcoxon t test. This is a non parametric test and can be used to test the unpaired indepenedent sample###
  ## This test is used by Seurat as default ###
    
    all <- rbind(Tukey$`as.factor(partition_centered_data_long$condition)`,c("","",dim(partition_centered_data)[[1]],paste("1-way Anova p-value=",summary(res.aov)[[1]][["Pr(>F)"]][[1]],sep="")))
    outfile = paste(outdir,"/",A, ".subcluster_", clust_name_id, ".txt", sep='')
    outfile2 = paste(outdir,"/",A, ".subcluster_", clust_name_id, ".Tukeyvalues.txt", sep='')
    outfile3 = paste(outdir,"/",A, ".subcluster_", clust_name_id, ".pairwise.wilcoxon.txt", sep='')
    write.table(partition_centered_data, file=outfile, quote=F, sep="\t", col.names=NA)
    write.table(all, file=outfile2, quote=F, sep="\t", col.names=NA)
    write.table(PT1, file=outfile3, quote = F, sep="\t", col.names = NA)
		
    ## plotting the gene clusters in the form of trend box plot ###
    data = partition_centered_data
    plot_label = paste(length(data[,1]), " regions in subcluster ", partition_in_cluster_order[i], sep='')
    #colrs=c('beige','bisque','bisque3','bisque4')
		colrs = c("cadetblue1")
		par(bty='n')
    #boxplot(t(data), type="l", lty=1, col=rgb(.3, .3, .3, .3), ylim=c(min(data)-0.5,max(data)+0.5), notch=TRUE, main=plot_label, xaxt='n', xlab='',ylab='Normalized expression')
    #matplot(t(data), type="l", lty=1, ylim=c(min(data)-0.5,max(data)+0.5), notch=TRUE, main=plot_label, xaxt='n', xlab='',ylab='Standardized expression', col=colrs)
		matplot(t(data), type="l", lty=1, ylim=c(-5,5), main=plot_label, xaxt='n', xlab='',ylab='Standardized expression', col=colrs)
    axis(side=1, at=1:ncol(data), labels=colnames(data), las=2) #http://rfunction.com/archives/1302
		abline(v=0,lty=2,col="gray")
    points(as.numeric(colMeans(data)), type='p', col="red", lwd=2, pch=19)
    #points(as.numeric(colMeans(data)), type='p', col=rep(brewer.pal(5, "Dark2"), each=3), lwd=2, pch=19)
    lines(as.numeric(colMeans(data)), type='c', col="black", lwd=3)
  }

  dev.off()
}

heatmapWithhc(filename)
