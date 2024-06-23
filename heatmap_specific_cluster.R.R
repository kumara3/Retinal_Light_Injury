library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
library(GetoptLong)
library(tools)
library(stringi)


setwd("D:\\data")
file_names <- list.files(".", pattern="\\.subcluster_[0-9]*\\.txt$")
add_name <- file_path_sans_ext(file_names[1])   # Define prefix for the heatmap
clus_name <- gsub("\\.heatmap.*","",add_name)

plots_increase = list()
plots_decrease = list()
for(i in 1:length(file_names)){
  df <- read.table(file_names[i], sep="\t", header = TRUE, row.names = 1)
  df_mat <- as.matrix(df)
  
  ## First trial, edit when required. ##
  #df4_mat <- df3_mat
  #df4_mat <- df4_mat[order(df4_mat[,4],decreasing = TRUE),]
  #df5_mat <- df4_mat[order(df4_mat[,3],decreasing = TRUE),]
  # select the genes which shows changes at Day5 #
  #increase <- head(rownames(df4_mat),10)
  #decrease <- tail(rownames(df4_mat),10)
  ##########################

  # select the genes which shows changes between Day5 and baseline
  select_genes_day5 = df_mat[,4] - df_mat[,1]
  cc_day5_increase <- names(select_genes_day5[order(select_genes_day5, decreasing = TRUE)])[1:5]  ### if the change is increase in day5

  cc_day5_decrease <- names(select_genes_day5[order(select_genes_day5, decreasing = FALSE)])[1:5]  ### if the change is decrease in day5
  
    # select the genes which shows changes between 48hrs and baseline
  select_genes_48hrs = df_mat[,3] - df_mat[,1]
  cc_48hrs_increase <- names(select_genes_48hrs[order(select_genes_48hrs, decreasing = TRUE)])[1:5]  ### if the change is increase in 48hrs

  cc_48hrs_decrease <- names(select_genes_48hrs[order(select_genes_48hrs, decreasing = FALSE)])[1:5]  ### if the change is decrease in 48hrs
  
  ##prepare for heatmap Day5
  ccl_day5_increase = rownames(df_mat) %in% cc_day5_increase 
  genes_repr_day5_increase = rownames(df_mat)[ccl_day5_increase]  

  ccl_day5_decrease = rownames(df_mat) %in% cc_day5_decrease 
  genes_repr_day5_decrease = rownames(df_mat)[ccl_day5_decrease]  
  
  
  ## prepare for heatmap 48hrs
  ccl_48hrs_increase = rownames(df_mat) %in% cc_48hrs_increase
  genes_repr_48hrs_increase = rownames(df_mat)[ccl_48hrs_increase]

  ccl_48hrs_decrease = rownames(df_mat) %in% cc_48hrs_decrease
  genes_repr_48hrs_decrease = rownames(df_mat)[ccl_48hrs_decrease]
  
  
  #gene_list <- c(cc_gene_increase,cc_gene_decrease)

  prefix <- file_path_sans_ext(file_names[i])   # Define prefix for the heatmap
  clus <- stri_sub(prefix, -1)
  
  plots_increase[[i]] = Heatmap(df_mat, col = colorRamp2(c(-2,-1,0,1,2), c("blue","lightblue3","floralwhite", "gold","gold4")), 
             name = "scaled_expr", column_title = qq("Expression for @{nrow(df_mat)} genes in cluster",clus),
             show_column_names = TRUE, width = unit(8, "cm"),show_row_dend = FALSE,
             heatmap_legend_param = list(at=c(-2,-1,0,1,2),color_bar="continuous",title = "Scaled Expression"),cluster_columns = FALSE,cluster_rows = TRUE)+
  Heatmap(ccl_day5_increase + 0, name = "Day5", col = c("0" = "white", "1" = "brown"), 
          show_heatmap_legend = FALSE, width = unit(5, "mm"))+
  Heatmap(ccl_48hrs_increase + 0, name = "48hrs", col = c("0" = "white", "1" = "red"), 
          show_heatmap_legend = FALSE, width = unit(5, "mm"))+
  rowAnnotation(mark = anno_mark(at = which(ccl_day5_increase|ccl_48hrs_increase ==TRUE),labels = rownames(df_mat)[ccl_day5_increase|ccl_48hrs_increase ==TRUE],
                                 labels_gp = gpar(fontsize = 10,fontface="bold"), padding = unit(0.25, "mm"), which = "row"))
  

  
  
  plots_decrease[[i]] = Heatmap(df_mat, col = colorRamp2(c(-2,-1,0,1,2), c("blue","lightblue3","floralwhite", "gold","gold4")), 
                                name = "scaled_expr", column_title = qq("Expression for @{nrow(df_mat)} genes in cluster",clus),
                                show_column_names = TRUE, width = unit(8, "cm"),show_row_dend = FALSE,
                                heatmap_legend_param = list(at=c(-2,-1,0,1,2),color_bar="continuous",title = "Scaled Expression"),cluster_columns = FALSE,cluster_rows = TRUE)+
    Heatmap(ccl_day5_decrease + 0, name = "Day5", col = c("0" = "white", "1" = "brown"), 
            show_heatmap_legend = FALSE, width = unit(5, "mm"))+
    Heatmap(ccl_48hrs_decrease + 0, name = "48hrs", col = c("0" = "white", "1" = "red"), 
            show_heatmap_legend = FALSE, width = unit(5, "mm"))+
    rowAnnotation(mark = anno_mark(at = which(ccl_day5_decrease|ccl_48hrs_decrease ==TRUE),labels = rownames(df_mat)[ccl_day5_decrease|ccl_48hrs_decrease ==TRUE],
                                   labels_gp = gpar(fontsize = 10,fontface="bold"), padding = unit(0.25, "mm"), which = "row"))
  
  
  
  }

temp <- paste(clus_name,".increase_48hrs_Day5",".pdf",sep="")
pdf(temp, height=8, width=8)
for(i in 1:length(plots_increase)){
  draw(plots_increase[[i]])
}
dev.off()


temp <- paste(clus_name,".decrease_48hrs_Day5",".pdf",sep="")
pdf(temp, height=8, width=8)
for(i in 1:length(plots_decrease)){
  draw(plots_decrease[[i]])
}
dev.off()
