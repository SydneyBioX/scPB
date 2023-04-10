library(Seurat)
library(preprocessCore)
library(scAB)


# data_dir = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/data/sce_sim_test.rds"
# results_dir = "results"
# results_file_name = "scAB_output.csv"
# condition_label = "Condition"
# sample_label = "Sample"

run_scAB <- function(data_dir, results_dir, 
                     results_file_name, condition_label, sample_label, 
                     celltype_label, 
                     celltype_prop, 
                     num_per_sample,
                     target_label,
                     ref_label) {
  
  
  sce_sim <- readRDS(data_dir) 
  counts(sce_sim) <- logcounts(sce_sim)
  sc_dataset <- as.Seurat(sce_sim)
  names(sc_dataset@assays) <- "RNA"
  DefaultAssay(sc_dataset) <- "RNA"
  sc_dataset <- Seurat::FindVariableFeatures(object = sc_dataset, nfeatures = 3000, 
                                             selection.method = "vst", verbose = TRUE)
  sc_dataset <- Seurat::ScaleData(object = sc_dataset, verbose = TRUE)
  sc_dataset <- Seurat::RunPCA(object = sc_dataset, features = VariableFeatures(sc_dataset), 
                               verbose = TRUE)
  sc_dataset <- Seurat::FindNeighbors(object = sc_dataset, dims = 1:40, 
                                      verbose = TRUE)
  sc_dataset <- Seurat::RunUMAP(object = sc_dataset, dims = 1:10, 
                                verbose = TRUE)
  
  
  
  # Create bulk data
  
  
  bulk_res <- constructBulkFromSC(sce_sim, sample_label, condition_label, 
                                  celltype_label = celltype_label, 
                                  celltype_prop = celltype_prop, 
                                  num_per_sample = num_per_sample)
  meta_bulk <- bulk_res$meta_bulk
  bulk_dataset <- bulk_res$bulk_dataset
  
  phenotype <- meta_bulk[, condition_label]
  names(phenotype) <- rownames(meta_bulk)
  
  print(table(phenotype, ifelse(phenotype == target_label, 0, 1)))
  phenotype <- ifelse(phenotype == target_label, 0, 1)
  
  
  scAB_data <- create_scAB(sc_dataset, bulk_dataset, phenotype, method = "binary")
  
  
  K <- select_K(scAB_data)
  K
  
  scAB_result <- scAB(Object=scAB_data, K=K)
  
  sc_dataset <- findSubset(sc_dataset, scAB_Object = scAB_result, tred = 2)
  
  
  # table(sc_dataset@meta.data$scAB_select, sc_dataset$Condition)
  
  # UMAP_scAB <- DimPlot(sc_dataset,group.by="scAB_select",cols=c("#80b1d3","red"),pt.size=0.001,order=c("scAB+ cells","Other cells"))
  # patchwork::wrap_plots(plots = list(UMAP_celltype,UMAP_scAB), ncol = 2)
  # 
  # g1 <- ggplot(sc_dataset@meta.data, aes(x = Subset1_loading, y = Subset2_loading, color = scAB_select)) +
  #   geom_point() +
  #   theme_bw()
  # 
  # g2 <- ggplot(sc_dataset@meta.data, aes(x = Subset1_loading, y = Subset2_loading, color = scAB_Subset1)) +
  #   geom_point() +
  #   theme_bw()
  # 
  # 
  # g3 <- ggplot(sc_dataset@meta.data, aes(x = Subset1_loading, y = Subset2_loading, color = scAB_Subset2)) +
  #   geom_point() +
  #   theme_bw()
  # ggpubr::ggarrange(g1, g2, g3, ncol = 2, nrow = 2, align = "hv")
  # 
  # scAB_loading_scores <- apply(sc_dataset@meta.data[, grep("loading", colnames(sc_dataset@meta.data))], 1, max)
  # ggplot(sc_dataset@meta.data, aes(x = scAB_loading_scores, y = Subset2_loading, color = scAB_select)) +
  #   geom_point() +
  #   theme_bw()
  # 
  # ggplot(sc_dataset@meta.data, aes(x = scAB_loading_scores, y = Subset1_loading, color = scAB_select)) +
  #   geom_point() +
  #   theme_bw()
  
  scAB_loading_scores <- apply(sc_dataset@meta.data[, grep("loading", colnames(sc_dataset@meta.data))], 1, max)
  write.csv(cbind(scAB_scores = scAB_loading_scores, 
                  scAB_select = sc_dataset@meta.data$scAB_select, 
                  sc_dataset@meta.data[, grep("loading", colnames(sc_dataset@meta.data))]), 
            file = file.path(results_dir, results_file_name), row.names = FALSE)
  return(scAB_loading_scores)
}
