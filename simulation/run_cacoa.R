library(cacoa)
library(coda.base)
library(Seurat)

run_cacoa <- function(data_dir, results_dir, 
                      results_file_name, 
                      condition_label,
                      sample_label,
                      celltype_label,
                      ref_label,
                      target_label) {
  sce_sim <- readRDS(data_dir) 
  
  sample_per_cell <- colData(sce_sim)[, sample_label]
  names(sample_per_cell) <- colnames(sce_sim)
  condition_group <- colData(sce_sim)[, condition_label]
  names(condition_group) <- colnames(sce_sim)
  names(sce_sim$Condition) <- colnames(sce_sim)
  celltype_per_cell <- colData(sce_sim)[, celltype_label]
  names(celltype_per_cell) <- colnames(sce_sim)
  #names(sce_sim$scClassify_tumour_prediction_coarse) <- colnames(sce_sim)
  
  meta_data <- colData(sce_sim)[, c(condition_label, sample_label)]
  meta_data <- unique(meta_data)
  meta_data <- data.frame(meta_data)
  colnames(meta_data) <- c("Condition", "Sample")
  sample.groups <- meta_data$Condition
  names(sample.groups) <- meta_data$Sample
  counts(sce_sim) <- logcounts(sce_sim)
  
  sc_dataset <- Seurat::as.Seurat(sce_sim)
  names(sc_dataset@assays) <- "RNA"
  DefaultAssay(sc_dataset) <- "RNA"
  sc_dataset <- Seurat::FindVariableFeatures(object = sc_dataset, nfeatures = 3000, 
                                             selection.method = "vst", verbose = TRUE)
  sc_dataset <- Seurat::ScaleData(object = sc_dataset, verbose = TRUE)
  sc_dataset <- Seurat::RunPCA(object = sc_dataset, features = VariableFeatures(sc_dataset), 
                               verbose = TRUE)
  sc_dataset <- Seurat::FindNeighbors(object = sc_dataset, dims = 1:40, 
                                      k.param = 100, # instead of default k = 20
                                      prune.SNN = 0, # instead of default prune.SNN = 1/15
                                      verbose = TRUE)
  sc_dataset <- Seurat::RunUMAP(object = sc_dataset, dims = 1:10, 
                                verbose = TRUE)
  
  cao <- cacoa::Cacoa$new(
    sc_dataset, 
    sample.groups = sample.groups, 
    sample.per.cell = sample_per_cell,
    cell.groups = celltype_per_cell,
    target.level = target_label, 
    ref.level = ref_label, 
    n.cores = 1, verbose = FALSE,
    embedding = reducedDim(sce_sim, "UMAP"),
    graph.name = "RNA_snn"
  )
  
  
  cao$estimateCellDensity(method = 'graph') # Use graph so that there is cell level score
  cao$estimateDiffCellDensity() # This is running the permutation test
  cao$estimateDiffCellDensity(type = "wilcox")
  cacoa_results <- cbind(pertmutation_raw = cao$test.results$cell.density$diff$permutation$raw,
                         pertmutation_adj = cao$test.results$cell.density$diff$permutation$adj,
                         wilcox_row = cao$test.results$cell.density$diff$wilcox$raw,
                         wilcox_adj = cao$test.results$cell.density$diff$wilcox$adj)
  rownames(cacoa_results) <- colnames(sc_dataset)
  saveRDS(cao$test.results, file = file.path(results_dir, gsub(".csv", ".rds", results_file_name)))
  write.csv(cacoa_results,
            file = file.path(results_dir, results_file_name), row.names = FALSE)
  return(cacoa_results)
}
