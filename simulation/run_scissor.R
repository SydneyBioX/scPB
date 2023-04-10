library(Scissor)
library(SingleCellExperiment)
source("utils_sim.R")

# data_dir = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/data/sce_sim_test.rds"
# results_dir = "results"
# results_file_name = "scissor_output.csv"
# condition_label = "Condition"
# sample_label = "Sample"




run_scissor <- function(data_dir, results_dir, 
                        results_file_name, condition_label, sample_label, 
                        celltype_label, 
                        celltype_prop, 
                        num_per_sample) {
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
  phenotype <- as.numeric(as.factor(phenotype)) - 1
  
  infos1 <- Scissor(as(bulk_dataset, "matrix"), sc_dataset, phenotype, 
                    alpha = 0.05, 
                    tag = names(table(meta_bulk[, condition_label])),
                    family = "binomial", 
                    Save_file = gsub(".csv", ".RData", file.path(results_dir, results_file_name)))
  
  write.csv(infos1$Coefs, file = file.path(results_dir, results_file_name), row.names = FALSE)
  return(infos1)
}

