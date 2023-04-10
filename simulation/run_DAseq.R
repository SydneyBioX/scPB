
library(DAseq)
library(BiocNeighbors)
library(BiocParallel)
# data_dir = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/data/sce_sim_test.rds"
# results_dir = "results"
# results_file_name = "DAseq_output.csv"
# condition_label = "Condition"


run_DAseq <- function(data_dir, results_dir, 
                      results_file_name, 
                      condition_label,
                      k.vector = seq(50, 500, 50)) {
  
  sce_sim <- readRDS(data_dir) 
  X.label.info <- colData(sce_sim)
  condition_group <- names(table(colData(sce_sim)[, condition_label]))
  condition <- X.label.info[, condition_label]
  labels.1 <- colnames(sce_sim)[condition == condition_group[1]]
  labels.2 <- colnames(sce_sim)[condition == condition_group[2]]
  
  da_cells <- getDAcells(
    X = reducedDim(sce_sim, "PCA"), # input merged dataset of interest after dimension reduction.
    cell.labels = colnames(sce_sim),
    labels.1 = labels.1,
    labels.2 = labels.2,
    k.vector = k.vector,
    plot.embedding = reducedDim(sce_sim, "UMAP"),
    do.plot = FALSE,
    save.knn = FALSE
  )
  
  

  write.csv(da_cells$da.pred,
            file = file.path(results_dir, results_file_name), row.names = FALSE)
  return(da_cells)
  
}
