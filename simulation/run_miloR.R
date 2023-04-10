library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
# 
# data_dir = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/data/sce_sim_test.rds"
# results_dir = "results"
# results_file_name = "miloR_output.csv"
# condition_label = "Condition"
# sample_label = "Sample"



run_miloR <- function(data_dir, results_dir, 
                      results_file_name, condition_label,
                      sample_label) {
  
  #Create a Milo object
  sce_sim <- readRDS(data_dir) 
  milo_object <- Milo(sce_sim)
  
  
  milo_object <- buildGraph(milo_object, k = 10, d = 30)
  
  
  milo_object <- makeNhoods(milo_object, prop = 0.1, k = 10, d=30, refined = TRUE)
  plotNhoodSizeHist(milo_object)
  
  
  milo_object <- countCells(milo_object, meta.data = data.frame(colData(milo_object)), sample = sample_label)
  
  ## Defining experimental design
  
  milo_design <- data.frame(colData(milo_object))[, c(sample_label, condition_label)]
  colnames(milo_design) <- c("Sample", "Condition")
  ## Convert batch info from integer to factor
  
  milo_design <- distinct(milo_design)
  rownames(milo_design) <- milo_design$Sample
  
  ## Computing neighbourhood connectivity
  
  milo_object <- calcNhoodDistance(milo_object, d=30, reduced.dim = "PCA")
  
  
  ## Testing
  
  da_results <- testNhoods(milo_object, design = ~ Condition, design.df = milo_design,
                           reduced.dim = "PCA")

  
  nhoods_mat <- milo_object@nhoods
  nhoods_mat <- apply(nhoods_mat, 2, function(x) x/sum(x))
  milo_scores <- da_results$logFC %*% t(nhoods_mat)
  milo_scores <- milo_scores[1, ]
  
  saveRDS(list(da_results = da_results, nhoods_mat = nhoods_mat), file = file.path(results_dir, gsub(".csv", ".rds", results_file_name)))
  write.csv(milo_scores,
            file = file.path(results_dir, results_file_name), row.names = FALSE)
  return(milo_scores)
}
