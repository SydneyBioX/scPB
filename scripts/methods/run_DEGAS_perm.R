
# Set the Python virtual environment
library(GenomeInfoDb)
# Load reticulate library
library(reticulate)

# Verify the Python configuration
py_config()


library(DEGAS)
library(Rtsne)
library(ggplot2)
library(data.table)
library(matrixStats)
library(SingleCellExperiment)
source("scripts/utils_sim.R")


initDEGAS()
set_seed_term(1)


# zscore normalization
normFunc <- function(x) {
  return((x - mean(x, na.rm = T)) / (sd(x, na.rm = T) + 1e-3))
}

# scaling from 0-1
scaleFunc <- function(x) {
  return((x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T) + 1e-3))
}

# Preprocess count data
normalizeScale <- function(X) {
  return(t(apply(t(apply(as.matrix(t(X)), 1, normFunc)), 1, scaleFunc)))
}


preprocessLogCounts <- function(X) {
  return(normalizeScale(1.5^X))
}

DEGAS <- function(sce_sim, bulk_dataset, meta_bulk, condition_label, celltype_label, results_dir, seed = NULL) {
  tmpDir <- file.path(results_dir, paste0("/DEGAS/temp_", seed, "/"))
  dir.create(tmpDir, recursive = TRUE)
  
  
  
  patient_matrix <- as.data.frame(model.matrix(~ meta_bulk[, condition_label] - 1))
  # Rename columns to match levels
  colnames(patient_matrix) <- gsub("meta_bulk\\[, condition_label\\]", "", colnames(patient_matrix))
  

  sc_data <- preprocessLogCounts(logcounts(sce_sim))
  
  scLab <- toOneHot(colData(sce_sim)[, celltype_label])
  bulk_dataset <- preprocessLogCounts(bulk_dataset)
  
  
  ccModel1 <- runCCMTLBag(
    scExp = sc_data, scLab = scLab,
    patExp = bulk_dataset, patLab = patient_matrix,
    tmpDir,
    model_type = "ClassClass",
    architecture = "DenseNet",
    FFdepth = 3,
    Bagdepth = 5
  )
  scpatPreds <- predClassBag(ccModel1, sc_data, "pat")
  colnames(scpatPreds) <- colnames(patient_matrix)
  unlink(tmpDir, recursive = TRUE)
  return(scpatPreds)
}


run_DEGAS <- function(data_dir, results_dir,
                      results_file_name, condition_label, sample_label,
                      celltype_label,
                      celltype_prop,
                      num_per_sample,
                      ref_label, target_label,
                      perm = FALSE,
                      nperm = 1000,
                      BPPARAM = SerialParam()) {
  sce_sim <- readRDS(data_dir)
  
  
  
  
  
  bulk_res <- constructBulkFromSC(sce_sim, sample_label, condition_label,
                                  celltype_label = celltype_label,
                                  num_per_sample = num_per_sample,
                                  celltype_prop = celltype_prop
  )
  meta_bulk <- bulk_res$meta_bulk
  bulk_dataset <- bulk_res$bulk_dataset
  
  meta_bulk[, condition_label] <- factor(meta_bulk[, condition_label], levels = c(ref_label, target_label))
  

  
  scpatPreds <- DEGAS(sce_sim, bulk_dataset, meta_bulk, condition_label, celltype_label, results_dir)
  
  if (perm) {
    perm_res <- BiocParallel::bplapply(1:nperm, function(x) {
      DEGAS(sce_sim, bulk_dataset, meta_bulk, sample(condition_label), celltype_label, results_dir, seed = x)
    }, BPPARAM = BPPARAM)
    write.csv(do.call(rbind, perm_res),
              file = file.path(results_dir, gsub(".csv", "_perm.csv", results_file_name)),
              row.names = FALSE
    )
  } 
  
  write.csv(scpatPreds, file = file.path(results_dir, results_file_name), row.names = FALSE)

  
  return(scpatPreds)
}


