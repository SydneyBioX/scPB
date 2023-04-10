

library(reticulate)
print(py_config())


library(DEGAS)
library(Rtsne)
library(ggplot2)
library(data.table)
library(matrixStats)
source("utils_sim.R")




# data_dir = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/data/sce_sim_test.rds"
# results_dir = "results"
# results_file_name = "DEGAS_output.csv"
# condition_label = "Condition"
# sample_label = "Sample"

# zscore normalization
normFunc <- function(x){return((x-mean(x, na.rm = T))/(sd(x, na.rm = T)+1e-3))}

# scaling from 0-1
scaleFunc <- function(x){return((x- min(x, na.rm = T)) /(max(x, na.rm = T)-min(x, na.rm = T)+1e-3))}

# Preprocess count data
normalizeScale <-function(X){
  return(t(apply(t(apply(as.matrix(t(X)),1,normFunc)),1,scaleFunc)))
}


preprocessLogCounts <- function(X){
  return(normalizeScale(1.5^X))  
}

DEGAS <- function(sce_sim, bulk_dataset, meta_bulk, condition_label, celltype_label) {
  
  patient_matrix <- as.data.frame(model.matrix(~meta_bulk[, condition_label] - 1))
  # Rename columns to match levels
  colnames(patient_matrix) <- gsub("meta_bulk\\[, condition_label\\]", "", colnames(patient_matrix))
  
  # Apply the function to each row of the data frame
  #sc_data <-  t(apply(logcounts(sce_sim), 1, range01))
  sc_data <- preprocessLogCounts(logcounts(sce_sim))
  #sc_matrix <- as.data.frame(model.matrix(~colData(sce_sim)[, celltype_label] - 1))
  #colnames(sc_matrix) <- gsub("meta_bulk\\[, celltype_label\\]", "", colnames(sc_matrix))
  
  print(dim(sc_data))
  scLab <- toOneHot(colData(sce_sim)[, celltype_label])
  #bulk_dataset <-  t(apply(bulk_dataset, 1, range01))
  bulk_dataset <- preprocessLogCounts(bulk_dataset)

  ccModel1 = runCCMTLBag(scExp = sc_data, scLab = scLab, 
                         patExp = bulk_dataset, patLab = patient_matrix,
                         tmpDir,
                         model_type = 'ClassClass',
                         architecture = 'DenseNet',
                         FFdepth = 3,
                         Bagdepth= 5)
  scpatPreds <- predClassBag(ccModel1, t(sc_data), 'pat')
  colnames(scpatPreds) <- colnames(patient_matrix)
  return(scpatPreds)
}


run_DEGAS <- function(data_dir, results_dir, 
                      results_file_name, condition_label, sample_label, 
                      celltype_label, 
                      celltype_prop, 
                      num_per_sample,
                      ref_label, target_label) {
  
  sce_sim <- readRDS(data_dir) 

  initDEGAS()
  tmpDir <-  file.path(results_dir, "/DEGAS/temp/")
  dir.create(tmpDir, recursive = TRUE)
  set_seed_term(1)
  
  
  
  bulk_res <- constructBulkFromSC(sce_sim, sample_label, condition_label, 
                                  celltype_label = celltype_label,
                                  num_per_sample = num_per_sample, 
                                  celltype_prop = celltype_prop) 
  meta_bulk <- bulk_res$meta_bulk
  bulk_dataset <- bulk_res$bulk_dataset
  
  meta_bulk[, condition_label] <- factor(meta_bulk[, condition_label], levels = c(ref_label, target_label))
  # 
  # patient_matrix <- as.data.frame(model.matrix(~meta_bulk[, condition_label] - 1))
  # 
  # # Rename columns to match levels
  # colnames(patient_matrix) <- gsub("meta_bulk\\[, condition_label\\]", "", colnames(patient_matrix))
  # 
  # 
  # sc_matrix <- as.data.frame(model.matrix(~colData(sce_sim)[, celltype_label] - 1))
  # colnames(sc_matrix) <- gsub("meta_bulk\\[, celltype_label\\]", "", colnames(sc_matrix))
  # 
  # print(dim(sc_matrix))
  # bulk_dataset <-  t(apply(bulk_dataset, 1, range01))
  # ccModel1 = runCCMTLBag(scExp = t(sc_data), scLab = sc_matrix, 
  #                        patExp = t(bulk_dataset), patLab = patient_matrix,
  #                        tmpDir,
  #                        model_type = 'ClassClass',
  #                        architecture = 'DenseNet',
  #                        FFdepth = 3,
  #                        Bagdepth= 5)
  # scpatPreds <- predClassBag(ccModel1, t(sc_data), 'pat')
  # colnames(scpatPreds) <- colnames(patient_matrix)
  # 
  
  scpatPreds <- DEGAS(sce_sim, bulk_dataset, meta_bulk, condition_label, celltype_label)
  write.csv(scpatPreds, file = file.path(results_dir, results_file_name), row.names = FALSE)
  
  return(scpatPreds)
  
  
}

