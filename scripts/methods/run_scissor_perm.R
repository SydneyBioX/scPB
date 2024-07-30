library(Scissor)
library(SingleCellExperiment)
source("scripts/utils_sim.R")


run_scissor <- function(data_dir, results_dir, 
                        results_file_name, condition_label, sample_label, 
                        celltype_label, 
                        celltype_prop, 
                        num_per_sample,
                        perm = FALSE,
                        nperm = 100,
                        BPPARAM = BiocParallel::SerialParam(),
                        subset = 1,
                        top_hvg = 2000,
                        profile_mem = FALSE) {
  sce_sim <- readRDS(data_dir) 
  counts(sce_sim) <- logcounts(sce_sim)

  if (subset < 1) {
    sce_sim <- sce_sim[, sample(ncol(sce_sim), round(ncol(sce_sim) * subset))]
  }

  sc_dataset <- as.Seurat(sce_sim)
  names(sc_dataset@assays) <- "RNA"
  DefaultAssay(sc_dataset) <- "RNA"
  sc_dataset <- Seurat::FindVariableFeatures(object = sc_dataset, nfeatures = top_hvg, 
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
                                  num_per_sample = num_per_sample,
                                  prop = 0.6)
  meta_bulk <- bulk_res$meta_bulk
  bulk_dataset <- bulk_res$bulk_dataset
  
  phenotype <- meta_bulk[, condition_label]
  names(phenotype) <- rownames(meta_bulk)
  phenotype <- as.numeric(as.factor(phenotype)) - 1
  
  if (profile_mem) {
    Rprof("scissor_mem.out", memory.profiling=T, interval=0.02)
    start_time <- Sys.time()
  }
  
  infos1 <- Scissor(as(bulk_dataset, "matrix"), sc_dataset, phenotype, 
                    alpha = 0, 
                    tag = names(table(meta_bulk[, condition_label])),
                    family = "binomial", 
                    Save_file = gsub(".csv", ".RData", file.path(results_dir, results_file_name)))
  
  if (profile_mem) {
    end_time <- Sys.time()
    Rprof(NULL)
    mem <- max(summaryRprof("scissor_mem.out", chunksize=20000, memory="both")$by.total$mem.total)
    saveRDS(c(time = end_time - start_time, mem = mem),
            file = file.path(results_dir, gsub("_res.csv", "_time.rds", results_file_name)))
  }
  
  if (perm) {

    perm_res <- BiocParallel::bplapply(1:nperm, function(x) {
      infos1 <- Scissor(as(bulk_dataset, "matrix"), sc_dataset, sample(phenotype), 
                        alpha = 0, 
                        tag = names(table(meta_bulk[, condition_label])),
                        family = "binomial", 
                        Save_file = gsub(".csv", ".RData", file.path(results_dir, "tmp.csv")))$Coefs
    }, BPPARAM = BPPARAM)
    
    
  }
  
  
  write.csv(do.call(rbind, perm_res),
            file = file.path(results_dir, gsub(".csv", "_perm.csv", results_file_name)), 
            row.names = FALSE)
  
  write.csv(infos1$Coefs, file = file.path(results_dir, results_file_name), row.names = FALSE)
  return(infos1)
}

