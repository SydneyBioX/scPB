
library(DAseq)
library(BiocNeighbors)
library(SingleCellExperiment)
library(BiocParallel)
library(scater)



run_DAseq <- function(data_dir, results_dir, 
                      results_file_name, 
                      condition_label,
                      k.vector = seq(50, 500, 50),
                      subset = 1,
                      nperm = 2,
                      BPPARAM = BiocParallel::SerialParam(),
                      top_hvg = NULL,
                      profile_mem = FALSE) {
  
  sce_sim <- readRDS(data_dir) 
  
  if (subset < 1) {
    idx <- sort(sample(ncol(sce_sim), round(ncol(sce_sim) * subset)))
    sce_sim <- sce_sim[, idx]
    
    write.csv(idx,
              file = file.path(results_dir, gsub(".csv", "_idx.csv", results_file_name)), row.names = FALSE)
  }
  
  X.label.info <- colData(sce_sim)
  condition_group <- names(table(colData(sce_sim)[, condition_label]))
  condition <- X.label.info[, condition_label]
  labels.1 <- colnames(sce_sim)[condition == condition_group[1]]
  labels.2 <- colnames(sce_sim)[condition == condition_group[2]]
  
  
  if (!is.null(top_hvg)) {
    decomp <- scran::modelGeneVar(sce_sim)
    hvg <- scran::getTopHVGs(decomp, n = top_hvg)
    sce_sim <- runPCA(sce_sim, BPPARAM = BPPARAM, subset_row = hvg)
    sce_sim <- runUMAP(sce_sim, dimred = "PCA", BPPARAM = BPPARAM)
  } 
  
  
  if (profile_mem) {
    Rprof("daseq_mem.out", memory.profiling=T, interval=0.02)
    start_time <- Sys.time()
  }
  
  da_cells <- getDAcells(
    X = reducedDim(sce_sim, "PCA"), # input merged dataset of interest after dimension reduction.
    cell.labels = colnames(sce_sim),
    labels.1 = labels.1,
    labels.2 = labels.2,
    k.vector = k.vector,
    n.rand = nperm,
    plot.embedding = reducedDim(sce_sim, "UMAP"),
    do.plot = FALSE,
    save.knn = FALSE,
    BPPARAM = BPPARAM
  )
  

  
  if (profile_mem) {
    end_time <- Sys.time()
    Rprof(NULL)
    mem <- max(summaryRprof("daseq_mem.out", chunksize=20000, memory="both")$by.total$mem.total)
    saveRDS(c(time = end_time - start_time, mem = mem),
            file = file.path(results_dir, gsub("_res.csv", "_time.rds", results_file_name)))
  }

  write.csv(da_cells$da.pred,
            file = file.path(results_dir, results_file_name), row.names = FALSE)
  
  write.csv(do.call(cbind, da_cells$rand.pred),
            file = file.path(results_dir, gsub(".csv", "_perm.csv", results_file_name)), row.names = FALSE)
  return(da_cells)
  
}
