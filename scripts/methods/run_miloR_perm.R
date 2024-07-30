library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)


perm_milo <- function(milo_object, sample_label, condition_label, reduced.dim) {
  
  meta <- data.frame(colData(milo_object))[, c(sample_label, condition_label)]
  meta_perm <- meta[sample(nrow(meta)), ]
  milo_object <- countCells(milo_object, meta.data = meta_perm, sample = sample_label)
  milo_design <- meta_perm[, c(sample_label, condition_label)]
  colnames(milo_design) <- c("Sample", "Condition")
  milo_design <- distinct(milo_design)
  rownames(milo_design) <- milo_design$Sample
  
  ## Testing
  da_results <- testNhoods(milo_object, design = ~ Condition, design.df = milo_design,
                           reduced.dim = reduced.dim)
  
  
  
  nhoods_mat <- milo_object@nhoods
  nhoods_mat <- apply(nhoods_mat, 2, function(x) x/sum(x))
  milo_scores <- da_results$logFC %*% t(nhoods_mat)
  return(milo_scores)
}
library(BiocParallel)


run_miloR <- function(data_dir, results_dir, 
                      results_file_name, condition_label,
                      sample_label,
                      perm = FALSE,
                      nperm = 1000,
                      BPPARAM = SerialParam(),
                      subset = 1,
                      top_hvg = NULL,
                      k = 20,
                      prop = 0.1,
                      reduced.dim = "PCA",
                      reduced.dim_dir = NULL,
                      profile_mem = FALSE) {
  
  #Create a Milo object
  sce_sim <- readRDS(data_dir) 
  colData(sce_sim)[, sample_label] <- as.character(colData(sce_sim)[, sample_label])
  if (subset < 1) {
    sce_sim <- sce_sim[, sample(ncol(sce_sim), round(ncol(sce_sim) * subset))]
  }
  
  if (!is.null(top_hvg)) {
    decomp <- scran::modelGeneVar(sce_sim, BPPARAM = BPPARAM)
    hvg <- scran::getTopHVGs(decomp, n = top_hvg)
    sce_sim <- runPCA(sce_sim, BPPARAM = BPPARAM, subset_row = hvg)
  }
  
  
  if (!is.null(reduced.dim_dir)) {
    reduced.dim_mat <- read.csv(reduced.dim_dir)
    reduced.dim_mat <- as.matrix(reduced.dim_mat)
    rownames(reduced.dim_mat) <- colnames(sce_sim)
    reducedDim(sce_sim, reduced.dim) <- reduced.dim_mat
  }
  
  if (profile_mem) {
    Rprof("milo_mem.out", memory.profiling=T, interval=0.02)
    start_time <- Sys.time()
  }
  
  milo_object <- Milo(sce_sim)
  
  
  milo_object <- buildGraph(milo_object, k = k, d = 30, BPPARAM = BPPARAM)
  milo_object <- makeNhoods(milo_object, prop = prop, k = k, d = 30, refined = TRUE)
  plotNhoodSizeHist(milo_object)
  
  
  milo_object <- countCells(milo_object, meta.data = data.frame(colData(milo_object)), sample = sample_label)
  
  ## Defining experimental design
  
  milo_design <- data.frame(colData(milo_object))[, c(sample_label, condition_label)]
  colnames(milo_design) <- c("Sample", "Condition")
  ## Convert batch info from integer to factor
  
  milo_design <- distinct(milo_design)
  rownames(milo_design) <- milo_design$Sample
  
  ## Computing neighbourhood connectivity
  
  milo_object <- calcNhoodDistance(milo_object, d=ncol(reducedDim(sce_sim, reduced.dim)), reduced.dim = reduced.dim)
  
  
  ## Testing
  
  da_results <- testNhoods(milo_object, design = ~ Condition, design.df = milo_design,
                           reduced.dim = reduced.dim)
  
  
  nhoods_mat <- milo_object@nhoods
  nhoods_mat <- apply(nhoods_mat, 2, function(x) x/sum(x))
  milo_scores <- da_results$logFC %*% t(nhoods_mat)
  
  if (profile_mem) {
    end_time <- Sys.time()
    Rprof(NULL)
    mem <- max(summaryRprof("milo_mem.out", chunksize=20000, memory="both")$by.total$mem.total)
    saveRDS(c(time = end_time - start_time, mem = mem),
            file = file.path(results_dir, gsub("_res.csv", "_time.rds", results_file_name)))
  }
  
  
  
  if (perm) {
    milo_scores_perm <- BiocParallel::bplapply(1:nperm, function(x) {
      perm_milo(milo_object = milo_object, sample_label, condition_label, reduced.dim)
    }, BPPARAM = BPPARAM)
    
  }
  
  pvals_cell <- apply(nhoods_mat, 1, function(x) {
    if (sum(x) != 0) {
      pvals <- da_results$SpatialFDR[x != 0]
      combine_pval <- metapod::combineParallelPValues(lapply(pvals, function(x) x[1]), method = "fisher")$p.value
    } else {
      combine_pval <- 1
    }
    
  })
  
  
  milo_scores <- milo_scores[1, ]
  
  saveRDS(list(da_results = da_results, nhoods_mat = nhoods_mat, pvals_cell = pvals_cell), 
          file = file.path(results_dir, gsub(".csv", ".rds", results_file_name)))
  write.csv(milo_scores,
            file = file.path(results_dir, results_file_name), row.names = FALSE)
  
  write.csv(do.call(rbind, milo_scores_perm),
            file = file.path(results_dir, gsub(".csv", "_perm.csv", results_file_name)), 
            row.names = FALSE)
  
  return(milo_scores)
}
