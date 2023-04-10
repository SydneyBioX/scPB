set.seed(2023)
library(Scissor)
library(SingleCellExperiment)
library(optparse)
library(insight)
source("utils_sim.R")
source("run_scissor.R")


option_list <- list(
  make_option(c("--data_dir"), type = "character", default = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/data_20230306", 
              help = ""),
  make_option(c("--data_file"), type = "character", default = "sce_sim_test2.rds", 
              help = ""),
  make_option(c("--dataset_name"), type = "character", default = "sim_test2", 
              help = ""),
  make_option(c("--sample_label"), type = "character", default = "library_id", 
              help = ""),
  make_option(c("--condition_label"), type = "character", default = "Condition", 
              help = ""),
  make_option(c("--celltype_label"), type = "character", default = "cell_type", 
              help = ""),
  make_option(c("--results_dir"), type = "character", default = "results", 
              help = ""),
  make_option(c("--results_file_name"), type = "character", default = "scissor_output_tmp.csv", 
              help = ""),
  # make_option(c("--UMAP"), type = "logical", default = TRUE, 
  #             help = "run UMAP"),
  make_option(c("--n_repeat"), type = "numeric", default = 10, 
              help = "The number of repeatition"),
  make_option(c("--perturb_celltype"), type = "character", default = "Microglia", 
              help = "")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)




sample_label = opt$sample_label
condition_label = opt$condition_label
celltype_label = opt$celltype_label
results_dir = opt$results_dir
results_file_name = opt$results_file_name
data_dir = opt$data_dir
data_file = opt$data_file
dataset_name = opt$dataset_name
n_repeat = opt$n_repeat
perturb_celltype = opt$perturb_celltype
num_per_sample_list = c(seq(1, 5, len = 5), 10, 20, 50, 100)

# 
# sample_label = "library_id"
# condition_label = "Condition"
# celltype_label = "cell_type"
# results_dir = "results"
# results_file_name = "scissor_output_tmp.csv"
# data_dir = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/data_20230306"
# data_file = "sce_sim_test2.rds"
# dataset_name = "sim_test2"
# num_per_sample_list = c(seq(1, 5, len = 5), 10, 20, 50, 100)
# n_repeat = 10
# 



sce_sim <- readRDS(file.path(data_dir, data_file))


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




#####################################   Pure   ##########################################################
print_color("====== Running pseudobulk Pure ======\n", "blue")


library(parallel)
set.seed(2023)
scissor_res_list <- mclapply(1:n_repeat, function(j) {
  print(j)
  
  scissor_res <- list()
  for (i in 1:length(num_per_sample_list)) {
    set.seed(j)
    print(i)
    bulk_res <- constructBulkFromSC(sce_sim[, sce_sim$cell_type == perturb_celltype],
                                    sample_label = sample_label,
                                    condition_label = condition_label,
                                    celltype_label = celltype_label,
                                    celltype_prop = FALSE,
                                    num_per_sample = num_per_sample_list[i],
                                    prop = 0.5)
    
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
    scissor_res[[i]] <-infos1$Coefs
    
  }
  return(scissor_res)
}, mc.cores = n_repeat)

scissor_res_list <- lapply(scissor_res_list, function(x) {
  names(x) <- num_per_sample_list
  x
})

scissor_res_list <- lapply(seq_along(num_per_sample_list), function(x) {
  lapply(scissor_res_list, function(l) {
    l[[x]]
  })
})
library(pROC)
auc_res <- lapply(scissor_res_list, function(l) {
  auc <- lapply(l, function(x) {
    idx <- sce_sim$Condition != "UN"
    roc_object <- roc(ifelse(sce_sim[, idx]$cell_type == "Microglia", 1, 0), x[idx])
    as.numeric(roc_object$auc)
  })
  unlist(auc)
})
saveRDS(scissor_res_list, file = file.path(results_dir, paste0(dataset_name, "_scissor_output_vary_nBulk_pure.rds")))
#
#
# auc <- lapply(scissor_res, function(x) {
#   idx <- sce_sim$Condition != "UN"
#   roc_object <- roc(ifelse(sce_sim[, idx]$cell_type == "Microglia", 1, 0), x[idx])
#   as.numeric(roc_object$auc)
# })
# unlist(auc)

#####################################   Same   ##########################################################
print_color("====== Running pseudobulk Same ======\n", "blue")


library(parallel)
set.seed(2023)
scissor_res_list <- mclapply(1:n_repeat, function(j) {
  print(j)

  scissor_res <- list()
  for (i in 1:length(num_per_sample_list)) {
    set.seed(j)
    print(i)
    bulk_res <- constructBulkFromSC(sce_sim,
                                    sample_label = sample_label,
                                    condition_label = condition_label,
                                    celltype_label = celltype_label,
                                    celltype_prop = FALSE,
                                    num_per_sample = num_per_sample_list[i],
                                    prop = 0.5)

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
    scissor_res[[i]] <-infos1$Coefs

  }
  return(scissor_res)
}, mc.cores = n_repeat)

scissor_res_list <- lapply(scissor_res_list, function(x) {
  names(x) <- num_per_sample_list
  x
})

scissor_res_list <- lapply(seq_along(num_per_sample_list), function(x) {
  lapply(scissor_res_list, function(l) {
    l[[x]]
  })
})
library(pROC)
auc_res <- lapply(scissor_res_list, function(l) {
  auc <- lapply(l, function(x) {
    idx <- sce_sim$Condition != "UN"
    roc_object <- roc(ifelse(sce_sim[, idx]$cell_type == "Microglia", 1, 0), x[idx])
    as.numeric(roc_object$auc)
  })
  unlist(auc)
})
saveRDS(scissor_res_list, file = file.path(results_dir, paste0(dataset_name, "_scissor_output_vary_nBulk.rds")))
#
#
# auc <- lapply(scissor_res, function(x) {
#   idx <- sce_sim$Condition != "UN"
#   roc_object <- roc(ifelse(sce_sim[, idx]$cell_type == "Microglia", 1, 0), x[idx])
#   as.numeric(roc_object$auc)
# })
# unlist(auc)

#####################################   Random   ############################################################
print_color("====== Running pseudobulk Random ======\n", "blue")



set.seed(2023)
scissor_res_list <- mclapply(1:n_repeat, function(j) {
  print(j)

  scissor_res <- list()
  for (i in 1:length(num_per_sample_list)) {
    set.seed(j)
    print(i)
    bulk_res <- constructBulkFromSC(sce_sim,
                                    sample_label = sample_label,
                                    condition_label = condition_label,
                                    celltype_label = celltype_label,
                                    celltype_prop = TRUE,
                                    num_per_sample = num_per_sample_list[i],
                                    prop = 0.5)

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
    scissor_res[[i]] <-infos1$Coefs

  }
  return(scissor_res)
}, mc.cores = n_repeat)

scissor_res_list <- lapply(scissor_res_list, function(x) {
  names(x) <- num_per_sample_list
  x
})

scissor_res_list <- lapply(seq_along(num_per_sample_list), function(x) {
  lapply(scissor_res_list, function(l) {
    l[[x]]
  })
})
library(pROC)
auc_res <- lapply(scissor_res_list, function(l) {
  auc <- lapply(l, function(x) {
    idx <- sce_sim$Condition != "UN"
    roc_object <- roc(ifelse(sce_sim[, idx]$cell_type == "Microglia", 1, 0), x[idx])
    as.numeric(roc_object$auc)
  })
  unlist(auc)
})



saveRDS(scissor_res_list, file = file.path(results_dir, paste0(dataset_name, "_scissor_output_vary_nBulk_varyCellProp.rds")))




#####################################   Confound   ##########################################################
print_color("====== Running pseudobulk Confound ======\n", "blue")


set.seed(2023)
scissor_res_list <- mclapply(1:n_repeat, function(j) {
  print(j)

  scissor_res <- list()
  for (i in 1:length(num_per_sample_list)) {
    set.seed(j)
    print(i)

    bulk_res_UN <- constructBulkFromSC(sce_sim[, sce_sim$Condition == "Control"],
                                       sample_label = sample_label,
                                       condition_label = condition_label,
                                       celltype_label = celltype_label,
                                       celltype_prop = TRUE,
                                       num_per_sample = num_per_sample_list[i],
                                       prop = 0.5,
                                       ctProp_bulk = c(0.4, 0.1, 0.2, 0.3)
    )
    meta_bulk_UN <- bulk_res_UN$meta_bulk
    bulk_dataset_UN <- bulk_res_UN$bulk_dataset

    bulk_res_TM <- constructBulkFromSC(sce_sim[, sce_sim$Condition != "Control"],
                                       sample_label = sample_label,
                                       condition_label = condition_label,
                                       celltype_label = celltype_label,
                                       celltype_prop = TRUE,
                                       num_per_sample = num_per_sample_list[i],
                                       prop = 0.5,
                                       ctProp_bulk = c(0.1, 0.4, 0.2, 0.3)
    )
    meta_bulk_TM <- bulk_res_TM$meta_bulk
    bulk_dataset_TM <- bulk_res_TM$bulk_dataset

    bulk_dataset <- cbind(bulk_dataset_UN, bulk_dataset_TM)
    meta_bulk <- rbind(meta_bulk_UN, meta_bulk_TM)


    phenotype <- meta_bulk[, condition_label]
    names(phenotype) <- rownames(meta_bulk)
    phenotype <- as.numeric(as.factor(phenotype)) - 1



    infos1 <- Scissor(as(bulk_dataset, "matrix"), sc_dataset, phenotype,
                      alpha = 0.05,
                      tag = names(table(meta_bulk[, condition_label])),
                      family = "binomial",
                      Save_file = gsub(".csv", ".RData", file.path(results_dir, results_file_name)))
    scissor_res[[i]] <-infos1$Coefs

  }
  return(scissor_res)
}, mc.cores = n_repeat)

scissor_res_list <- lapply(scissor_res_list, function(x) {
  names(x) <- num_per_sample_list
  x
})

scissor_res_list <- lapply(seq_along(num_per_sample_list), function(x) {
  lapply(scissor_res_list, function(l) {
    l[[x]]
  })
})
library(pROC)
auc_res <- lapply(scissor_res_list, function(l) {
  auc <- lapply(l, function(x) {
    idx <- sce_sim$Condition != "UN"
    roc_object <- roc(ifelse(sce_sim[, idx]$cell_type == "Microglia", 1, 0), x[idx])
    as.numeric(roc_object$auc)
  })
  unlist(auc)
})
boxplot(do.call(cbind, auc_res))


saveRDS(scissor_res_list, file = file.path(results_dir, paste0(dataset_name, "_scissor_output_vary_nBulk_varyCellProp_confound.rds")))





