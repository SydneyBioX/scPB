library(inline)
openblas.set.num.threads <- cfunction( signature(ipt="integer"),
                                       body = 'openblas_set_num_threads(*ipt);',
                                       otherdefs = c ('extern void openblas_set_num_threads(int);'),
                                       libargs = c ('/usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3'),
                                       language = "C",
                                       convention = ".C"
)
openblas.set.num.threads(1)
library(SingleCellExperiment)
library(scater)

dir <- "/albona/nobackup2/biostat/datasets/singlecell/scPB/datasets/GSE165003_courtine_lumbar"
sce <- readRDS(file.path(dir, "sce_GSE165003_courtine_lumbar.rds"))
table(sce$cell_type, sce$library_id)


# We will subset the data for simulation
# Here, we will remove Neurons and Oligo

sce <- sce[, !sce$cell_type %in% c("Neurons", "Oligodendrocytes")]
sce



# Next, we will remove the genes that are not very expressed... (expressed at least 5% of cells)

sce <- sce[rowSums(counts(sce)) != 0, ]
sce


exprs_pct <- rowMeans(counts(sce) != 0)
summary((exprs_pct))
keep_genes <- names(exprs_pct[exprs_pct > 0.05])
length(keep_genes)

sce <- sce[keep_genes, ]

set.seed(2023)
sce <- runPCA(sce)
sce <- runUMAP(sce, min_dist = 0.3)

g1 <- plotUMAP(sce[, sample(ncol(sce))], colour_by = "cell_type")
g2 <- plotUMAP(sce[, sample(ncol(sce))], colour_by = "replicate")
g3 <- plotUMAP(sce[, sample(ncol(sce))], colour_by = "label")
g4 <- plotUMAP(sce[, sample(ncol(sce))], colour_by = "library_id")
library(ggpubr)
ggarrange(g1, g2, g3, g4, ncol = 2, nrow = 2, align = "hv")


sce$Condition <- sce$label

saveRDS(sce, "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/model/courtine_lumbar_scDesign3_sce.rds")

library(scDesign3)
library(SingleCellExperiment)
library(data.table)

gc(reset = TRUE)

library(scater)

# 
# sce <- logNormCounts(sce)
# sce <- runPCA(sce)
# sce <- runUMAP(sce)
# plotUMAP(sce, colour_by = "Condition")
# plotUMAP(sce, colour_by = "cell_type")

print(sce)
data <- construct_data(
  sce = sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = c("Condition", "library_id"),
  corr_by = "cell_type"
)
saveRDS(data, "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/model//courtine_lumbar_scDesign3_data.rds")
n_cores <- 25


print("marginal")
marginal <- fit_marginal(
  data = data,
  predictor = "gene",
  mu_formula = "cell_type * Condition + library_id",
  sigma_formula = 1,
  family_use = "nb",
  n_cores = n_cores,
  usebam = FALSE
)

saveRDS(marginal, file = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/model//courtine_lumbar_scDesign3_marginal_nb.rds")
print("copula")
copula <- fit_copula(
  sce = sce,
  assay_use = "counts",
  marginal_list = marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 1,
  new_covariate = NULL,
  input_data = data$dat
)

saveRDS(copula, file = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/model//courtine_lumbar_scDesign3_copula_gaussian.rds")
# 
print("copula")
copula <- fit_copula(
  sce = sce,
  assay_use = "counts",
  marginal_list = marginal,
  family_use = "nb",
  copula = "vine",
  n_cores = 1,
  new_covariate = NULL,
  input_data = data$dat
)

saveRDS(copula, file = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/model//courtine_lumbar_scDesign3_copula_vine.rds")
# 
# marginal <- readRDS("/dski/nobackup/yingxinl/melanoma/results/scDesign3_marginal_nb_courtine_lumbar.rds")
# copula <- readRDS("/dski/nobackup/yingxinl/melanoma/results/scDesign3_copula_vine_courtine_lumbar.rds")
# data <- readRDS("/dski/nobackup/yingxinl/melanoma/results/scDesign3_data_courtine_lumbar.rds")
# 
# new_covariate <- data.frame(Condition = rep(names(table(data$dat$Condition)), each = 1500),
#                             cell_type = rep(names(table(data$dat$cell_type)), 1000))
# n_cores <- 25
# para <- extract_para(
#   sce = sce,
#   marginal_list = marginal,
#   n_cores = n_cores,
#   family_use = "nb",
#   new_covariate = new_covariate
# )
# 
# 
# new_covariate$corr_group <- new_covariate$cell_type
# new_data <- scDesign3::simu_new(
#   sce = sce[keep_genes, ],
#   mean_mat = para$mean_mat,
#   sigma_mat = para$sigma_mat,
#   zero_mat = para$zero_mat,
#   copula_list = copula$copula_list,
#   family_use = "nb",
#   input_data = data$dat,
#   new_covariate = new_covariate,
#   n_cores = n_cores)
# 
# 
# library(scater)
# sce_sim <- SingleCellExperiment(assay = list(counts = new_data),
#                                 colData = new_covariate)
# 
# sce_sim <- logNormCounts(sce_sim)
# sce_sim <- runPCA(sce_sim)
# sce_sim <- runUMAP(sce_sim)
# 
# plotPCA(sce_sim, colour_by = "Condition")
# plotUMAP(sce_sim, colour_by = "cell_type")
# 
# g1 <- plotUMAP(sce[, sample(ncol(sce))], colour_by = "cell_type")
# g2 <- plotUMAP(sce[, sample(ncol(sce))], colour_by = "replicate")
# g3 <- plotUMAP(sce[, sample(ncol(sce))], colour_by = "label")
# g4 <- plotUMAP(sce[, sample(ncol(sce))], colour_by = "library_id")
# library(ggpubr)
# ggarrange(g1, g2, g3, g4, ncol = 2, nrow = 2, align = "hv")
# ggsave(file.path(dir, "UMAP.pdf"), width = 15, height = 12)
# saveRDS(sce, file = file.path(dir, "sce_GSE165003_courtine_lumbar.rds"))
# 

# saveRDS(sce, file = "data/K562_CRISPR_RNA_seq/sce_K562_CRISPR_conn.rds")
# # 
# 
# 
# n_cells_per_celltype <- 1000
# new_covariate <- data.frame(Condition = rep("DMSO", n_cells_per_celltype * 3),
#                             cell_type = rep(names(table(data$dat$cell_type)), n_cells_per_celltype))
# n_cores <- 25
# para <- extract_para(
#   sce = sce[keep_genes, ],
#   marginal_list = marginal,
#   n_cores = n_cores,
#   family_use = "nb",
#   new_covariate = new_covariate
# )
# 
# #data$corr_group <- 1
# new_covariate$corr_group <- new_covariate$cell_type
# new_data <- scDesign3::simu_new(
#   sce = sce[keep_genes, ],
#   mean_mat = para$mean_mat,
#   sigma_mat = para$sigma_mat,
#   zero_mat = para$zero_mat,
#   copula_list = copula$copula_list,
#   family_use = "nb",
#   input_data = data$dat,
#   new_covariate = new_covariate,
#   n_cores = n_cores)
# 
# # mean_eff <- ifelse(seq_len(ncol(para$mean_mat)) %in% nc_genes, 0, 1)
# # mean_eff <- matrix(rep(mean_eff, nrow(para$mean_mat)), nrow(para$mean_mat), ncol(para$mean_mat), byrow = TRUE)
# # new_data_DT <- scDesign3::simu_new(
# #   sce = sce[keep_genes, ],
# #   mean_mat = para$mean_mat + mean_eff,
# #   sigma_mat = para$sigma_mat,
# #   zero_mat = para$zero_mat,
# #   copula_list = copula$copula_list,
# #   family_use = "nb",
# #   input_data = data$dat,
# #   new_covariate = new_covariate,
# #   n_cores = n_cores)
# nc_genes <- sort(sample(ncol(para$mean_mat), 100))
# mean_eff <- ifelse(seq_len(ncol(para$mean_mat)) %in% nc_genes, 0, 1)
# mean_eff <- sapply(mean_eff, function(x) {
#   if (x == 0) {
#     rnorm(nrow(para$mean_mat), x, 0.1)
#   } else {
#     rnorm(nrow(para$mean_mat), x, 0.1)
#   }
# })
# new_data <- normalizeCounts(new_data)
# new_data_DT <- new_data + t(mean_eff)
# #new_data_DT[new_data_DT < 0] <- 0
# summary(rowMeans(new_data_DT)[nc_genes] - rowMeans(new_data)[nc_genes])
# summary(rowMeans(new_data_DT)[-nc_genes] - rowMeans(new_data)[-nc_genes])
# 
# library(scater)
# new_covariate_DT <- new_covariate
# new_covariate_DT$Condition <- "DT"
# sce_sim <- SingleCellExperiment(assay = list(logcounts = cbind(new_data, new_data_DT)),
#                                 colData = rbind(new_covariate, new_covariate_DT))
# 
# # sce_sim <- logNormCounts(sce_sim)
# sce_sim <- runPCA(sce_sim)
# sce_sim <- runUMAP(sce_sim)
# 
# 
# g1 <- plotUMAP(sce_sim, colour_by = "Condition") + theme(aspect.ratio = 1)
# g2 <- plotUMAP(sce_sim, colour_by = "cell_type") + theme(aspect.ratio = 1)
# 
# grid.arrange(g1, g2, ncol = 2)
# 
# library(reticulate)
# options(warn=1)
# ot <- import("ot")
# print(reticulate::py_config())
# 
# pca <- reducedDim(sce_sim, "PCA")
# pca_DMSO <- pca[sce_sim$Condition == "DMSO", ]
# pca_DT <- pca[sce_sim$Condition == "DT", ]
# 
# exprs_DMSO <- logcounts(sce_sim)[, sce_sim$Condition == "DMSO"]
# exprs_DT <- logcounts(sce_sim)[, sce_sim$Condition == "DT"]
# 
# 
# a = rep(1, nrow(pca_DMSO))
# b = rep(1, nrow(pca_DT))
# M = ot$dist(pca_DMSO, pca_DT)
# gamma = ot$emd(a, b, M)
# 
# 
# gamma_normalised <- apply(gamma, 1, function(x) x/sum(x))
# gamma_normalised <- t(gamma_normalised)
# exprs_gamma_DMSO <- exprs_DT %*% t(gamma_normalised)
# 
# rank_pct <- mean(diag(apply(gamma, 1, rank))/ncol(gamma))
# print(rank_pct)
# diff_gamma_DMSO <- (exprs_gamma_DMSO - exprs_DMSO)
# summary(rowMeans(exprs_gamma_DMSO - exprs_DMSO)[nc_genes])
# summary(rowMeans(exprs_gamma_DMSO - exprs_DMSO)[-nc_genes])
# 
# 
# pvalue <- apply(diff_gamma_DMSO, 1, function(x) t.test(x)$p.value)
# pvalue <- p.adjust(pvalue, method = "BH")
# pvalue_nc <- p.adjust(pvalue[nc_genes], method = "BH")
# print(mean(pvalue_nc < 0.05, na.rm = TRUE))
# print(mean(pvalue[-nc_genes] < 0.05, na.rm = TRUE))
# 
# 
# 
# 
# reg <- c(1, 10, 100, 1000)
# reg_m <- c(1, 10, 100, 1000, 5000, 10000)
# 
# grid_search <- expand.grid(reg, reg_m)
# 
# for (k in 1:nrow(grid_search)) {
#   print("-------------------------")
#   print(grid_search[k, ])
#   a = rep(1, nrow(pca_DMSO))
#   b = rep(1, nrow(pca_DT))
#   M = ot$dist(pca_DMSO, pca_DT)
#   gamma = ot$unbalanced$sinkhorn_unbalanced(a, b, M, 
#                                             reg = grid_search[k, 1], 
#                                             reg_m = grid_search[k, 2], 
#                                             log = FALSE)
#   gamma = ot$emd(a, b, M)
#   
#   
#   gamma_normalised <- apply(gamma, 1, function(x) x/sum(x))
#   gamma_normalised <- t(gamma_normalised)
#   exprs_gamma_DMSO <- exprs_DT %*% t(gamma_normalised)
#   
#   rank_pct <- mean(diag(apply(gamma, 1, rank))/ncol(gamma))
#   print(rank_pct)
#   diff_gamma_DMSO <- (exprs_gamma_DMSO - exprs_DMSO)
#   diff_gamma_DMSO
# 
#   
#   pvalue <- apply(diff_gamma_DMSO, 1, function(x) t.test(x)$p.value)
#   pvalue <- p.adjust(pvalue, method = "BH")
#   pvalue_nc <- p.adjust(pvalue[nc_genes], method = "BH")
#   print(mean(pvalue_nc < 0.05, na.rm = TRUE))
#   print(mean(pvalue[-nc_genes] < 0.05, na.rm = TRUE))
#   
#   
#   
# }


