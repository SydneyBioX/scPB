library(scDesign3)
library(SingleCellExperiment)
sce <- readRDS("/dski/nobackup/yingxinl/melanoma/data/sce_melanoma_DMSO_patient_SCC16-0069.rds")
zero_exprs <- rowMeans(counts(sce) == 0)
sce <- sce[, sce$scClassify_tumour_prediction_coarse != "unassigned"]


melanoma_markers <- readxl::read_xlsx("/dski/nobackup/yingxinl/melanoma/mmc4-2.xlsx", skip = 1)
melanoma_markers <- split(melanoma_markers$Gene, melanoma_markers$Signature)
Nazarian_markers <- read.csv("/dski/nobackup/yingxinl/melanoma/results/Nazarian_mapk_signature.csv", header = FALSE)
Nazarian_markers <- intersect(Nazarian_markers$V1, rownames(sce))
#cc_markers <- Seurat::cc.genes.updated.2019
hallmarkList <- readRDS("/dski/nobackup/yingxinl/melanoma/results/hallmarkList.rds")
nc_genes <- c(hallmarkList$HALLMARK_FATTY_ACID_METABOLISM,
              hallmarkList$HALLMARK_BILE_ACID_METABOLISM,
              hallmarkList$HALLMARK_CHOLESTEROL_HOMEOSTASIS)
keep_genes <- intersect(c(unlist(melanoma_markers),
                          Nazarian_markers,
                          #unlist(cc_markers),
                          nc_genes),
                        rownames(sce))
sce <- sce[keep_genes, ]
gc(reset = TRUE)
set.seed(2023)
sce <- sce[, sample(ncol(sce), ncol(sce) * 0.2)]
library(scater)

# 
# sce <- logNormCounts(sce)
# sce <- runPCA(sce)
# sce <- runUMAP(sce)
# plotUMAP(sce, colour_by = "Condition")
# plotUMAP(sce, colour_by = "scClassify_tumour_prediction_coarse")

table(sce$scClassify_tumour_prediction_coarse)
table(sce$Condition)

data <- construct_data(
  sce = sce,
  assay_use = "counts",
  celltype = "scClassify_tumour_prediction_coarse",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = "Condition",
  corr_by = "scClassify_tumour_prediction_coarse"
)
saveRDS(data, "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/model/melanoma_20220410/scDesign3_data_SCC16-0069.rds")
n_cores <- 25


print("marginal")
marginal <- fit_marginal(
  data = data,
  predictor = "gene",
  mu_formula = "scClassify_tumour_prediction_coarse * Condition",
  sigma_formula = 1,
  family_use = "nb",
  n_cores = n_cores,
  usebam = FALSE
)

saveRDS(marginal, file = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/model/melanoma_20220410/scDesign3_marginal_nb_SCC16-0069.rds")
print("copula")
copula <- fit_copula(
  sce = sce,
  assay_use = "counts",
  marginal_list = marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = n_cores,
  new_covariate = NULL,
  input_data = data$dat
)

saveRDS(copula, file = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/model/melanoma_20220410/scDesign3_copula_gaussian_SCC16-0069.rds")
saveRDS(sce, file = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/model/melanoma_20220410/scDesign3_sce_real.rds")
print("copula")
copula <- fit_copula(
  sce = sce,
  assay_use = "counts",
  marginal_list = marginal,
  family_use = "nb",
  copula = "vine",
  n_cores = n_cores,
  new_covariate = NULL,
  input_data = data$dat
)

saveRDS(copula, file = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/model/melanoma_20220410/scDesign3_copula_vine_SCC16-0069.rds")