library(Augur)
library(SingleCellExperiment)
sce_sim <- readRDS("data/sce_sim_test.rds")
sce_sim$cell_type <- sce_sim$scClassify_tumour_prediction_coarse
augur <- calculate_auc(sce_sim, label_col = "Condition", cell_type_col = "cell_type", n_threads = 1)
augur$AUC