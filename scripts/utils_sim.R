
.libPaths("/gpfs/gibbs/project/zhao/yl2687/R/4.3")

#.libPaths("/dskh/nobackup/biostat/projects/singlecell/scPB/Rpackages")
gen_sim_composition_wPatEff <- function(data, marginal, copula, sce_real, 
                                        cell_type_label, sample_label, condition_label, n_per_celltype,
                                        ref_label = NULL, target_label = NULL,
                                        sigma_prop = 0.01,
                                        type = "discrete",
                                        start_celltype = NULL,
                                        end_celltype = NULL) {
  
  
  n_celltypes <- length(table(data$dat[, cell_type_label]))
  n_library_id <- length(table(data$dat[, sample_label]))
  
  
  
  meta <- unique(data$dat[, c(sample_label, condition_label)])
  colnames(n_per_celltype) <- unique(meta[, condition_label])
  
  n_per_sample <- lapply(meta$library_id, function(i) {
    # cnt <- n_per_celltype[, meta[ meta$library_id == i, "Condition"]]
    # res <- rmultinom(1,
    #                  size = sum(cnt),
    #                  prob = cnt/sum(cnt))
    cov_mat <- diag(rep(sigma_prop, n_celltypes))
    cnt <- n_per_celltype[, meta[ meta$library_id == i, "Condition"]]
    mu <- log(cnt/sum(cnt))
    probs <- mvtnorm::rmvnorm(1, mean = mu, sigma = cov_mat)
    probs <- exp(probs)/sum(exp(probs))
    res <- rmultinom(1,
                     size = sum(cnt),
                     prob = probs)
    
    rownames(res) <- names(table(data$dat[, cell_type_label]))
    return(res)
  })
  
  n_per_sample <- do.call(cbind, n_per_sample)
  
  colnames(n_per_sample) <- meta$library_id
  print(n_per_sample)
  
  
  new_covariate <- lapply(colnames(n_per_sample), function(x) {
    df <- data.frame(rep(x, sum(n_per_sample[, x])),
                     rep(names(n_per_sample[, x]), n_per_sample[, x]))
    colnames(df) <- c(sample_label, cell_type_label)
    df
  })
  
  new_covariate <- do.call(rbind, new_covariate)
  new_covariate[, condition_label] <- mapvalues(new_covariate[, sample_label],  
                                                from = unique(data$dat[, c(sample_label, condition_label)])$library_id, 
                                                to = unique(data$dat[, c("library_id", "Condition")])$Condition)
  
  new_covariate$corr_group <- new_covariate[, cell_type_label]
  new_covariate$Condition2 <- new_covariate[, condition_label] # Condition 2 is the real condition label
  new_covariate$library_id2 <- new_covariate[, sample_label]
  
  idx1 <- new_covariate$Condition2 %in% target_label
  new_covariate$Condition[idx1] <- ref_label
  mapping <- matrix(names(table(new_covariate[, sample_label])), ncol = 2)
  
  new_covariate[, sample_label][idx1] <- mapvalues(new_covariate[, sample_label][idx1],
                                                   from = mapping[, 2],
                                                   to = mapping[, 1])
  
  
  
  # Here we remove some data from the UN
  
  n_cores <- 10
  para <- extract_para(
    sce = sce_real,
    marginal_list = marginal,
    n_cores = n_cores,
    family_use = "nb",
    new_covariate = new_covariate
  )
  new_data <- scDesign3::simu_new(
    sce = sce_real,
    mean_mat = para$mean_mat,
    sigma_mat = para$sigma_mat,
    zero_mat = para$zero_mat,
    copula_list = copula$copula_list,
    family_use = "nb",
    input_data = data$dat,
    new_covariate = new_covariate,
    n_cores = n_cores)
  
  new_data <- normalizeCounts(new_data)
  sce_sim <- SingleCellExperiment(assay = list(logcounts = cbind(new_data)),
                                  colData = rbind(new_covariate))
  
  if (type == "trajectory") {
    
    set.seed(2023)
    sce_sim <- runPCA(sce_sim)
    sce_sim <- runUMAP(sce_sim, min_dist = 0.3, verbose = TRUE, dimred = "PCA")
    plotPCA(sce_sim, colour_by = cell_type_label)
    plotUMAP(sce_sim, colour_by = cell_type_label)
    sce_slingshot <- slingshot::slingshot(sce_sim, 
                                          reducedDim = "PCA")
    
    sce_sim$pseudotime <- sce_slingshot$slingPseudotime_1
    
    agg_time <- aggregate(sce_sim$pseudotime, list(colData(sce_sim)[, cell_type_label]), mean)
    if (agg_time[agg_time$Group.1 == end_celltype, "x"] < agg_time[agg_time$Group.1 == start_celltype, "x"]) {
      sce_sim$pseudotime <- max(sce_sim$pseudotime) - sce_sim$pseudotime
    }
    
    plotUMAP(sce_sim, colour_by = "pseudotime")
    
    # Normalise the pseudotime
    sce_sim$pseudotime_prob <- minMax(qnorm(rank(sce_sim$pseudotime)/(length(sce_sim$pseudotime) + 1)))
    plotUMAP(sce_sim, colour_by = "pseudotime_prob")
    
    cond_label <- rbinom(ncol(sce_sim), 1, sce_sim$pseudotime_prob)
    sce_sim$Condition2 <- ifelse(cond_label == 1, target_label, ref_label)
    print("Trajectory")
    print(table(sce_sim$Condition2))
  }
  
  
  colData(sce_sim)[, condition_label] <- sce_sim$Condition2
  
  
  if (type == "trajectory") {
    colData(sce_sim)[, sample_label] <- paste(sce_sim$Condition2, sample(n_library_id/2, ncol(sce_sim), replace = TRUE), sep = '_')
  } else {
    colData(sce_sim)[, sample_label] <- sce_sim$library_id2
  }
  
  # sce_sim <- logNormCounts(sce_sim)
  sce_sim <- runPCA(sce_sim)
  sce_sim <- runUMAP(sce_sim, min_dist = 0.3)
  
  
  colData(sce_sim)[, condition_label] <- ifelse(colData(sce_sim)[, condition_label] == ref_label, "Control", "Treatment")
  colData(sce_sim)[, condition_label] <- factor(colData(sce_sim)[, condition_label], levels = c("Control", "Treatment"))
  
  
  
  print(table(colData(sce_sim)[, condition_label], colData(sce_sim)[, cell_type_label]))
  
  return(sce_sim)
}



gen_sim_combine <- function(data, marginal, copula, sce_real, cell_type_label, 
                            sample_label, condition_label, n_per_celltype,
                            DS_prop = 0.05, DS_logFC = 1.2, sigma_prop = 0.01,
                            ref_label = NULL, target_label = NULL, target_celltype = NULL) {
  
  n_celltypes <- length(table(data$dat[, cell_type_label]))
  n_library_id <- length(table(data$dat[, sample_label]))
  
  
  meta <- unique(data$dat[, c(sample_label, condition_label)])
  colnames(n_per_celltype) <- unique(meta[, condition_label])
  
  
  
  
  n_per_sample <- lapply(meta$library_id, function(i) {
    # cnt <- n_per_celltype[, meta[ meta$library_id == i, "Condition"]]
    # res <- rmultinom(1,
    #                  size = sum(cnt),
    #                  prob = cnt/sum(cnt))
    cov_mat <- diag(rep(sigma_prop, n_celltypes))
    cnt <- n_per_celltype[, meta[ meta$library_id == i, "Condition"]]
    mu <- log(cnt/sum(cnt))
    probs <- mvtnorm::rmvnorm(1, mean = mu, sigma = cov_mat)
    probs <- exp(probs)/sum(exp(probs))
    res <- rmultinom(1,
                     size = sum(cnt),
                     prob = probs)
    
    rownames(res) <- names(table(data$dat[, cell_type_label]))
    return(res)
  })
  
  n_per_sample <- do.call(cbind, n_per_sample)
  
  colnames(n_per_sample) <- meta$library_id
  print(n_per_sample)
  
  
  
  new_covariate <- lapply(colnames(n_per_sample), function(x) {
    df <- data.frame(rep(x, sum(n_per_sample[, x])),
                     rep(names(n_per_sample[, x]), n_per_sample[, x]))
    colnames(df) <- c(sample_label, cell_type_label)
    df
  })
  
  new_covariate <- do.call(rbind, new_covariate)
  new_covariate[, condition_label] <- mapvalues(new_covariate[, sample_label],  
                                                from = unique(data$dat[, c(sample_label, condition_label)])$library_id, 
                                                to = unique(data$dat[, c("library_id", "Condition")])$Condition)
  
  new_covariate$corr_group <- new_covariate[, cell_type_label]
  new_covariate$Condition2 <- new_covariate[, condition_label] # Condition 2 is the real condition label
  new_covariate$library_id2 <- new_covariate[, sample_label]
  
  
  idx1 <- new_covariate$Condition2 %in% target_label
  new_covariate$Condition[idx1] <- ref_label
  mapping <- matrix(names(table(new_covariate[, sample_label])), ncol = 2)
  
  new_covariate[, sample_label][idx1] <- mapvalues(new_covariate[, sample_label][idx1],
                                                   from = mapping[, 2],
                                                   to = mapping[, 1])
  
  n_genes <- nrow(sce_real)
  DS_idx <- sort(sample(n_genes, round(n_genes * DS_prop)))
  n_DS <- length(DS_idx)
  print(n_DS)
  lfc <- rep(0, n_genes)
  signs <- sample(c(-1, 1), n_DS, TRUE)
  lfc[DS_idx] <- rlnorm(n_DS, log(DS_logFC), 0.43) * signs
  fc <- 2^lfc
  
  # Here we remove some data from the UN
  
  n_cores <- 10
  
  para <- extract_para(
    sce = sce_real,
    marginal_list = marginal,
    n_cores = n_cores,
    family_use = "nb",
    new_covariate = new_covariate
    # data = data$dat
  )
  
  mean_mat <- para$mean_mat
  fc_mat <- matrix(1, ncol = ncol(mean_mat), nrow = nrow(mean_mat))
  DS_cell_idx <- which(new_covariate$Condition2 == target_label & new_covariate[, cell_type_label] == target_celltype)
  
  fc_mat[DS_cell_idx, ] <- matrix(rep(fc, length(DS_cell_idx)), nrow = length(DS_cell_idx), byrow = TRUE)
  mean_mat <- mean_mat * fc_mat
  
  
  new_data <- scDesign3::simu_new(
    sce = sce_real,
    mean_mat = mean_mat,
    sigma_mat = para$sigma_mat,
    zero_mat = para$zero_mat,
    copula_list = copula$copula_list,
    family_use = "nb",
    input_data = data$dat,
    new_covariate = new_covariate,
    n_cores = n_cores)
  
  new_data <- normalizeCounts(new_data)
  sce_sim <- SingleCellExperiment(assay = list(logcounts = cbind(new_data)),
                                  colData = rbind(new_covariate))
  colData(sce_sim)[, condition_label] <- sce_sim$Condition2
  colData(sce_sim)[, sample_label] <- sce_sim$library_id2
  # sce_sim <- logNormCounts(sce_sim)
  sce_sim <- runPCA(sce_sim)
  sce_sim <- runUMAP(sce_sim, min_dist = 0.3)
  colData(sce_sim)[, condition_label] <- ifelse(colData(sce_sim)[, condition_label] == ref_label, "Control", "Treatment")
  colData(sce_sim)[, condition_label] <- factor(colData(sce_sim)[, condition_label], levels = c("Control", "Treatment"))
  
  print(table(colData(sce_sim)[, condition_label], colData(sce_sim)[, cell_type_label]))
  
  return(sce_sim)
}










gen_sim_expression <- function(data, marginal, copula, sce_real, cell_type_label, 
                               sample_label, condition_label, n_per_celltype,
                               DS_prop = 0.05, DS_logFC = 1.2,
                               ref_label = NULL, target_label = NULL, target_celltype = NULL) {
  
  n_celltypes <- length(table(data$dat[, cell_type_label]))
  n_library_id <- length(table(data$dat[, sample_label]))
  
  
  meta <- unique(data$dat[, c(sample_label, condition_label)])
  colnames(n_per_celltype) <- unique(meta[, condition_label])
  new_covariate <- lapply(colnames(n_per_celltype), function(x) {
    df <- data.frame(rep(meta[meta$Condition == x, 1], sum(n_per_celltype[, x])),
                     rep(rep(names(table(data$dat[, cell_type_label])), n_per_celltype[, x]), 
                         length(meta[meta[, condition_label] == x, 1])))
    colnames(df) <- c(sample_label, cell_type_label)
    df
  })
  
  new_covariate <- do.call(rbind, new_covariate)
  new_covariate[, condition_label] <- mapvalues(new_covariate[, sample_label],  
                                                from = unique(data$dat[, c(sample_label, condition_label)])$library_id, 
                                                to = unique(data$dat[, c(sample_label, condition_label)])$Condition)
  
  new_covariate$corr_group <- new_covariate[, cell_type_label]
  new_covariate$Condition2 <- new_covariate[, condition_label] # Condition 2 is the real condition label
  new_covariate$library_id2 <- new_covariate[, sample_label]
  
  idx1 <- new_covariate$Condition2 %in% target_label
  new_covariate$Condition[idx1] <- ref_label
  mapping <- matrix(names(table(new_covariate[, sample_label])), ncol = 2)
  
  new_covariate[, sample_label][idx1] <- mapvalues(new_covariate[, sample_label][idx1],
                                                   from = mapping[, 2],
                                                   to = mapping[, 1])
  
  n_genes <- nrow(sce_real)
  DS_idx <- sort(sample(n_genes, round(n_genes * DS_prop)))
  n_DS <- length(DS_idx)
  print(n_DS)
  lfc <- rep(0, n_genes)
  signs <- sample(c(-1, 1), n_DS, TRUE)
  lfc[DS_idx] <- rlnorm(n_DS, log(DS_logFC), 0.43) * signs
  fc <- 2^lfc
  
  # Here we remove some data from the UN
  
  n_cores <- 10
  
  para <- extract_para(
    sce = sce_real,
    marginal_list = marginal,
    n_cores = n_cores,
    family_use = "nb",
    new_covariate = new_covariate
    # data = data$dat
  )
  
  mean_mat <- para$mean_mat
  fc_mat <- matrix(1, ncol = ncol(mean_mat), nrow = nrow(mean_mat))
  DS_cell_idx <- which(new_covariate$Condition2 == target_label & new_covariate[, cell_type_label] == target_celltype)
  
  fc_mat[DS_cell_idx, ] <- matrix(rep(fc, length(DS_cell_idx)), nrow = length(DS_cell_idx), byrow = TRUE)
  mean_mat <- mean_mat * fc_mat
  
  
  new_data <- scDesign3::simu_new(
    sce = sce_real,
    mean_mat = mean_mat,
    sigma_mat = para$sigma_mat,
    zero_mat = para$zero_mat,
    copula_list = copula$copula_list,
    family_use = "nb",
    input_data = data$dat,
    new_covariate = new_covariate,
    n_cores = n_cores)
  
  new_data <- normalizeCounts(new_data)
  sce_sim <- SingleCellExperiment(assay = list(logcounts = cbind(new_data)),
                                  colData = rbind(new_covariate))
  colData(sce_sim)[, condition_label] <- sce_sim$Condition2
  colData(sce_sim)[, sample_label] <- sce_sim$library_id2
  # sce_sim <- logNormCounts(sce_sim)
  sce_sim <- runPCA(sce_sim)
  sce_sim <- runUMAP(sce_sim, min_dist = 0.3)
  colData(sce_sim)[, condition_label] <- ifelse(colData(sce_sim)[, condition_label] == ref_label, "Control", "Treatment")
  colData(sce_sim)[, condition_label] <- factor(colData(sce_sim)[, condition_label], levels = c("Control", "Treatment"))
  
  print(table(colData(sce_sim)[, condition_label], colData(sce_sim)[, cell_type_label]))
  
  return(sce_sim)
}













plotSim <- function(sce_sim, cell_type_label, condition_label, sample_label) {
  df_toPlot <- moon::makeMoonDF(sce_sim)
  
  g4 <- ggplot(df_toPlot, aes(x = df_toPlot[, condition_label], fill = df_toPlot[, cell_type_label])) +
    geom_bar() +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +coord_flip() + theme(aspect.ratio = 0.2)
  
  g5 <- ggplot(df_toPlot, aes(x = df_toPlot[, sample_label], fill = df_toPlot[, cell_type_label])) +
    geom_bar() +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +coord_flip() + theme(aspect.ratio = 0.2)
  
  g1 <- ggplot(df_toPlot, aes(x = UMAP1, y = UMAP2, color = df_toPlot[, condition_label])) +
    geom_point(alpha = 0.5) + theme_bw() +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette = "Set2") 
  
  g2 <- ggplot(df_toPlot, aes(x = UMAP1, y = UMAP2, color = df_toPlot[, cell_type_label])) +
    geom_point(alpha = 0.5) +   theme_bw() +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette = "Set1") 
  
  
  g3 <- ggplot(df_toPlot, aes(x = UMAP1, y = UMAP2, color = df_toPlot[, sample_label])) +
    geom_point(alpha = 0.5) + theme_bw() +
    theme(aspect.ratio = 1) +
    scale_color_brewer(palette = "Dark2") 
  
  g1combine <- ggarrange(g1, g2, g3, ncol = 3, nrow = 1, align = "hv")
  g2combine <- ggarrange(g4, g5, ncol = 1, nrow = 2, align = "hv")
  g <- ggarrange(g1combine, g2combine, ncol = 1, nrow = 2, align = "hv")
  return(g)
}




getAggMat <- function(mat, celltype) {
  cellType_tab <- table(droplevels(factor(celltype)))
  cellTypes_n_mat <- matrix(rep(cellType_tab, 
                                nrow(mat)), 
                            nrow = length(cellType_tab), byrow = FALSE)
  
  aggMat <- Matrix.utils::aggregate.Matrix(t(as(mat, "dgCMatrix")),
                                           groupings = celltype,  fun = "sum")
  
  aggMat <- aggMat/cellTypes_n_mat
  return(aggMat)
}

rSampleProp <- function(num_celltypes){
  samp <- runif(num_celltypes)
  samp/sum(samp)
}

constructBulkFromSC <- function(sce_sim, sample_label, condition_label, celltype_label = NULL,
                                num_per_sample = 1, prop = 1, celltype_prop = FALSE, 
                                #ct_diff = FALSE,
                                ctProp_bulk = NULL) {
  
  if (num_per_sample == 1) {
    prop = 1
    celltype_prop = FALSE
  }
  samples <- colData(sce_sim)[, sample_label]
  
  if (celltype_prop) {
    celltype <- colData(sce_sim)[, celltype_label]
    celltype_group <- unique(celltype)
    if (is.null(ctProp_bulk)) {
      ctProp_list <- lapply(seq_len(num_per_sample), function(x) rSampleProp(length(celltype_group)))
    } else {
      ctProp_list <- lapply(seq_len(num_per_sample), function(x) ctProp_bulk)
    }
    
    
  } else {
    ctProp_list <- NULL
  }
  
  
  
  
  bulk_dataset <- lapply(seq_len(num_per_sample), function(x) {
    if (celltype_prop) {
      # Sample cell type proportion for each pseudo-bulk
      ctProp <- ctProp_list[[x]]
      idx <- sapply(1:length(celltype_group), function(i) {
        sample(which(celltype == celltype_group[i]), sum(celltype == celltype_group[i])*ctProp[i])
      })
      idx <- sort(unlist(idx))
    } else {
      idx <- sort(sample(ncol(sce_sim), ncol(sce_sim) * prop))
    }
    res <- getAggMat(logcounts(sce_sim)[, idx], samples[idx])
    res <- t(res)
  })
  bulk_dataset <- do.call(cbind, bulk_dataset)
  
  meta_bulk <- unique(colData(sce_sim)[, c(sample_label, condition_label)])
  rownames(meta_bulk) <- meta_bulk[, sample_label]
  meta_bulk <- meta_bulk[colnames(bulk_dataset), ]
  return(list(bulk_dataset = bulk_dataset, meta_bulk = meta_bulk, ctProp_list = ctProp_list))
}


readMethodResults <- function(sce_sim, results_dir, dataset_name) {
  # Meld
  
  meld_res <- read.csv(file.path(results_dir, paste(dataset_name, "MELD_output.csv", sep = "_")))
  dim(meld_res)
  sce_sim$MELD_res <- meld_res[, 2]
  
  
  cna_res <- read.csv(file.path(results_dir, paste(dataset_name, "CNA_output.csv", sep = "_")))
  dim(cna_res)
  sce_sim$CNA_res <- cna_res[, 1]
  
  
  scissor_res <- read.csv(file.path(results_dir, paste(dataset_name, "scissor_output.csv", sep = "_")))
  sce_sim$scissor_res <- scissor_res[, 1]
  
  
  scAB_res <- read.csv(file.path(results_dir, paste(dataset_name, "scAB_output.csv", sep = "_")))
  sce_sim$scAB_res <- scAB_res[, 1]
  
  
  
  DAseq_res <- read.csv(file.path(results_dir, paste(dataset_name, "DAseq_output.csv", sep = "_")))
  sce_sim$DAseq_res <- DAseq_res[, 1]
  
  
  miloR_res <- read.csv(file.path(results_dir, paste(dataset_name, "miloR_output.csv", sep = "_")))
  sce_sim$miloR_res <- miloR_res[, 1]
  
  
  
  DEGAS_res <- read.csv(file.path(results_dir, paste(dataset_name, "DEGAS_output.csv", sep = "_")))
  sce_sim$DEGAS_res <- DEGAS_res[, 2]
  
  
  cacoa_res <- read.csv(file.path(results_dir, paste(dataset_name, "cacoa_output.csv", sep = "_")))
  sce_sim$cacoa_res <- cacoa_res[, 2]
  
  hidden_res <- read.csv(file.path(results_dir, paste(dataset_name, "Hidden_output.csv", sep = "_")))
  sce_sim$hidden_res <- hidden_res[, 1]
  
  pencil_res <- read.csv(file.path(results_dir, paste(dataset_name, "PENCIL_output.csv", sep = "_")))
  sce_sim$pencil_res <- pencil_res[, 1]
  return(sce_sim)
}



minMax <- function(x) {
  (x - min(x))/(max(x) - min(x))
}

