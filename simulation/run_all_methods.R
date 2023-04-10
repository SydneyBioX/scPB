source("run_scissor.R")
source("run_scAB.R")
source("run_miloR.R")
source("run_DEGAS.R")
source("run_DAseq.R")
source("run_cacoa.R")
# data_dir = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/data/sce_sim_test.rds"
# adata_file = gsub("\\.rds", "\\.h5ad", data_dir)
# results_dir = "results"
# dataset_name = "sim_test"
# condition_label = "Condition"
# sample_label = "Sample"
# celltype_label = "scClassify_tumour_prediction_coarse"
# num_per_sample = 10

run_all_method <- function(data_dir = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/data/sce_sim_test.rds",
                           adata_file = gsub("\\.rds", "\\.h5ad", data_dir),
                           results_dir = "results",
                           dataset_name = "sim_test",
                           condition_label = "Condition",
                           sample_label = "Sample",
                           celltype_label = "scClassify_tumour_prediction_coarse",
                           num_per_sample = 10,
                           celltype_prop = TRUE,
                           ref_label = "",
                           target_label = "",
                           run_methods = "all") {
  
  if (run_methods == "all" | grepl("cacoa", run_methods)) {
    
    print("=============== cacoa ===============")
    cacoa_res <- run_cacoa(data_dir, results_dir, 
                           results_file_name = paste(dataset_name, "cacoa_output.csv", sep = "_"),
                           condition_label,
                           sample_label,
                           celltype_label,
                           ref_label,
                           target_label)
    
  }
  
  if (run_methods == "all" | grepl("scissor", run_methods)) {
    
    print("=============== scissor ===============")
    scissor_res <- run_scissor(data_dir, results_dir,
                               results_file_name = paste(dataset_name, "scissor_output.csv", sep = "_"),
                               condition_label,
                               sample_label,
                               celltype_label,
                               celltype_prop = celltype_prop,
                               num_per_sample = num_per_sample)
  }
  
  
  if (run_methods == "all" | grepl("scAB", run_methods)) {
    
    print("=============== scAB ===============")
    scAB_res <- run_scAB(data_dir, results_dir,
                         results_file_name = paste(dataset_name, "scAB_output.csv", sep = "_"),
                         condition_label, sample_label,
                         celltype_label,
                         celltype_prop = celltype_prop,
                         num_per_sample = num_per_sample,
                         target_label = target_label,
                         ref_label = ref_label)
  }
  
  
  if (run_methods == "all" | grepl("DEGAS", run_methods)) {
    
    print("=============== DEGAS ===============")
    DEGAS_res <- run_DEGAS(data_dir, results_dir,
                           results_file_name = paste(dataset_name, "DEGAS_output.csv", sep = "_"),
                           condition_label, sample_label,
                           celltype_label,
                           celltype_prop = celltype_prop,
                           num_per_sample = num_per_sample,
                           target_label = target_label,
                           ref_label = ref_label)
  }
  
  
  if (run_methods == "all" | grepl("miloR", run_methods)) {
    
    print("=============== miloR ===============")
    miloR_res <- run_miloR(data_dir, results_dir,
                           results_file_name = paste(dataset_name, "miloR_output.csv", sep = "_"),
                           condition_label,
                           sample_label)
  }
  
  
  if (run_methods == "all" | grepl("DAseq", run_methods)) {
    
    print("=============== DAseq ===============")
    DAseq_res <- run_DAseq(data_dir, results_dir,
                           results_file_name = paste(dataset_name, "DAseq_output.csv", sep = "_"),
                           condition_label)
  }
  
  if (run_methods == "all" | grepl("CNA", run_methods)) {
    
    print("=============== CNA ===============")
    cna_arg <- paste("/usr/bin/python3", "run_CNA.py", "--adata_file", adata_file,
                     "--results_dir", results_dir,
                     "--condition_label", condition_label,
                     "--sample_label", sample_label,
                     "--results_file_name", paste(dataset_name, "CNA_output.csv", sep = "_"), 
                     "--ref_label", ref_label,
                     "--target_label", target_label,
                     sep = " ")
    system(cna_arg)
  }
  
  if (run_methods == "all" | grepl("Hidden", run_methods)) {
    
    print("=============== Hidden ===============")
    cna_arg <- paste("/usr/bin/python3", "run_Hidden.py", "--adata_file", adata_file,
                     "--results_dir", results_dir,
                     "--condition_label", condition_label,
                     "--results_file_name", paste(dataset_name, "Hidden_output.csv", sep = "_"), 
                     "--ref_label", ref_label,
                     "--target_label", target_label,
                     sep = " ")
    system(cna_arg)
  }
  
  if (run_methods == "all" | grepl("MELD", run_methods)) {
    
    print("=============== MELD ===============")
    meld_arg <- paste("/usr/bin/python3", "run_MELD.py", "--adata_file", adata_file,
                      "--results_dir", results_dir,
                      "--condition_label", condition_label,
                      "--results_file_name", paste(dataset_name, "MELD_output.csv", sep = "_"), sep = " ")
    system(meld_arg)
  }
  
  
  return(NULL)
  
  
}