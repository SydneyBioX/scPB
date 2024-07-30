import sys
sys.path.append('codes/PENCIL')

from pencil import *
import numpy as np
import anndata as adata
from argparse import ArgumentParser
import random
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import multiprocessing
import re
parser = ArgumentParser()
parser.add_argument("--adata_file", type=str, help="input h5ad file paths")
parser.add_argument("--results_dir", type=str, help="results output dir")
parser.add_argument("--condition_label", type=str, help="data set name")
parser.add_argument("--results_file_name", type=str, help="data set name")
parser.add_argument("--ref_label", type=str, help="Control")
parser.add_argument("--target_label", type=str, help="Treatment")
parser.add_argument("--perm", type=int, default=1)
parser.add_argument("--num_threads", type=int, help="Number of cores used")

args = parser.parse_args()
print(args)

adata_file = args.adata_file
results_dir = args.results_dir
condition_label = args.condition_label
results_file_name = args.results_file_name
ref_label = args.ref_label
target_label = args.target_label
num_threads = args.num_threads

if args.perm == 1:
  perm = True
if args.perm == 0:
  perm = False

multiprocessing.set_start_method('spawn', force=True)

def perm_pencil(expression_data, phenotype_labels, class_names, model):
    phenotype_labels = np.random.permutation(phenotype_labels)
    with mlflow.start_run():
        pred_labels, confidence, pred_labels_orig = model.fit_transform(
                expression_data, phenotype_labels,
                class_names = class_names,
                plot_show = False,
                use_cuda = True,
                #shuffle_rate=1/5,
                #lambda_L1=1e-5, 
                #lambda_L2=0.0,
                shuffle_rate=1/4,
                lambda_L1=1e-5, 
                lambda_L2=1e-3, 
                test = True,
                lr = 0.01,
                umap_plot = False
        )
    return pred_labels_orig[:, 1]



def run_PENCIL(data_file, condition_label = "Condition",
              results_dir = "results", 
              results_file_name = "PENCIL_output.csv",
              ref_label = "", target_label = "",
              perm = False,
              nperm = 100,
              num_threads = 8):
    
    random.seed(2023)
    data = adata.read_h5ad(data_file)
    print(data)


    # prepare data source
    sc.pp.highly_variable_genes(data, n_top_genes = 2000)
    data = data[:,data.var.highly_variable]
    expression_data = data.X
    phenotype_labels = data.obs[condition_label].replace(
            {ref_label: 0, target_label: 1}).astype(int)
    class_names = [ref_label, target_label]

    # init a pencil model
    model = Pencil(mode='multi-classification', select_genes=True, mlflow_record=True)

    # run

    with mlflow.start_run():
            pred_labels, confidence, pred_labels_orig = model.fit_transform(
                    expression_data, phenotype_labels,
                    class_names = class_names,
                    plot_show = False,
                    use_cuda = True,
                    #shuffle_rate=1/5,
                    #lambda_L1=1e-5, 
                    #lambda_L2=0.0,
                    shuffle_rate=1/4,
                    lambda_L1=1e-5, 
                    lambda_L2=1e-3, 
                    test = True,
                    lr = 0.01,
                    umap_plot = False
            )
        # gene_weights = model.gene_weights(plot=True)
    pred_labels = np.squeeze(pred_labels)
    confidence = np.squeeze(confidence)
    # res = pd.DataFrame(data = {"pred_labels": pred_labels[:,0], "confidence": confidence})
    res = pd.DataFrame(data = {"pred_labels": pred_labels, "confidence": confidence})
    res.to_csv(results_dir + '/' + results_file_name, index=False)

    if perm:
        res_list = []

        # # Create a multiprocessing Pool with the desired number of processes
        with multiprocessing.Pool(processes=num_threads) as pool:
            results = pool.starmap(
                perm_pencil, [(expression_data, phenotype_labels, class_names, model, ) for _ in range(nperm)])

        # Collect the results
        res_list.extend(results)
        print(res_list)

        pd.DataFrame(res_list).to_csv(results_dir + '/' +
                        re.sub(".csv", "_perm.csv", results_file_name), index=False)
        return res, res_list

    else:
        return res





run_PENCIL(adata_file, condition_label=condition_label, 
           results_dir=results_dir, 
           results_file_name=results_file_name,
           ref_label=ref_label, target_label=target_label,
           perm=True,
           nperm=100,
           num_threads=4)
