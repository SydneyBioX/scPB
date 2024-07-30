import sys


from argparse import ArgumentParser
import anndata as adata
import scipy as sp
from scipy.stats import describe as describe_stats
from sklearn.linear_model import LogisticRegression
from sklearn.cluster import KMeans
from scipy.special import logit
import numpy as np
import pandas as pd
import sys
import re
import multiprocessing
import scanpy as sc
import tracemalloc
import time


module_path = 'codes/LabelCorrection//'

# Add the directory to the system path
sys.path.append(module_path)
import hiddensc


RAND_SEED = 28
CASE_COND = 1  # define that case condition is 1


parser = ArgumentParser()
parser.add_argument("--adata_file", type=str, help="input h5ad file paths")
parser.add_argument("--results_dir", type=str, help="results output dir")
parser.add_argument("--condition_label", type=str, help="data set name")
parser.add_argument("--results_file_name", type=str, help="data set name")
parser.add_argument("--ref_label", type=str, help="Control")
parser.add_argument("--target_label", type=str, help="Treatment")
parser.add_argument("--num_threads", type=int, help="Number of cores used")
parser.add_argument("--subset", type=float, default=1)
parser.add_argument("--hvg", type=int, default=None)
parser.add_argument("--profile_mem", type=int, default=0)
parser.add_argument("--perm", type=int, default=1)

args = parser.parse_args()
print(args)

adata_file = args.adata_file
results_dir = args.results_dir
condition_label = args.condition_label
results_file_name = args.results_file_name
ref_label = args.ref_label
target_label = args.target_label
num_threads = args.num_threads
subset = args.subset
hvg = args.hvg
if args.profile_mem == 1:
  profile_mem = True
if args.profile_mem == 0:
  profile_mem = False
  
if args.perm == 1:
  perm = True
if args.perm == 0:
  perm = False



def perm_hidden(data, optimal_pcs):
    y = np.random.permutation(data.obs['status'])
    x_pca = hiddensc.models.get_pca(data, n_comps=optimal_num_pcs_ks)
    p_hat, p_labels = hiddensc.models.logistic_predictions(x_pca, y, 1, rand_state=1)
    return p_hat


def run_Hidden(data_file, condition_label = "Condition", results_dir = "results",
               results_file_name = "Hidden_output.csv", 
               ref_label = "", target_label = "",
               perm = False,
               nperm = 100,
               num_threads = 8,
               subset = 1,
               hvg = None,
               profile_mem = False):

    data = adata.read_h5ad(adata_file)
    data.obs['status'] = data.obs[condition_label].replace(
        {ref_label: 0, target_label: 1}).astype(int)
    
    if subset < 1:
      n_cells = data.X.shape[0]
      selected_cols = np.random.choice(n_cells, size=round(n_cells * subset), replace=False)
      data = data[selected_cols]

    if hvg is not None:
        sc.pp.highly_variable_genes(data, n_top_genes = hvg)
        data = data[:,data.var.highly_variable]

    

    if profile_mem:
      tracemalloc.start()
      start = time.time()

    num_pcs, ks, ks_pval = hiddensc.models.determine_pcs_heuristic_ks(adata=data, orig_label="status", max_pcs=60)
    optimal_num_pcs_ks = num_pcs[np.argmax(ks)]
    feats = {}
    x_pca = hiddensc.models.get_pca(data, n_comps=optimal_num_pcs_ks)

    y = (data.obs['status'].values == 1).astype(np.int32)

    p_hat, p_labels = hiddensc.models.logistic_predictions(x_pca, data.obs['status'], 1, rand_state)

    output = pd.DataFrame({
    'p_hat': p_hat,
    'p_labels': p_labels})




    if profile_mem:
      end = time.time()
      current, peak = tracemalloc.get_traced_memory()
      print('\033[93m' + f"Current memory usage is {current / 10**3}KB; Peak was {peak / 10**3}KB; Diff = {(peak - current) / 10**3}KB" + "\033[0m")
      tracemalloc.stop()
      elapsed = end - start
      res = np.array([elapsed, current/ 10**3, peak/ 10**3])
      pd.DataFrame(res).to_csv(results_dir + '/' + re.sub("res.csv", "time.csv", results_file_name), index=False)
    

    output.to_csv(results_dir + '/' + results_file_name, index=False)

    if perm:
        res_list = []

        # # Create a multiprocessing Pool with the desired number of processes
        with multiprocessing.Pool(processes=num_threads) as pool:
            results = pool.starmap(
                perm_hidden, [(data, optimal_num_pcs_ks, ) for _ in range(nperm)])

        # Collect the results
        res_list.extend(results)

        pd.DataFrame(res_list).to_csv(results_dir + '/' +
                        re.sub(".csv", "_perm.csv", results_file_name), index=False)
        return output, res_list

    else:
        return output

run_Hidden(adata_file, condition_label=condition_label,
           results_dir=results_dir,
           results_file_name=results_file_name,
           ref_label=ref_label,
           target_label=target_label,
           perm=perm,
           nperm=1000,
           num_threads=num_threads,
           hvg = hvg,
           profile_mem=profile_mem)
