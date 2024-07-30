import sys
import random
from argparse import ArgumentParser
import anndata as adata
import os
import meld
from sklearn.decomposition import PCA
import numpy as np
import gzip
import pandas as pd
import scanpy as sc

import re
import multiprocessing
from tqdm import tqdm
import tracemalloc
import time 


parser = ArgumentParser()
parser.add_argument("--adata_file", type=str, help="input h5ad file paths")
parser.add_argument("--results_dir", type=str, help="results output dir")
parser.add_argument("--condition_label", type=str, help="data set name")
parser.add_argument("--results_file_name", type=str, help="data set name")
parser.add_argument("--subset", type=float, default=1)
parser.add_argument("--hvg", type=int, default=None)
parser.add_argument("--nperm", type=int, default=1000)
parser.add_argument("--ncores", type=float, default=10)
parser.add_argument("--profile_mem", type=int, default=0)
parser.add_argument("--perm", type=int, default=1)
args = parser.parse_args()
print(args)

adata_file = args.adata_file
results_dir = args.results_dir
condition_label = args.condition_label
results_file_name = args.results_file_name
subset = args.subset
hvg = args.hvg
num_threads = args.ncores
nperm = args.nperm

if args.profile_mem == 1:
  profile_mem = True
if args.profile_mem == 0:
  profile_mem = False


if args.perm == 1:
  perm = True
if args.perm == 0:
  perm = False
  
# Define a helper function for processing each iteration
def perm_meld(condition, adata):
    perm_condition = np.random.permutation(condition)
    sample_densities = meld.MELD(verbose = 0).fit_transform(adata, perm_condition)
    sample_likelihoods_perm = meld.utils.normalize_densities(sample_densities)
    return np.array(sample_likelihoods_perm.iloc[:, 1])


  
def run_meld(adata_file, condition_label="Condition",
             results_dir="results",
             results_file_name="MELD_output.csv",
             perm = False,
             nperm = 100,
             num_threads = 8,
             subset = 1,
             hvg=None,
             profile_mem = False):
    # random.seed(2023)
    data = adata.read_h5ad(adata_file)
    print(data)

    
  
    if subset < 1:
      n_cells = data.X.shape[0]
      selected_cols = np.random.choice(n_cells, size=round(n_cells * subset), replace=False)
      data = data[selected_cols]
      pd.DataFrame(selected_cols).to_csv(results_dir + '/' + re.sub(".csv", "_idx.csv", results_file_name), index=False)
    
    if hvg is not None:
      sc.pp.highly_variable_genes(data, n_top_genes = hvg)
      data = data[:,data.var.highly_variable]


    condition = data.obs[condition_label]

    # Estimate density of each sample over the graph
    if profile_mem:
      tracemalloc.start()
      start = time.time()
      
    sample_densities = meld.MELD().fit_transform(data, condition)

    # Normalize densities to calculate sample likelihoods
    sample_likelihoods = meld.utils.normalize_densities(sample_densities)
    
    if profile_mem:
      end = time.time()
      current, peak = tracemalloc.get_traced_memory()
      print('\033[93m' + f"Current memory usage is {current / 10**3}KB; Peak was {peak / 10**3}KB; Diff = {(peak - current) / 10**3}KB" + "\033[0m")
      tracemalloc.stop()
      elapsed = end - start
      res = np.array([elapsed, current/ 10**3, peak/ 10**3])
      pd.DataFrame(res).to_csv(results_dir + '/' + re.sub("res.csv", "time.csv", results_file_name), index=False)
    
    
    sample_likelihoods.to_csv(results_dir + '/' + results_file_name, index=False)

    if perm:
        res_list = []
          


        # # Create a multiprocessing Pool with the desired number of processes
        with multiprocessing.Pool(processes=num_threads) as pool:
            results = list(tqdm(pool.starmap(perm_meld, [(condition, data, ) for _ in range(nperm)]), total=nperm))


        # Collect the results
        res_list.extend(results)

        
        pd.DataFrame(res_list).to_csv(results_dir + '/' + re.sub(".csv", "_perm.csv", results_file_name), index=False)
        return sample_likelihoods, res_list
    
    else:
        return sample_likelihoods



sample_likelihoods, perm_res = run_meld(adata_file, condition_label=condition_label, 
results_dir = results_dir,
results_file_name = results_file_name,
perm=perm, nperm = nperm, subset = subset, hvg = hvg, num_threads = int(num_threads),
profile_mem=profile_mem)
