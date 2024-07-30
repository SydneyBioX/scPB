import sys
import random
from argparse import ArgumentParser
import anndata as adata
import scanpy as sc
import multianndata as mad
import os
import cna
import numpy as np
import gzip
import pandas as pd
import sys
import re
import multiprocessing
import tracemalloc
import time

parser = ArgumentParser()
parser.add_argument("--adata_file", type=str, help="input h5ad file paths")
parser.add_argument("--results_dir", type=str, help="results output dir")
parser.add_argument("--condition_label", type=str, help="data set name")
parser.add_argument("--sample_label", type=str, help="data set name")
parser.add_argument("--results_file_name", type=str, help="data set name")
parser.add_argument("--ref_label", type=str, help="Control")
parser.add_argument("--target_label", type=str, help="Treatment")
parser.add_argument("--n_neighbors", type=int, default=15)
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
sample_label = args.sample_label
ref_label = args.ref_label
target_label = args.target_label
n_neighbors = args.n_neighbors
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

cna_results_dir = results_dir + "/CNA_perm"



def perm_cna(d, sample_label, i, cna_results_dir):
    np.random.seed(i)
    perm_idx = np.random.permutation(range(d.shape[0]))

    d.obs["Condition"] = np.array(d.obs["Condition"])[perm_idx]
    d.obs[sample_label] = np.array(d.obs[sample_label])[perm_idx]
    d.obs_to_sample(['Condition'])
    res = cna.tl.association(d, d.samplem.Condition, force_recompute = True)
    return res.ncorrs

def run_CNA(data_file, condition_label="Condition", sample_label="Sample",
            results_dir="results",
            results_file_name="CNA_output.csv",
            ref_label="", target_label="",
            perm=False,
            nperm=100,
            num_threads=8,
            n_neighbors=15,
            subset = 1,
            hvg=None,
            profile_mem=False):

    data = adata.read_h5ad(adata_file)
    print(data)


    if subset < 1:
      n_cells = data.X.shape[0]
      selected_cols = np.random.choice(n_cells, size=round(n_cells * subset), replace=False)
      data = data[selected_cols]
    
    if hvg is not None:
      sc.pp.highly_variable_genes(data, n_top_genes = hvg)
      data = data[:,data.var.highly_variable]

    # we can create a MultiAnnData object from these
    d = mad.MultiAnnData(X=data.X,            # expression matrix
                         obs=data.obs,
                         sampleid=sample_label)        # cell metadata matrix
    d.obs["Condition"] = d.obs[condition_label].replace(
        {ref_label: 0, target_label: 1}).astype(int)
    d.obs_to_sample(['Condition'])
    


      
    sc.pp.neighbors(d, n_neighbors = n_neighbors)
    sc.tl.umap(d)
    if profile_mem:
      tracemalloc.start()
      start = time.time()
      
    #d.obs_to_sample(['Condition'])
    res = cna.tl.association(d, d.samplem.Condition)
    
    if profile_mem:
      end = time.time()
      current, peak = tracemalloc.get_traced_memory()
      print('\033[93m' + f"Current memory usage is {current / 10**3}KB; Peak was {peak / 10**3}KB; Diff = {(peak - current) / 10**3}KB" + "\033[0m")
      tracemalloc.stop()
      elapsed = end - start
      mem_res = np.array([elapsed, current/ 10**3, peak/ 10**3])
      pd.DataFrame(mem_res).to_csv(results_dir + '/' + re.sub("res.csv", "time.csv", results_file_name), index=False)
    
    output = res.ncorrs

    output = pd.DataFrame({'correlation': output})
    #output_perm = pd.DataFrame(res.nullncorrs)

    output.to_csv(results_dir + '/' + results_file_name, index=False)
    pd.DataFrame(selected_cols).to_csv(results_dir + '/' + re.sub(".csv", "_idx.csv", results_file_name), index=False)


    if perm:
        res_list = []

        # # Create a multiprocessing Pool with the desired number of processes
        with multiprocessing.Pool(processes=num_threads) as pool:
            results = pool.starmap(
                perm_cna, [(d, sample_label, i, cna_results_dir) for i in range(nperm)])

        # Collect the results
        res_list.extend(results)

        pd.DataFrame(res_list).to_csv(results_dir + '/' +
                                      re.sub(".csv", "_perm.csv", results_file_name), index=False)
        return output, res_list

    else:
        return output





run_CNA(adata_file, condition_label=condition_label, sample_label=sample_label,
        results_dir=results_dir,
        results_file_name=results_file_name,
        ref_label=ref_label,
        target_label=target_label,
        perm = perm,
        nperm = 1000,
        num_threads = 10,
        n_neighbors = n_neighbors,
        hvg = hvg,
        subset = subset,
        profile_mem=profile_mem)
