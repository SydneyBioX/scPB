import sys


sys.path.append('/albona/nobackup2/yingxinl/miniconda3/lib/python3.9/site-packages/')
print (sys.path)
import pandas as pd
import gzip
import numpy as np
import cna
import os
import multianndata as mad
import scanpy as sc
import anndata as adata
from argparse import ArgumentParser
import random

parser = ArgumentParser()
parser.add_argument("--adata_file", type=str, help="input h5ad file paths")
parser.add_argument("--results_dir", type=str, help="results output dir")
parser.add_argument("--condition_label", type=str, help="data set name")
parser.add_argument("--sample_label", type=str, help="data set name")
parser.add_argument("--results_file_name", type=str, help="data set name")
parser.add_argument("--ref_label", type=str, help="Control")
parser.add_argument("--target_label", type=str, help="Treatment")

args = parser.parse_args()
print(args)

adata_file = args.adata_file
results_dir = args.results_dir
condition_label = args.condition_label
results_file_name = args.results_file_name
sample_label = args.sample_label
ref_label = args.ref_label
target_label = args.target_label

def run_CNA(data_file, condition_label = "Condition", sample_label = "Sample",
            results_dir = "results", 
            results_file_name = "CNA_output.csv",
            ref_label = "", target_label = ""):
    
    data = adata.read_h5ad(adata_file)
    print(data)
    # we can create a MultiAnnData object from these
    d = mad.MultiAnnData(X=data.X,            # expression matrix
                        obs=data.obs,
                        sampleid=sample_label)        # cell metadata matrix
    d.obs["Condition"] = d.obs[condition_label].replace(
        {ref_label: 0, target_label: 1}).astype(int)
    d.obs_to_sample(['Condition'])
    sc.pp.neighbors(d)
    sc.tl.umap(d) 
    d.obs_to_sample(['Condition'])
    res = cna.tl.association(d, d.samplem.Condition)

    output = res.ncorrs


    output = pd.DataFrame({'correlation': output})

    output.to_csv(results_dir + '/' + results_file_name, index=False)



# adata_file = "data/sce_sim_test.h5ad"
# condition_label = "Condition"
# results_dir = "results"
# results_file_name = "CNA_output.csv"


run_CNA(adata_file, condition_label = condition_label, sample_label = sample_label,
        results_dir = results_dir, 
        results_file_name = results_file_name,
        ref_label = ref_label,
        target_label = target_label)
