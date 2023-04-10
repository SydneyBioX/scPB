import sys
sys.path.append('/albona/nobackup2/yingxinl/miniconda3/lib/python3.9/site-packages/')
import pandas as pd
import gzip
import numpy as np
from sklearn.decomposition import PCA
import meld
import os
import anndata as adata
from argparse import ArgumentParser
import random


# os.chdir('/dski/nobackup/biostat/datasets/singlecell/scPB/datasets/MELD_dataset')
# filename = 'normalised_count.csv.gz'

# # Use pandas to read the csv.gz file
# df = pd.read_csv(gzip.open(filename), sep=',')

# # Do some processing on the data
# print(df.head())

# # make first column as the genes names
# df = df.set_index(df.columns[0])

# # df contains 19678 genes and 5754 cells 
# # perform library size normalization
# scale_factor = np.median(np.sum(df, axis=0))
# df = df / np.sum(df, axis=0) * scale_factor
# df = np.sqrt(df)

parser = ArgumentParser()
parser.add_argument("--adata_file", type=str, help="input h5ad file paths")
parser.add_argument("--results_dir", type=str, help="results output dir")
parser.add_argument("--condition_label", type=str, help="data set name")
parser.add_argument("--results_file_name", type=str, help="data set name")


args = parser.parse_args()
print(args)

adata_file = args.adata_file
results_dir = args.results_dir
condition_label = args.condition_label
results_file_name = args.results_file_name



def run_MELD (adata_file, condition_label = "Condition", 
              results_dir = "results", 
              results_file_name = "MELD_output.csv"):
    random.seed(2023)
    data = adata.read_h5ad(adata_file)
    print(data)

    # We will skip the size-factor normalisation above
    df = data.X
    # Perform PCA and extract the first 100 principal components
    pca = PCA(n_components=100)
    pcs = pca.fit_transform(df)   
    # pcs is a 6000 (cells) by 100 (PCs) matrix
    print(pcs.shape)
    condition = data.obs[condition_label]

    # Estimate density of each sample over the graph
    sample_densities = meld.MELD().fit_transform(pcs, condition)

    # Normalize densities to calculate sample likelihoods
    sample_likelihoods = meld.utils.normalize_densities(sample_densities)


    sample_likelihoods.to_csv(results_dir + '/' + results_file_name, index=False)




# adata_file = "data/sce_sim_test.h5ad"
# condition_label = "Condition"
# results_dir = "results"
# results_file_name = "MELD_output.csv"

run_MELD(adata_file, condition_label = condition_label, 
         results_dir = results_dir, 
         results_file_name = results_file_name)
