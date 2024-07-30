import os
import numpy as np
import anndata
import oor_benchmark
import scanpy as sc
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from oor_benchmark.methods import scArches_milo
from oor_benchmark.methods import scArches_mappingQC
from oor_benchmark.methods import scVI_milo
import pandas as pd
import re



parser = ArgumentParser()
parser.add_argument("--adata_file", type=str, help="input h5ad file paths")
parser.add_argument("--results_dir", type=str, help="results output dir")
parser.add_argument("--condition_label", type=str, help="data set name")
parser.add_argument("--sample_label", type=str, help="sample name")
parser.add_argument("--results_file_name", type=str, help="results file")
parser.add_argument("--ref_label", type=str, help="Control")
parser.add_argument("--target_label", type=str, help="Treatment")

args = parser.parse_args()
print(args)

adata_file = args.adata_file
results_dir = args.results_dir
condition_label = args.condition_label
results_file_name = args.results_file_name
ref_label = args.ref_label
target_label = args.target_label
sample_label = args.sample_label


def run_CR(adata_file, condition_label="Condition", 
            results_dir="results",
            sample_label="",
            results_file_name="CR_output.csv",
            ref_label="", target_label=""):

            adata = anndata.read_h5ad(adata_file)
            adata.obs["dataset_group"] = adata.obs[condition_label].replace({ref_label: "ctrl", target_label: "query"})
            adata.obs["sample_id"] = adata.obs[sample_label]
            adata.obs["cell_annotation"] = adata.obs["major.celltype"]

            adata_scArches = scArches_milo.scArches_ctrl_milo_ctrl(adata, harmonize_output=True, outdir = "tmp")
            #adata_scArches = scArches_mappingQC.scArches_ctrl_mappingQClabels(adata_scArches, harmonize_output=True, outdir= "tmp")
            sc.pp.neighbors(adata_scArches, use_rep="X_scVI")
            sc.tl.umap(adata_scArches)

            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5), gridspec_kw={'wspace':0.7})
            sc.pl.umap(adata_scArches, color='mappingQC_labels', ax=ax1, title="mappingQC_labels", show=False)
            sc.pl.umap(adata_scArches, color='cell_annotation', ax=ax2, title="cell_annotation", show=False)
            sc.pl.umap(adata_scArches, color='dataset_group', ax=ax3, title="dataset_group", show=False)

            fig.savefig(results_dir + '/' + re.sub(".csv", "_scarches_umap_plots.pdf", results_file_name))


            adata_scvi = scVI_milo.scVI_ctrl_milo_ctrl(adata, harmonize_output=True, outdir = "tmp")
            #adata_scvi = scArches_mappingQC.scArches_ctrl_mappingQClabels(adata_scvi, harmonize_output=True, outdir= "tmp")
            sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
            sc.tl.umap(adata_scvi)

            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5), gridspec_kw={'wspace':0.7})
            sc.pl.umap(adata_scvi, color='mappingQC_labels', ax=ax1, title="mappingQC_labels", show=False)
            sc.pl.umap(adata_scvi, color='cell_annotation', ax=ax2, title="cell_annotation", show=False)
            sc.pl.umap(adata_scvi, color='dataset_group', ax=ax3, title="dataset_group", show=False)

            fig.savefig(results_dir + '/' + re.sub(".csv", "_scvi_umap_plots.pdf", results_file_name))


            scArches_rep = adata_scArches.obsm["X_scVI"]
            scVI_rep = adata_scvi.obsm["X_scVI"]

            pd.DataFrame(scArches_rep).to_csv(results_dir + '/' + re.sub(".csv",
                                                    "_scArches_rep.csv", results_file_name), index=False)
            pd.DataFrame(scVI_rep).to_csv(results_dir + '/' + re.sub(".csv",
                                                    "_scVI_rep.csv", results_file_name), index=False)



run_CR(adata_file, condition_label=condition_label, 
results_dir=results_dir,
sample_label=sample_label,
results_file_name=results_file_name,
ref_label=ref_label, target_label=target_label)                                                           
