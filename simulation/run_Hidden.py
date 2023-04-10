import sys


sys.path.append('/albona/nobackup2/yingxinl/miniconda3/lib/python3.9/site-packages/')
print (sys.path)
RAND_SEED = 28
CASE_COND = 1 #define that case condition is 1 

import pandas as pd
import numpy as np
from scipy.special import logit
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from scipy.stats import describe as describe_stats
import scipy as sp
import anndata as adata
from argparse import ArgumentParser

def standalone_logistic(X, y):
    clf = LogisticRegression(random_state=RAND_SEED, penalty='none').fit(X, y)
    predicted_label = clf.predict(X)
    predicted_prob = clf.predict_proba(X)
    return predicted_prob[:,1]

def PCA_logistic_kmeans(adata, num_pcs):
    p_hat = standalone_logistic(adata.obsm['PCA'].iloc[:, 0:num_pcs], adata.obs['status'].values)
    df_p_hat = pd.DataFrame()
    df_p_hat['p_hat'] = p_hat
    df_p_hat['status'] = adata.obs['status'].values
    adata.obs['p_hat'] = p_hat
    logit_p_hat = logit(p_hat)
    df_p_hat['logit_p_hat'] = logit_p_hat
    adata.obs['logit_p_hat'] = logit_p_hat
    kmeans_case = KMeans(n_clusters=2, random_state=0).fit(pd.DataFrame(logit(p_hat))[(adata.obs['status']==CASE_COND).values])
    mean_p_hat_kmeans_label0 = np.mean(p_hat[(adata.obs['status']==CASE_COND).values][kmeans_case.labels_==0]) 
    mean_p_hat_kmeans_label1 = np.mean(p_hat[(adata.obs['status']==CASE_COND).values][kmeans_case.labels_==1])
    zero_lab_has_lower_mean = mean_p_hat_kmeans_label0 < mean_p_hat_kmeans_label1
    df_p_hat_clust_case = df_p_hat.copy()
    df_p_hat_clust_case['kmeans'] = 0
    df_p_hat_clust_case['kmeans'][(adata.obs['status']==CASE_COND).values] = [1 if x==int(zero_lab_has_lower_mean) else 0 for x in kmeans_case.labels_]
    return p_hat, df_p_hat_clust_case['kmeans']


def choose_optimal_PC(adata):
    num_pcs = 2
    num_pcs_vec = []
    ks_vec = []
    ks_pval_vec = []
    num_pcs_ks_dict = {'num_pcs':num_pcs_vec, 'ks':ks_vec, 'ks_pval':ks_pval_vec}
    for num_pcs in np.arange(2, 61, 1):
        try:
            adata_ks = adata.copy()
        
            adata_ks.obs['status'] = adata_ks.obs['status'].astype('int').values
            p_hat, new_labels = PCA_logistic_kmeans(adata_ks, num_pcs)
            adata_ks.obs['new_labels'] = new_labels.values
            
            conditions = [(adata_ks.obs['status']==0), 
                        (adata_ks.obs['status']==1) & (adata_ks.obs['new_labels']==0), 
                        (adata_ks.obs['status']==1) & (adata_ks.obs['new_labels']==1)]
            values = ['control', 'case_new0', 'case_new1']
            adata_ks.obs['three_labels_batch_newlabels'] = np.select(conditions, values)
            # Heuristic decision criterion for NUM_PCS: choose the NUM_PCS that maximizes the Kolmogorov-Smirnov test statistic 
            # comparing the sample distribution of p_hat for HiDDEN label 0 and HiDDEN label 1 within the case condition
            ks_stat, ks_pval = sp.stats.ks_2samp(adata_ks.obs['p_hat'][adata_ks.obs['three_labels_batch_newlabels'].isin(['control'])].values, 
                                            adata_ks.obs['p_hat'][adata_ks.obs['three_labels_batch_newlabels'].isin(['case_new1', 'case_new0'])].values, 
                                            alternative='greater')
                                                
            num_pcs_vec.append(num_pcs)
            ks_vec.append(ks_stat)
            ks_pval_vec.append(ks_pval)
            del adata_ks
        except:
            print(f"An error occurred while processing num_pcs={num_pcs}")
    return num_pcs_ks_dict



parser = ArgumentParser()
parser.add_argument("--adata_file", type=str, help="input h5ad file paths")
parser.add_argument("--results_dir", type=str, help="results output dir")
parser.add_argument("--condition_label", type=str, help="data set name")
parser.add_argument("--results_file_name", type=str, help="data set name")
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

def run_Hidden(data_file, condition_label = "Condition", results_dir = "results", 
results_file_name = "Hidden_output.csv", ref_label = "", target_label = ""):
  
  data = adata.read_h5ad(adata_file)
  data.obs['status'] = data.obs[condition_label].replace(
        {ref_label: 0, target_label: 1}).astype(int)
  num_pcs_ks_dict = choose_optimal_PC(data)
  optimal_pcs = num_pcs_ks_dict['num_pcs'][np.argmax(num_pcs_ks_dict['ks'])]
  
  # run with optimal PCs
  p_hat, new_labels = PCA_logistic_kmeans(data, 
                                          num_pcs=optimal_pcs)
  
  output = pd.DataFrame({'p_hat': p_hat, "new_labels" : new_labels})
  output.to_csv(results_dir + '/' + results_file_name, index=False)
                                          


run_Hidden(adata_file, condition_label = condition_label, 
          results_dir = results_dir, 
          results_file_name = results_file_name,
          ref_label = ref_label,
          target_label = target_label)

                                          
