#scanpy tutorial day 3
#load required libraries
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

#set up directory
large_root = r"C:\Users\bsparta\Desktop\scRNAseq_workshop"
results_file = large_root + "\pbmc3k.h5ad"

#load data
print("Reading data...")
adata = sc.read_10x_mtx(large_root, var_names = 'gene_symbols',
                        cache=True)

# annotate mitochondrial genes as 'mt' and calculate qc metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# actually do the filtering by slicing the AnnData object
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# normalize counts to 10,000 reads per cell..
sc.pp.normalize_total(adata, target_sum=1e4)

#log n+1 transform
sc.pp.log1p(adata)

# identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

#plot dispersion vs expression
sc.pl.highly_variable_genes(adata)

#filter to only include HVGs
adata.raw = adata
adata = adata[:, adata.var.highly_variable]

#scale to unit variance for PCA, clip values exceding var 10
sc.pp.scale(adata, max_value=10)

# reduce the dimensionality of the data with PCA
sc.tl.pca(adata, svd_solver='arpack')

#scatter in the PCA coordinates, won't be used later on, only for vis
sc.pl.pca(adata, color='CST3')

#save the result
adata.write(results_file)

#computing the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

#calculating umap coordinates
sc.tl.umap(adata)

############ day 3 new material

#clustering the neighborhood graph
# run leiden clustering
sc.tl.leiden(adata)

#plot the clusters
sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])

#save results
adata.write(results_file)

#compute PAGA graph
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color=['leiden', 'CST3', 'NKG7'])

#recompute embedding using PAGA initialization
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['leiden', 'CST3', 'NKG7'], legend_loc='on data')

#############

#find marker genes
#compute ranking for differential genes in each cluster
#use raw data
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

sc.settings.verbosity = 2  # reduce the verbosity

#compute using  wilcoxon rank-sum
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
adata.write(results_file)

#define list for reference
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']

adata = sc.read(results_file)

#show the results
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)

#get a table with scores and groups
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(5)

#compare to a single cluster

sc.tl.rank_genes_groups(adata, 'leiden', groups=['0'], reference='1', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['0'], n_genes=20)

#plot distributions with violin plot
sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)

#compare cluster 0 to rest of groups
adata = sc.read(results_file)
sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)

#compare a certain gene across groups
sc.pl.violin(adata, ['CST3', 'NKG7', 'PPBP'], groupby='leiden')

#label the cell types
new_cluster_names = [
    'extra cluster','CD4 T', 'CD14 Monocytes',
    'B', 'CD8 T',
    'NK', 'FCGR3A Monocytes',
    'Dendritic', 'Megakaryocytes']
adata.rename_categories('leiden', new_cluster_names)
sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False, save='.pdf')

#visualize the marker genes in a dotplot
sc.pl.dotplot(adata, marker_genes, groupby='leiden');

#visualize marker gene distributions in a dot plot
sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', rotation=90);

#save results
adata.write(results_file, compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading

# if you want to export features to csv

# Export single fields of the annotation of observations
# adata.obs[['n_counts', 'louvain_groups']].to_csv(
#     './write/pbmc3k_corrected_louvain_groups.csv')

# Export single columns of the multidimensional annotation
# adata.obsm.to_df()[['X_pca1', 'X_pca2']].to_csv(
#     './write/pbmc3k_corrected_X_pca.csv')

# Or export everything except the data using `.write_csvs`.
# Set `skip_data=False` if you also want to export the data.
# adata.write_csvs(results_file[:-5], )







