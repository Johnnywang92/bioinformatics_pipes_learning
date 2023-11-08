# scanpy tutorial 1

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

print('loading libs')

# verbosity: errors(0), warnings(1), hints(3)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=80, facecolor='white')

#set up directory
large_root = '/Users/yifanwang/workspace/single_cell_course/datasets/filtered_gene_bc_matrices/hg19'
result_file = large_root + 'pbmc3k.h5ad'

#load datasets
print('reading dataset...')
adata = sc.read_10x_mtx(large_root,var_names='gene_symbols', cache=True)

#print data
adata

#show genes that yield highest fraction of counts
sc.pl.highest_expr_genes(adata, n_top=20)

#annotate mitochondrial genes as 'mt' and calculate qc metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#plot a violin plot with some qc metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts','pct_counts_mt'], jitter=0.4,multi_panel=True)


# plot a histogram of a single variable
sns.histplot(adata.obs['n_genes_by_counts'])

# plot a scatter plot of QC metrics
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color='n_genes_by_counts')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt')

# save plot as png or pdf
ax = plt.gca()
plt.savefig(large_root + "/" + "QCbefore.png")

sns.jointplot(
    x = "n_genes_by_counts",
    y= "total_counts",
    data = adata.obs,
    kind = 'scatter',
    s=5
)

# save plot
ax = plt.gca()
plt.savefig(large_root + "/" + "QCjoint.png")

#basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# actually do the filtering by slicing the AnnData object
bdata = adata
bdata = bdata[bdata.obs.n_genes_by_counts < 2200,:]
bdata = bdata[bdata.obs.pct_counts_mt < 5,:]

# plot a scatter plot of filtered data

sc.pl.scatter(bdata, x='total_counts', y='pct_counts_mt', color='n_genes_by_counts')
sc.pl.scatter(bdata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt')

#save plot
ax = plt.gca()
plt.savefig(large_root + "/" + "QCafter.png")