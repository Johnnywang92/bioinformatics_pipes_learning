#scanpy tutorial day 2

#remember to install kneed, sklearn

#load required libraries
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn
from sklearn.neighbors import NearestNeighbors #pip install scikit-learn
# from sklearn.neighbors import DistanceMetric
from sklearn.metrics import DistanceMetric
from sklearn.decomposition import PCA
from kneed import KneeLocator as kl #pip install kneed
# verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

#set up directory
large_root = r"/Users/yifanwang/workspace/single_cell_course/datasets/filtered_gene_bc_matrices/hg19"
results_file = large_root + "\pbmc3k.h5ad"

#load data
print("Reading data...")
adata = sc.read_10x_mtx(large_root, var_names = 'gene_symbols',
                        cache=True)

# annotate mitochondrial genes as 'mt' and calculate qc metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

QC_data = adata.copy()
#basic filtering
sc.pp.filter_cells(QC_data, min_genes=200)
sc.pp.filter_genes(QC_data, min_cells=3)

# actually do the filtering by slicing the AnnData object
QC_data = QC_data[QC_data.obs.n_genes_by_counts < 2200, :]
QC_data = QC_data[QC_data.obs.pct_counts_mt < 5, :]

norm_data = QC_data.copy()
# normalize counts to 10,000 reads per cell..
sc.pp.normalize_total(norm_data, target_sum=1e4)

#log n+1 transform
sc.pp.log1p(norm_data)

# identify highly variable genes
sc.pp.highly_variable_genes(norm_data, min_mean=0.0125, max_mean=3, min_disp=0.5)

#plot dispersion vs expression
sc.pl.highly_variable_genes(norm_data)

HVG_data = norm_data.copy()
#filter to only include HVGs
HVG_data.raw = HVG_data # freeze state of HVG_data by saving as raw attribute
HVG_data = HVG_data[:, HVG_data.var.highly_variable]

pcHVG_data = HVG_data.copy()
#scale to unit variance for PCA, clip values exceding var 10
sc.pp.scale(pcHVG_data, max_value=10)

# reduce the dimensionality of the data using PCA
sc.tl.pca(pcHVG_data, svd_solver='arpack')

#save the result
pcHVG_data.write(results_file)
#print
pcHVG_data

#scatter in the PCA coordinates of HVGs, won't be used later on, only for vis
sc.pl.pca(pcHVG_data, color='CST3')

#plot the loadings, or the contribution of each gene to the PCs
sc.pl.pca_loadings(pcHVG_data)

#scatter in the PCA coordinates of all genes
sc.pp.scale(norm_data, max_value=10)
sc.tl.pca(norm_data, svd_solver='arpack')
sc.pl.pca(norm_data, color='CST3')
sc.pl.pca_loadings(norm_data)

#image you have two samples, and you want to check how similar the PCs are
#how can we quantify the similarity of gene contributions to PCs?

#how many PCs should we use to cluster?
# visualize scree plot
sc.pl.pca_variance_ratio(pcHVG_data, log=True)

#define new function that finds the elbow dimension of the scree plot
def PCA_Elbow_fit(data):
    model = PCA().fit(data)
    explained_variance = model.explained_variance_ratio_
    pcs = list(range(1,explained_variance.shape[0]+1))#enumerate(explained_variance,1)
    klm = kl(pcs, explained_variance, S=1.0, curve='convex', direction='decreasing')
    pcs_used = klm.knee
    pc_list = list(range(1,pcs_used+1))
    new_data = PCA(n_components= pcs_used,svd_solver= "arpack").fit_transform(data)
    return pcs_used, new_data, pcs, explained_variance,pc_list

#extract cell x gene data from scanpy annData and create new pandas dataframe
label = "HVG_data"
new_frame = pd.DataFrame(HVG_data.X.toarray(), index=HVG_data.obs_names, columns=HVG_data.var_names)
pandas_data=new_frame.values
#new_frame.to_csv(output_path)

#execute PCA elbow fitting function
dim,new_matrix,pc_ax,pc_ay,col_labels=PCA_Elbow_fit(pandas_data)
print(dim)

#create new data frame of cells x PC values
columns=["PC_"+str(i) for i in col_labels]
output_path = large_root+ "_PCA_"+ label + str(dim)+ ".csv"
PC_frame= pd.DataFrame(new_matrix, index= new_frame.index.values.tolist(), columns = columns)
print(PC_frame.head())
print(PC_frame.shape)
#PC_frame.to_csv(output_path)

fig= plt.figure(figsize = (6,6))
ax1 = fig.add_subplot(111)
ax1.set_title("Skree Plot for "+label+" (Elbow at "+str(dim)+")")
ax1.plot(pc_ax, pc_ay)
ax1.set_xlabel("PC Number")
ax1.set_ylabel("Explained Variance Ratio")
fig.tight_layout()
#plt.savefig(large_root+"/_PCA_"+str(dim)+".png")
#fig.clear()

#back to scanpy
#computing the neighborhood graph
sc.pp.neighbors(pcHVG_data, n_neighbors=10, n_pcs=dim)

#visualizing data with umap
sc.tl.umap(pcHVG_data)

#visualize marker gene expression
sc.pl.umap(pcHVG_data, color=['CST3', 'NKG7', 'PPBP'])

#visualize data with t-SNE
sc.tl.tsne(pcHVG_data, n_pcs = dim)
sc.pl.tsne(pcHVG_data, color=['CST3', 'NKG7', 'PPBP'])

# read raw cell by count matrix in pandas
# what is the fraction of non-zero entries in the matrix?



# we created the following annData objects:
# adata - raw cell by gene count matrix
# QC_data - QC thresholds applied to remove low quality cells and sparsely expressed genes
# norm_data - CPM normalized and log+1 transformed data
# HVG_data - highly variable genes only, cell by CPM & log+1 transformed gene counts
# pcHVG_data - HVG_data with dimensions reduced using PCA


# how many genes were removed by the QC filters? how many cells were removed?



# how many different mitochondrial genes were observed?



# is the PC elbow dimension the same for HVGs and raw data? what about differently filtered or transformed data?



# how will the umap change if we use non-reduced data?

#computing the neighborhood graph on raw data, no PCs
sc.pp.neighbors(pcHVG_data, n_neighbors=10, n_pcs=0)
sc.tl.umap(pcHVG_data)
sc.pl.umap(pcHVG_data, color=['CST3', 'NKG7', 'PPBP'])

#how will the umap change if we use different numbers of PCs?

#computing the neighborhood graph, 2 PCs
sc.pp.neighbors(pcHVG_data, n_neighbors=10, n_pcs=2)
sc.tl.umap(pcHVG_data)
sc.pl.umap(pcHVG_data, color=['CST3', 'NKG7', 'PPBP'])

#computing the neighborhood graph, 10 PCs
sc.pp.neighbors(pcHVG_data, n_neighbors=10, n_pcs=10)
sc.tl.umap(pcHVG_data)
sc.pl.umap(pcHVG_data, color=['CST3', 'NKG7', 'PPBP'])

#computing the neighborhood graph, 50 PCs
sc.pp.neighbors(pcHVG_data, n_neighbors=10, n_pcs=50)
sc.tl.umap(pcHVG_data)
sc.pl.umap(pcHVG_data, color=['CST3', 'NKG7', 'PPBP'])

#what if we change the number of neighbors?

#2 neighbors
sc.pp.neighbors(pcHVG_data, n_neighbors=2, n_pcs=20)
sc.tl.umap(pcHVG_data)
sc.pl.umap(pcHVG_data, color=['CST3', 'NKG7', 'PPBP'])

#30 neighbors
sc.pp.neighbors(pcHVG_data, n_neighbors=30, n_pcs=20)
sc.tl.umap(pcHVG_data)
sc.pl.umap(pcHVG_data, color=['CST3', 'NKG7', 'PPBP'])

#100 neighbors
sc.pp.neighbors(pcHVG_data, n_neighbors=100, n_pcs=20)
sc.tl.umap(pcHVG_data)
sc.pl.umap(pcHVG_data, color=['CST3', 'NKG7', 'PPBP'])
































#












