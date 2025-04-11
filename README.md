# A Step-by-Step Guide to Remapping Unannotated Follicular Cells from FlyCellAtlas
This guide explains the methodology and provides annotated code for reclassifying unannotated ovarian follicular cells in Drosophila melanogaster using single-cell RNA-seq data obtained from FlyCellAtlas (Li et al., 2022; Jia et al., 2022). The objective is to identify novel subpopulations corresponding to distinct stages of follicular development by reprocessing and reclustering the unannotated cells.

1. Prerequisites and Setup
Before you begin, ensure you have the following Python packages installed:

Scanpy for single-cell data analysis

Pandas and NumPy for data manipulation

Anndata for handling single-cell data objects

You can install these packages via pip (if not already installed):


pip install scanpy pandas numpy anndata
2. Data Loading and Initial Filtering
First, the script loads the single-cell dataset from an .h5ad file. The dataset contains various cells from Drosophila ovaries. In this example, we isolate the subset of cells with an annotation labeled as “unannotated”:


import scanpy as sc
import pandas as pd
import numpy as np

# Load the original anndata object (e.g., .h5ad file)
adata = sc.read_h5ad("new_fly_single/s_fca_biohub_ovary_10x.h5ad")

# Filter cells that are marked as "unannotated"
unannotated_cells = adata[adata.obs['annotation'] == 'unannotated'].copy()
print(unannotated_cells)
Explanation:
This step isolates cells without a predefined annotation for further analysis.

3. Dimensionality Reduction and UMAP Construction
The next steps involve reducing the dimensionality of the data using PCA and constructing a UMAP embedding. This helps visualize the global structure of the unannotated cells:


# Perform PCA on the unannotated cells
sc.tl.pca(unannotated_cells, svd_solver='arpack')

# Compute the neighborhood graph using PCA results
sc.pp.neighbors(unannotated_cells, n_neighbors=15, n_pcs=20)

# Compute the UMAP embedding for visualization
sc.tl.umap(unannotated_cells)
Explanation:

PCA (Principal Component Analysis): Reduces data complexity by projecting it onto a limited set of principal components.

Neighbors computation: Uses the PCA-reduced space to construct a graph representation of cell-to-cell similarity.

UMAP: A non-linear dimensionality reduction technique that preserves both local and global data structure in a low-dimensional space.

4. Clustering Using the Leiden Algorithm
Clustering is performed on the unannotated cells to delineate distinct subpopulations. The Leiden algorithm is used here with an adjustable resolution parameter:


# Perform clustering using the Leiden algorithm (adjust the resolution as needed)
sc.tl.leiden(unannotated_cells, resolution=1.0)
sc.pl.umap(unannotated_cells, color=['leiden'], title="UMAP Unannotated Subset (by Leiden)")
Explanation:
A resolution of 1.0 is initially chosen, though later in the pipeline a resolution of 0.3 is also applied for a more conservative split. Clustering subdivides the continuous cellular space into discrete groups based on gene expression similarities.

5. Identification of Marker Genes per Cluster
Once clusters are defined, the next step is to identify key marker genes responsible for the observed differences. Scanpy’s built-in function extracts the most significant differentially expressed genes for each cluster:

# Identify marker genes for each cluster using the Wilcoxon method
sc.tl.rank_genes_groups(
    unannotated_cells, 
    groupby='leiden', 
    method='wilcoxon', 
    use_raw=False  # Adjust this if you use an alternative data layer
)

# Visualize the top 10 marker genes per cluster
sc.pl.rank_genes_groups(unannotated_cells, n_genes=10, sharey=False)

# Extract the results into a DataFrame and save to CSV for further inspection
df_markers = sc.get.rank_genes_groups_df(unannotated_cells, group=None)
df_markers.to_csv("markers_all_clusters-leiden.csv", index=False)
Explanation:
This step produces a ranked list of marker genes per cluster that you can compare against prior literature (e.g., Tootle et al., 2011; Zartman et al., 2009) to link clusters with developmental stages.

6. Remapping and Merging of Annotations
A second part of the analysis involves reassigning new annotations to the unannotated cells based on the Leiden clusters. The code below demonstrates how to map cluster IDs to biologically interpretable stage labels (e.g., Stage 9, Stage 10A, Stage 10B/11, adipocytes):

# Using a pre-defined mapping dictionary for unannotated cells:
mapping = {
    "0": "Stage 9",
    "2": "Stage 10A",
    "6": "Stage 10B/11",
    "4": "adipocytes"
}

# Generate a sub-clustering annotation in the unannotated cells
sc.tl.leiden(unannotated_cells, resolution=0.3)
unannotated_cells.obs['sub_leiden'] = unannotated_cells.obs['leiden']

# Propagate the sub_clustering results to the full dataset (adata) for visualization
adata.obs.loc[unannotated_cells.obs_names, 'sub_leiden'] = unannotated_cells.obs['sub_leiden']
sc.pl.umap(adata, color='sub_leiden', title='Global UMAP colored by sub_leiden')

# Create a new annotation field and map the new labels
adata.obs['annotation_merged'] = adata.obs['annotation'].astype(str)
mask_unann = adata.obs['annotation'] == 'unannotated'
adata.obs['sub_leiden'] = adata.obs['sub_leiden'].astype(str)

adata.obs.loc[mask_unann, 'annotation_merged'] = adata.obs.loc[mask_unann, 'sub_leiden'].apply(
    lambda x: mapping.get(x, "unannotated")
)

# Inspect the final annotation distribution
print("Final distribution of annotation_merged:")
print(adata.obs['annotation_merged'].value_counts())

# Visualize the global UMAP with the updated annotations
sc.pl.umap(adata, color='annotation_merged', legend_loc="right margin", title="UMAP with updated annotations")
Explanation:

A dictionary (mapping) is used to convert cluster labels into biologically relevant stages.

The new annotations are merged back into the original dataset for a global view.

Final UMAP plots help verify that the remapping reflects the expression patterns observed.

7. Checking Marker Gene Expression
To validate the new annotations, expression levels of known marker genes are plotted. This step cross-checks whether the predicted stages correspond with expected transcriptional profiles:

# Generate a dot plot of marker gene expression across the new annotated groups
sc.pl.dotplot(
    adata,
    var_names=['Cad74A', 'Cad87A', 'Cad88C', 'Cad99C', 'bond', 'CG5326', 'PH4alphaEFB', 'ndl', 'wbl', 'Ilp6'],
    groupby="annotation_merged",
    dendrogram=True,
    title="Expression of Marker Genes by group"
)
Explanation:

The dot plot visualizes the expression of key genes (e.g., the Cadherin family and others) in relation to the newly assigned follicular stages.

Notice that the annotation “adipocytes” was assigned solely based on high Ilp6 expression, following the previous annotation by Wu-Min Deng.

8. Finalizing the Analysis
Lastly, the updated anndata object is saved to an .h5ad file for downstream analyses or sharing with collaborators:

python
Copiar
# Save the updated anndata object
adata.write("ovary_flycellatlas.h5ad")
