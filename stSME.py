import stlearn as st
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc
st.settings.set_figure_params(dpi=180)
# most of the code is from stSME tutorial

# specify PATH to data
BASE_PATH = Path("/mnt/CMML_2/DLPFC/")

# spot tile is the intermediate result of image pre-processing
TILE_PATH = Path("/mnt/CMML_2/tiles")
TILE_PATH.mkdir(parents=True, exist_ok=True)

# output path
OUT_PATH = Path("/mnt/CMML_2")
OUT_PATH.mkdir(parents=True, exist_ok=True)

#input data
data = st.Read10X(BASE_PATH)

# pre-processing for gene count table
st.pp.filter_genes(data,min_cells=1)
st.pp.normalize_total(data)
st.pp.log1p(data)

# pre-processing for spot image
st.pp.tiling(data, TILE_PATH)

# this step uses deep learning model to extract high-level features from tile images
# may need few minutes to be completed
st.pp.extract_feature(data)

# run PCA for gene expression data
st.em.run_pca(data,n_comps=50)

data_SME = data.copy()
# apply stSME to normalise log transformed data
st.spatial.SME.SME_normalize(data_SME, use_data="raw")
data_SME.X = data_SME.obsm['raw_SME_normalized']
st.pp.scale(data_SME)
st.em.run_pca(data_SME,n_comps=50)

# louvain clustering on stSME normalised data
st.pp.neighbors(data_SME,n_neighbors=15,use_rep='X_pca')
st.tl.clustering.louvain(data_SME, resolution=1.0)

# data_SME.write("data_sme_dlpfc.h5ad")
# data_SME = sc.read_h5ad("/mnt/CMML_2/data_sme_dlpfc.h5ad")


# 7 color-Viridis palette, reordered to match the expected brain layers
palette = [
    "#FFD700",  # 7
    "#443983",  # 1
    "#35b779",  # 4
    "#31688e",  # 2
    "#21918c",  # 3
    "#90d743",  # 5
    '#f5e801',  # 6
    "#440154",  # 0
]

# mapping color into category levels
data_SME.uns['louvain_colors'] = palette 

# visualize clustering results
sc.pl.spatial(
    data_SME,
    color    = 'louvain',
    img_key  = None,     # no image, make plot more clear
    size     = 1.5,
    alpha    = 0.8,
    cmap     = 'viridis',
    show     = True
)

plt.savefig(
    "/mnt/CMML_2/5.png",
    dpi = 300,
    bbox_inches = 'tight'
)

#save the clustering results, prepared for benchmarking in R
data_SME.obs['louvain'].to_csv("/mnt/CMML_2/stSME_labels.csv", header=True)
