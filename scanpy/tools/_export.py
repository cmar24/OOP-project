import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scanpy.plotting._dinamicplot import top_expressed_genes

def export_data_to_csv(filename, pbmc, info):
    """
    Export data from the complex AnnData data structure to a CSV file.

    Parameters:
        filename (str): The name of the CSV file to be exported.
        pbmc (AnnData): An AnnData object containing the single-cell data.
        info (str): Type of information to export. Choose from:
            - 'raw': Export raw expression data.
            - 'normalized': Export normalized expression data.
            - 'variable_genes': Export information about variable genes and related information.
            - 'umap': Export UMAP coordinates.
            - 'tsne': Export tSNE coordinates.
            - 'pca': Export PCA coordinates.

    Raises:
        ValueError: If the 'info' parameter is not one of the valid options.

    Returns:
        None
    """
    if info == "raw":
        data = pd.DataFrame(pbmc.layers['raw'].todense(), columns=pbmc.var_names)
    elif info == "normalized":
        data = pd.DataFrame(pbmc.X.todense(), columns=pbmc.var_names)
    elif info == "variable genes":
        data = pbmc.var[['gene_names', 'highly_variable', 'means', 'dispersions', 'dispersions_norm']]
    elif info == "umap":
        ids = pbmc.obs_names
        x_coord = pbmc.obsm["X_umap"][:, 0]
        y_coord = pbmc.obsm["X_umap"][:, 1]
        data = pd.DataFrame({
        "ID": ids,
        "UMAP_X": x_coord,
        "UMAP_Y": y_coord
    })
    elif info == "tsne":
        ids = pbmc.obs_names
        x_coord = pbmc.obsm["X_tsne"][:, 0]
        y_coord = pbmc.obsm["X_tsne"][:, 1]
        data = pd.DataFrame({
        "ID": ids,
        "tSNE_X": x_coord,
        "tSNE_Y": y_coord})
    elif info == "pca":
        ids = pbmc.obs_names
        x_coord = pbmc.obsm["X_pca"][:, 0]
        y_coord = pbmc.obsm["X_pca"][:, 1]
        data = pd.DataFrame({
        "ID": ids,
        "PCA_X": x_coord,
        "PCA_Y": y_coord})
    else:
        raise ValueError("Invalid info type. Choose 'raw', 'normalized', 'variable genes', 'umap', 'pca', or 'tsne'.")

    data.to_csv(filename + ".csv", index=False)
    print("Exported UMAP data to CSV file:", filename + ".csv")

def export_scatter_data_to_csv(filename, pbmc, plot, num_top_genes=5):
    """
    Export data from a scatter plot to a CSV file.

    Parameters:
        filename (str): The name of the CSV file to be exported.
        pbmc (anndata.AnnData): Anndata object containing the single-cell data.
        plot (str): Type of plot to export data from. Choose from:
            - 'umap': Export data from a UMAP plot.
            - 'tsne': Export data from a tSNE plot.
            - 'pca': Export data from a PCA plot.
        num_top_genes (int): Number of top expressed genes to include in the exported data. Default is 5.

    Raises:
        ValueError: If the 'plot' parameter is not one of the valid options ('umap', 'tsne', 'pca').

    Returns:
        None
    """
    ids = pbmc.obs_names

    if plot == "umap":
        x_coord = pbmc.obsm["X_umap"][:, 0]
        y_coord = pbmc.obsm["X_umap"][:, 1]
    elif plot == "tsne":
        x_coord = pbmc.obsm["X_tsne"][:, 0]
        y_coord = pbmc.obsm["X_tsne"][:, 1]
    elif plot == "pca":
        x_coord = pbmc.obsm["X_pca"][:, 0]
        y_coord = pbmc.obsm["X_pca"][:, 1]
    else:
        raise ValueError("Invalid plot type. Choose 'umap', 'pca', r 'tsne'.")

    top_genes = []
    for x, y in zip(x_coord, y_coord):
        top_genes.append(top_expressed_genes(pbmc, x, y, plot, num_top_genes))

    umap_data = pd.DataFrame({
        "ID": ids,
        "UMAP_X": x_coord,
        "UMAP_Y": y_coord,
        "Top_Expressed_Genes": top_genes
    })

    umap_data.to_csv(filename + ".csv", index=False)
    print("Exported UMAP data to CSV file:", filename + ".csv")