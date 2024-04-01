import numpy as np

def find_top_expressed_genes(x_cord, y_cord, pbmc, num_genes=5):
    """
    Find the top expressed genes closest to the given coordinates in a single-cell RNA-seq dataset.

    Parameters:
    - x_cord (float): The x-coordinate.
    - y_cord (float): The y-coordinate.
    - pbmc (anndata.AnnData): An AnnData object containing single-cell RNA-seq data.
    - num_genes (int, optional): The number of top expressed genes to return. Default is 5.

    Returns:
    - top_genes (list): A list of gene names corresponding to the top expressed genes, with the most expressed genes first.
    """
    closest_distance = float('inf')
    closest_index = None
    for idx, (x, y) in enumerate(pbmc.obsm["X_umap"]):
        distance = np.sqrt((x - x_cord)**2 + (y - y_cord)**2)
        if distance < closest_distance:
            closest_distance = distance
            closest_index = idx

    list_counts = pbmc.X[closest_index]

    sorted_indices = np.argsort(list_counts)
    top_indices = sorted_indices[-num_genes:][::-1]  # Reverse the order to get most expressed genes first

    top_genes = []
    for i in top_indices:
        gene_name = pbmc.var.iloc[i].name
        top_genes.append(gene_name)

    return top_genes


