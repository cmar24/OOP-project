import numpy as np
import pandas as pd
import hvplot.pandas
import ipywidgets as widgets
from IPython.display import display, clear_output

def top_expressed_genes(adata, x_coord, y_coord, graph_type, n_genes=5):
    """
    Function to retrieve top expressed genes in a cell based on its UMAP/PCA/t-SNE coordinates.

    Parameters:
    - adata (anndata.AnnData): Annotated data object containing the PBMC data.
    - x_coord (float): X-coordinate of the map.
    - y_coord (float): Y-coordinate of the map.
    - graph_type (str): Type of graph to which the coordinate correspond ('umap', 'pca', or 'tsne').
    - n_genes (int): Number of top expressed genes to retrieve.

    Returns:
    - top_genes (list): List of top expressed genes in the specified cell.
    """

    if graph_type == "umap":
        nearest_cell_index = ((adata.obsm['X_umap'] - np.array([x_coord, y_coord]))**2).sum(axis=1).argmin()
    elif graph_type == "pca":
        nearest_cell_index = ((adata.obsm['X_pca'] - np.array([x_coord, y_coord]))**2).sum(axis=1).argmin()
    elif graph_type == "tsne":
        nearest_cell_index = ((adata.obsm['X_tsne'] - np.array([x_coord, y_coord]))**2).sum(axis=1).argmin()
    else:
        raise ValueError("Invalid graph type. Choose 'umap', 'pca', or 'tsne'.")

    cell_expression = adata.X[nearest_cell_index].toarray().flatten()

    sorted_genes_indices = cell_expression.argsort()[::-1][:n_genes]

    top_genes = [adata.var_names[i] for i in sorted_genes_indices]

    return top_genes


def umap_din(pbmc, visualization_option):
    """
    Generate a UMAP scatter plot based on the specified visualization option.

    Parameters:
    - pbmc (anndata.AnnData): The AnnData object containing the data.
    - visualization_option (str): The visualization option to choose from:
        - "clusters": Plot based on cluster labels.
        - "n_counts": Plot based on the total counts of each cell.
        - "bulk_labels": Plot based on labels for bulk cell types.
        - "explore_genes": Plot based on marker genes. Includes a menu to select genes more easily.

    Returns:
    - scatter_plot (hvplot): The generated scatter plot using hvplot.
    """

    # Extract coordinates data from the AnnData object
    umap_x = pbmc.obsm["X_umap"][:, 0]
    umap_y = pbmc.obsm["X_umap"][:, 1]
    
    if visualization_option == "clusters":
        if "louvain" in pbmc.obs.columns:
            labels = pbmc.obs["louvain"]
        elif "louvain_cell_types" in pbmc.obs.columns:
            labels = pbmc.obs["louvain_cell_types"]
        else:
            raise ValueError("No cluster labels found in AnnData object.")
            
        # Create a DataFrame for hvplot
        data = pd.DataFrame({"UMAP_X": umap_x, "UMAP_Y": umap_y, "Labels": labels})
        top_genes = data.apply(lambda row: ', '.join(top_expressed_genes(pbmc,row["UMAP_X"], row["UMAP_Y"], "umap")), axis=1)
        data["Top_5_Genes"] = top_genes

        scatter_plot = data.hvplot.scatter(x="UMAP_X", y="UMAP_Y", by="Labels", hover_cols=["Top_5_Genes"])
            
    elif visualization_option == "n_counts":
        labels = pbmc.obs["n_counts"]
        data = pd.DataFrame({"UMAP_X": umap_x, "UMAP_Y": umap_y, "Counts": labels})
        # we add a gradient color bar
        scatter_plot = data.hvplot.scatter(x="UMAP_X", y="UMAP_Y", c="Counts", cmap="viridis", colorbar=True, hover_cols=["Counts"])
            
            
    elif visualization_option == "bulk_labels":
        labels = pbmc.obs["bulk_labels"]
        data = pd.DataFrame({"UMAP_X": umap_x, "UMAP_Y": umap_y, "Labels": labels})
        # Create and display the scatter plot using hvplot.scatter
        scatter_plot = data.hvplot.scatter(x="UMAP_X", y="UMAP_Y", by="Labels")
        
    elif visualization_option == "explore_genes":
        # Define a function to update scatter plot based on selected gene
        def update_umap(change):
            selected_gene = change.new
            sc.pl.umap(pbmc, color=selected_gene)

        # Create dropdown widget for gene selection
        gene_dropdown = widgets.Dropdown(
            options=pbmc.var_names.tolist(),
            description='Select Gene:',
            disabled=False,
        )

        # Attach the update function to the dropdown's event
        gene_dropdown.observe(update_umap, names='value')

        # Display dropdown widget
        display(gene_dropdown)
        # Return None as scatter_plot since it will be dynamically updated
        scatter_plot = None
        
    else:
        raise ValueError("Invalid visualization option. Please choose from 'clusters', 'n_counts', 'bulk_labels', or 'explore_genes'.")

    return scatter_plot

def tsne_din(pbmc, visualization_option):
    """
    Generate a t-SNE scatter plot based on the specified visualization option.

    Parameters:
    - pbmc (anndata.AnnData): The AnnData object containing the data.
    - visualization_option (str): The visualization option to choose from:
        - "clusters": Plot based on cluster labels.
        - "n_counts": Plot based on the total counts of each cell.
        - "bulk_labels": Plot based on labels for bulk cell types.
        - "explore_genes": Plot based on marker genes. Includes a menu to select genes more easily.

    Returns:
    - scatter_plot (hvplot): The generated scatter plot using hvplot.
    """

    # Extract coordinates data from the AnnData object
    tsne_x = pbmc.obsm["X_tsne"][:, 0]
    tsne_y = pbmc.obsm["X_tsne"][:, 1]
    
    if visualization_option == "clusters":
        if "louvain" in pbmc.obs.columns:
            labels = pbmc.obs["louvain"]
        elif "louvain_cell_types" in pbmc.obs.columns:
            labels = pbmc.obs["louvain_cell_types"]
        else:
            raise ValueError("No cluster labels found in AnnData object.")
            
        # Create a DataFrame for hvplot
        data = pd.DataFrame({"tSNE_X": tsne_x, "tSNE_Y": tsne_y, "Labels": labels})
        top_genes = data.apply(lambda row: ', '.join(top_expressed_genes(pbmc,row["tSNE_X"], row["tSNE_Y"], "tsne")), axis=1)
        data["Top_5_Genes"] = top_genes

        scatter_plot = data.hvplot.scatter(x="tSNE_X", y="tSNE_Y", by="Labels", hover_cols=["Top_5_Genes"])
            
    elif visualization_option == "n_counts":
        labels = pbmc.obs["n_counts"]
        data = pd.DataFrame({"tSNE_X": tsne_x, "tSNE_Y": tsne_y, "Counts": labels})
        # we add a gradient color bar
        scatter_plot = data.hvplot.scatter(x="tSNE_X", y="tSNE_Y", c="Counts", cmap="viridis", colorbar=True, hover_cols=["Counts"])
            
            
    elif visualization_option == "bulk_labels":
        labels = pbmc.obs["bulk_labels"]
        data = pd.DataFrame({"tSNE_X": tsne_x, "tSNE_Y": tsne_y, "Labels": labels})
        # Create and display the scatter plot using hvplot.scatter
        scatter_plot = data.hvplot.scatter(x="tSNE_X", y="tSNE_Y", by="Labels")
        
    elif visualization_option == "explore_genes":
        # Define a function to update scatter plot based on selected gene
        def update_tsne(change):
            selected_gene = change.new
            sc.pl.tsne(pbmc, color=selected_gene)

        # Create dropdown widget for gene selection
        gene_dropdown = widgets.Dropdown(
            options=pbmc.var_names.tolist(),
            description='Select Gene:',
            disabled=False,
        )

        # Attach the update function to the dropdown's event
        gene_dropdown.observe(update_tsne, names='value')

        # Display dropdown widget
        display(gene_dropdown)
        # Return None as scatter_plot since it will be dynamically updated
        scatter_plot = None
        
    else:
        raise ValueError("Invalid visualization option. Please choose from 'clusters', 'n_counts', 'bulk_labels', or 'explore_genes'.")

    return scatter_plot