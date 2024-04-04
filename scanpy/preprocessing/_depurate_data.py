import matplotlib.pyplot as plt
import pandas as pd
from ._pca import pca
from ..neighbors import neighbors
from ..tools._umap import umap
from ..tools._tsne import tsne

def dimensionality_reduction(pbmc, input_data):
    """
    Perform dimensionality reduction based on input_data.

    Parameters:
        pbmc (anndata.AnnData): Anndata object containing the data.
        input_data (str): Either 'highly_variable' or 'cell_types'.

    Returns:
        anndata.AnnData: Anndata object with dimensionality reduction results.
    """
    pbmc_red = pbmc.copy()
    
    if input_data == 'highly_variable':
        # the plot helps the scientitst decide if it is useful to delete some data
        value_counts = pbmc.var["highly_variable"].value_counts()
        value_counts_df = pd.DataFrame(value_counts).reset_index()
        value_counts_df.columns = ['highly_variable', 'count']
        highly_variable_genes = pbmc.var.index[pbmc.var['highly_variable']].tolist()
        not_highly_variable_genes = pbmc.var.index[~pbmc.var['highly_variable']].tolist()
        value_counts_df['genes'] = ['-'.join(highly_variable_genes) if x else '-'.join(not_highly_variable_genes) for x in value_counts_df['highly_variable']]
        
        # Generating bar plot using Matplotlib
        plt.figure(figsize=(5,4))
        plt.bar(value_counts_df['highly_variable'], value_counts_df['count'])
        plt.ylabel('Count')
        plt.title('Highly Variable Genes')
        plt.xticks(value_counts_df['highly_variable'], ['Not Highly Variable', 'Highly Variable'])
        plt.show()
        
        # the scientist decides if he/she eliminates the data
        delete_not_hv_genes = input("Do you want to delete not highly variable genes? (yes/no): ")
        if delete_not_hv_genes.lower() == 'yes':
            pbmc_red = pbmc_red[:, pbmc_red.var['highly_variable']]
        elif delete_not_hv_genes.lower() != 'no':
            raise ValueError("Invalid input. Please enter 'yes' or 'no'.")
        
        # recalculate the diferent matrix with the data left
        pca(pbmc_red)
        neighbors(pbmc_red)
        umap(pbmc_red)
        tsne(pbmc_red)
        
    elif input_data == 'cell_types':
        # Display bar plot of cell type counts
        cell_type_counts = pbmc.obs['louvain_cell_types'].value_counts()
        cell_type_counts.plot(kind='bar', color='skyblue')
        plt.ylabel('Count')
        plt.title('Cell Type Counts')
        plt.xticks(rotation=90)
        plt.show()
        
        # Display list of cell types
        cell_types = pbmc.obs['louvain_cell_types'].unique()
        print("Cell Types:")
        for idx, cell_type in enumerate(cell_types, start=1):
            print(f"{idx}. {cell_type}")
        
        delete_cell_type_idx = input(f"Which cell type do you want to delete? (Enter the corresponding number) [1-{len(cell_types)}]: ")
        try:
            delete_cell_type_idx = int(delete_cell_type_idx)
            if delete_cell_type_idx < 1 or delete_cell_type_idx > len(cell_types):
                raise ValueError("Invalid cell type index.")
        except ValueError:
            raise ValueError("Invalid input. Please enter a valid number.")
        
        delete_cell_type = cell_types[delete_cell_type_idx - 1]
        print(f"Deleting cell type: {delete_cell_type}")
        
        pbmc_red = pbmc_red[pbmc_red.obs['louvain_cell_types'] != delete_cell_type]
        
        pca(pbmc_red)
        neighbors(pbmc_red)
        umap(pbmc_red)
        tsne(pbmc_red)
    else:
        raise ValueError("Invalid input_data. Please provide either 'highly_variable' or 'cell_types'.")
    
    return pbmc_red