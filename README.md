# Improved scanpy package
Scanpy is a scalable toolkit for analyzing single-cell gene expression data built jointly with anndata. It includes preprocessing, visualization, clustering, trajectory inference and differential expression testing. 

This is a modified version of the scanpy package from scverse, incorporating enhancements and changes, including dynamic UMAP plotting using the Holoviews package.


 ## Installation
 To install this modified version of the scanpy package, follow these steps:

1. Clone this repository to your local machine:
 ```bash
git clone https://github.com/cmar24/OOP-project.git 
 ```
2. Navigate to the cloned directory:
 ```bash
 cd OOP-project
 ```
3. Install the package using pip:
```bash
pip install .
```
# Usage
After installation, you can use the modified scanpy package in your Python environment. Here's a quick overview of the new implemented function. For further detail consult the tutorial.

The scanpy tutorial with its original functions can be found on the following link: https://scanpy-tutorials.readthedocs.io/en/latest/index.html 

1. `umap_din`: Generates a UMAP scatter plot based on the specified visualization option. It accepts an AnnData object containing the data and a visualization option parameter, allowing the user to choose from different options such as plotting based on cluster labels, total counts of each cell, labels for bulk cell types, or marker genes. The function returns the generated scatter plot using hvplot. Additionally, for the "clusters" option, it includes the top 5 expressed genes as hover information for each point on the plot.

2. `tsne_din`: Similar to `umap_din`, this function generates a t-SNE scatter plot based on the specified visualization option. It accepts an AnnData object containing the data and a visualization option parameter. The user can choose from different options such as plotting based on cluster labels, total counts of each cell, labels for bulk cell types, or marker genes. The function returns the generated scatter plot using hvplot. Additionally, for the "clusters" option, it includes the top 5 expressed genes as hover information for each point on the plot.

3. `dimensionality_reduction`: Performs data reduction based on the input data attributes, allowing the user to filter data based on either highly variable genes or cell types. It recalculates the PCA, UMAP, and tSNE coordinates based on the filtered data.

4. `export_data_to_csv`: Exports data from the complex AnnData data structure to a CSV file based on the specified information, including raw or normalized expression data, information about variable genes, or coordinates from PCA, UMAP, or tSNE plots.

5. `export_scatter_data_to_csv`: Exports data from a scatter plot to a CSV file, including coordinates and top expressed genes for each point on the plot. The type of plot (UMAP, tSNE, or PCA) and the number of top expressed genes to include can be specified.

These functions provide useful utilities for analyzing and exporting single-cell RNA-seq data, allowing users to perform dimensionality reduction, filter data based on specific attributes, and export data for further analysis.

