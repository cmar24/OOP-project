# Improved scanpy package
Scanpy is a scalable toolkit for analyzing single-cell gene expression data built jointly with anndata. It includes preprocessing, visualization, clustering, trajectory inference and differential expression testing. 

This is a modified version of the scanpy package from scverse, incorporating enhancements and changes, including dynamic UMAP and tSNE plotting using hvplot and ipywidgets packages, as well as a tool for reducing the dimentionality of the data and functions to make the exportation of data easier for non-expert users. 

 ## Installation
 To install this modified version of the scanpy package, we recommend that you create a python virtual environment to install and use the library, following this steps:

1. Create a virtual environment and activate it: 
```bash
python -m venv name_of_venv

name_of_venv\Scripts\ativate #windows
source name_of_venv/bin/activate #linux
```

2. Clone this repository to your local machine:
 ```bash
git clone https://github.com/cmar24/OOP-project.git 
 ```
3. Navigate to the cloned directory:
 ```bash
 cd OOP-project
 ```
4. Install the package using pip:
```bash
pip install .
```
Now the new package has been installed.
# Usage
After installation, you can use the modified scanpy package in your Python environment you just created. 

The scanpy tutorial with its original functions can be found on the following link: https://scanpy-tutorials.readthedocs.io/en/latest/index.html. 

Here's a quick overview of the new implemented functions. For further detail consult the tutorial, available in the git repository.

1. `top_expressed_genes`: Returns the list of the most expressed genes in a cell when given the coordinates of that cell from a umap or tsne. UMAP scatter plot based on the specified visualization option. Here is how to call the function:
```python
sc.pl.top_expressed_genes()
```

2. `umap_din`: Generates a UMAP scatter plot based on the specified visualization option. It accepts an AnnData object containing the data and a visualization option parameter, allowing the user to choose from different options such as plotting based on cluster labels, total counts of each cell, labels for bulk cell types, or marker genes. The function returns the generated scatter plot using hvplot. Additionally, for the "clusters" option, it includes the top 5 expressed genes as hover information for each point on the plot. Here is how to call the function:
```python
sc.pl.umap_din()
```

3. `tsne_din`: Similar to `umap_din`, this function generates a t-SNE scatter plot based on the specified visualization option. It accepts an AnnData object containing the data and a visualization option parameter. The user can choose from different options such as plotting based on cluster labels, total counts of each cell, labels for bulk cell types, or marker genes. The function returns the generated scatter plot using hvplot. Additionally, for the "clusters" option, it includes the top 5 expressed genes as hover information for each point on the plot. Here is how to call the function:
```python
sc.pl.tsne_din()
```

4. `dimensionality_reduction`: Performs data reduction based on the input data attributes, allowing the user to filter data based on either highly variable genes or cell types. It recalculates the PCA, UMAP, and tSNE coordinates based on the filtered data. Here is how to call the function:
```python
sc.pp.dimentionality_reduction()
```

5. `export_data_to_csv`: Exports data from the complex AnnData data structure to a CSV file based on the specified information, including raw or normalized expression data, information about variable genes, or coordinates from PCA, UMAP, or tSNE plots. Here is how to call the function:
```python
sc.tl.export_data_to_csv()
```

6. `export_scatter_data_to_csv`: Exports data from a scatter plot to a CSV file, including coordinates and top expressed genes for each point on the plot. The type of plot (UMAP, tSNE, or PCA) and the number of top expressed genes to include can be specified. Here is how to call the function:
```python
sc.tl.export_scatter_data_to_csv()
```

These functions provide useful utilities for analyzing and exporting single-cell RNA-seq data, allowing users to perform dimensionality reduction, filter data based on specific attributes, and export data for further analysis.

