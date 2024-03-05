import pandas as pd
import hvplot.pandas

def my_plot_function():
    print("Hello world")

def umap_din(pbmc):
    # Extraer datos relevantes del objeto AnnData
    umap_x = pbmc.obsm["X_umap"][:, 0]
    umap_y = pbmc.obsm["X_umap"][:, 1]
    louvain_labels = pbmc.obs["louvain"]

    # Crear un DataFrame para hvplot
    data = pd.DataFrame({"UMAP_X": umap_x, "UMAP_Y": umap_y, "Louvain": louvain_labels})

    # Crear y mostrar el scatter plot utilizando hvplot.scatter
    scatter_plot = data.hvplot.scatter(x="UMAP_X", y="UMAP_Y", by="Louvain")
    # hvplot.show(scatter_plot)
    return scatter_plot
