from scanpy.plotting.plot import umap_din

from scanpy import *

pbmc = datasets.pbmc68k_reduced()

print(pbmc)

umap_din(pbmc)








