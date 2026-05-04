import pandas as pd

rna = pd.read_csv("results/regulons/tf_activity_rna_by_celltype.csv")
atac = pd.read_csv("results/regulons/atac_motif_by_celltype.csv")

print("RNA:", sorted(rna["cell_type"].astype(str).str.strip().unique()))
print("ATAC:", sorted(atac["cell_type"].astype(str).str.strip().unique()))