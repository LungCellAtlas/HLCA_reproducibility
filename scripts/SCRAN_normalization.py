import scanpy as sc
import sys
sys.path.append("./")
import scib_excerpts

adata = sc.read("../data/HLCA_core_h5ads/HLCA_v1_intermediates/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_RAW_filt_ann.h5ad")
# make smaller for testing:
# adata = adata[:2000,:10000].copy()
# normalize 
adata_norm = scib_excerpts.SCRAN_normalize(adata, log_transform=False)
# store
adata_norm.write("../data/HLCA_core_h5ads/HLCA_v1_intermediates/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_SCRAN_normalized.h5ad")
