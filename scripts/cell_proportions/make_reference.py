#!/usr/bin/env python
"""
make_reference_v6.py
Build the MuSiC reference for run_music2_4.R from the FULL HLCA parenchyma file,
with a custom cell-type column `celltype_v6`.

Why: the old reference (ad_HLCA_lung_parenchyma_ref_MUSIC.h5ad) dropped cells
whose ann_level_3 was "None", which removed plain smooth muscle (1,939 cells).

Gene handling is the same as data_prep_RNA_deconv_MUSIC-1.ipynb:
  raw counts from layers["counts"]; var index = ensembl_id without version;
  keep genes present in the PRESTO bulk raw-count file; duplicated IDs: keep last.

celltype_v6 (first matching rule wins):
  Alveolar macrophages      ann_level_4 == "Alveolar macrophages"
  Interstitial macrophages  ann_level_4 == "Interstitial macrophages"
  (dropped)                 ann_level_3 == "Macrophages" with other level 4 (Unknown)
  Lymphoid                  ann_level_2 == "Lymphoid" and ann_level_3 in T cell lineage /
                            Innate lymphoid cell NK / B cell lineage (B incl. plasma cells)
  Smooth muscle             ann_level_2 == "Smooth muscle" and ann_level_3 in
                            "None" / "Smooth muscle FAM83D+"
                            ("SM activated stress response" keeps its own label)
  Mesothelium               ann_level_2 == "Mesothelium"
  <ann_level_3>             everything else with a real level-3 label
                            (Fibroblasts, EC arterial, EC venous, Basal, AT1, AT2, ...)
  (dropped)                 ann_level_3 "None"/"Unknown" not covered above
The donor/cell-count filter in run_music2_4.R then decides which labels are used.

Run (in the env with anndata/h5py/scipy, e.g. envRpy):
  python make_reference_v6.py
Needs RAM for the selected counts matrix (~5-6 GB for ~350k cells x ~38k genes).
"""
import os
import numpy as np
import pandas as pd
import h5py
import scipy.sparse as sp
import anndata as ad

BASE = "/home/martin/Desktop/MyPHD/Projects/PRESTO/DNA_Methylation/Analysis_Martin/cell_proportions"
FULL = os.path.join(BASE, "scratch/Music/data/scRNAseq_reference/ad_HLCA_lung_parenchyma.h5ad")
BULK = os.path.join(BASE, "music2/data/RNAseq_PRESTO/PRESTO - mRNA Gene Expression (RNA Seq)/20200427_103972-001_rawcounts.txt")
OLD_REF = os.path.join(BASE, "music2/data/scRNAseq_reference/ad_HLCA_lung_parenchyma_ref_MUSIC.h5ad")   # only for a gene-set check
OUT = os.path.join(BASE, "music2/data/scRNAseq_reference/ad_HLCA_lung_parenchyma_ref_MUSIC_v6.h5ad")
CHUNK = 20000


def dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a], dtype=object)


def read_col(group, name):
    g = group[name]
    if isinstance(g, h5py.Group) and "values" in g:    # nullable string array (newer anndata)
        out = dec(g["values"][:])
        if "mask" in g:
            out[g["mask"][:].astype(bool)] = None
        return out
    if isinstance(g, h5py.Group):                      # categorical
        cats = dec(g["categories"][:])
        codes = g["codes"][:]
        out = np.empty(len(codes), dtype=object)
        out[codes >= 0] = cats[codes[codes >= 0]]
        out[codes < 0] = None
        return out
    return dec(g[:])


def index_name(group):
    return group.attrs.get("_index", "_index")


def label_v6(l2, l3, l4):
    lab = np.array([x if x not in ("None", "Unknown", None) else None for x in l3], dtype=object)
    lab[l3 == "Macrophages"] = None
    lab[l4 == "Alveolar macrophages"] = "Alveolar macrophages"
    lab[l4 == "Interstitial macrophages"] = "Interstitial macrophages"
    lymph = (l2 == "Lymphoid") & np.isin(l3, ["T cell lineage", "Innate lymphoid cell NK", "B cell lineage"])
    lab[lymph] = "Lymphoid"
    lab[(l2 == "Smooth muscle") & np.isin(l3, ["None", "Smooth muscle FAM83D+"])] = "Smooth muscle"
    lab[l2 == "Mesothelium"] = "Mesothelium"
    return lab


def main():
    print("Reading bulk gene IDs:", BULK)
    bulk_genes = pd.read_csv(BULK, sep="\t", usecols=[0]).iloc[:, 0].astype(str)
    bulk_genes = set(bulk_genes.str.replace(r"\.\d+$", "", regex=True))

    f = h5py.File(FULL, "r")
    obs, var = f["obs"], f["var"]
    cells = read_col(obs, index_name(obs))
    l2, l3, l4 = (read_col(obs, k) for k in ("ann_level_2", "ann_level_3", "ann_level_4"))
    sample = read_col(obs, "sample")
    lab = label_v6(l2, l3, l4)

    # genes: as in the prep notebook
    ens = pd.Series(read_col(var, "ensembl_id")).astype(str).str.replace(r"\.\d+$", "", regex=True)
    in_bulk = ens.isin(bulk_genes).to_numpy()
    # notebook order: subset to bulk genes first, then drop duplicates keeping the last
    sub = ens[in_bulk]
    gidx = sub.index[~sub.duplicated(keep="last")].to_numpy()
    print(f"Genes: {len(ens)} in HLCA -> {len(gidx)} shared with bulk (duplicates: keep last)")

    keep = np.array([x is not None for x in lab])
    print(f"Cells: {len(lab)} -> {keep.sum()} with a celltype_v6 label")

    # counts in row chunks (raw counts layer)
    C = f["layers"]["counts"]
    indptr = C["indptr"][:]
    n = len(indptr) - 1
    G = int(C.attrs["shape"][1])
    blocks = []
    for s in range(0, n, CHUNK):
        e = min(n, s + CHUNK)
        a, b = indptr[s], indptr[e]
        M = sp.csr_matrix((C["data"][a:b], C["indices"][a:b], indptr[s:e + 1] - a), shape=(e - s, G))
        blocks.append(M[keep[s:e]][:, gidx].astype(np.float32))
        print(f"  read {e}/{n} cells", flush=True)
    X = sp.vstack(blocks, format="csr")
    f.close()

    obs_df = pd.DataFrame({
        "sample": pd.Categorical(sample[keep]),
        "ann_level_2": pd.Categorical(l2[keep]),
        "ann_level_3": pd.Categorical(l3[keep]),
        "ann_level_4": pd.Categorical(l4[keep]),
        "celltype_v6": pd.Categorical(lab[keep]),
    }, index=pd.Index(cells[keep].astype(str), name=None))
    var_df = pd.DataFrame(index=pd.Index(ens.iloc[gidx].to_numpy(), name=None))
    adata = ad.AnnData(X=X, obs=obs_df, var=var_df)

    # summary
    tab = pd.crosstab(obs_df["celltype_v6"], obs_df["sample"])
    summ = pd.DataFrame({"n_cells": tab.sum(axis=1), "donors_ge10": (tab >= 10).sum(axis=1),
                         "top_donor_share": (tab.max(axis=1) / tab.sum(axis=1)).round(2),
                         "top_donor": tab.idxmax(axis=1)}).sort_values("n_cells", ascending=False)
    print("\ncelltype_v6 summary (filter in run_music2_4.R: >=100 cells, >=20 donors with >=10 cells, top donor <= 0.30):")
    print(summ.to_string())
    summ.to_csv(OUT.replace(".h5ad", "_celltype_summary.csv"))

    if os.path.exists(OLD_REF):
        with h5py.File(OLD_REF, "r") as g:
            old = set(read_col(g["var"], index_name(g["var"])))
        new = set(var_df.index)
        print(f"\nGene check vs old MuSiC reference: old {len(old)}, new {len(new)}, shared {len(old & new)}")

    print("\nWriting", OUT)
    adata.write_h5ad(OUT)
    print("Done:", adata.shape)


if __name__ == "__main__":
    main()
