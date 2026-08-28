#!/usr/bin/env python
"""Produce data/ovarian_annotated.h5ad — the day-2 starting point.

INSTRUCTOR SCRIPT. Runs the notebook 02–03 pipeline headlessly so that every
student on day 2 begins from an identical, correctly annotated object regardless
of how far they got on day 1.

    python scripts/make_annotated_checkpoint.py

IMPORTANT: the automatic annotation below assigns each cluster the marker
signature it scores highest on. That is a draft. Open the result, check it
against the markers and against the tissue, and correct the labels by hand
before you ship it — students will treat these labels as ground truth.

Edit MARKERS for a different tissue.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

# Primary markers plus alternates. A signature is a claim about a cell type, not
# about specific genes: if COL1A1 is absent from the panel, another fibrillar
# collagen makes the same claim. Alternates are only used when too few primaries
# survive the panel check.
MARKERS = {
    "Tumour / epithelial": (
        ["EPCAM", "CP", "LAPTM4B", "NDRG1", "TFPI2", "H19", "UCHL1", "PLXNB1", "CD47"],
        ["KRT8", "KRT18", "KRT19", "MUC1", "CLDN4"]),
    "Tumour, proliferating": (
        ["TOP2A", "BIRC5", "SMC4", "HNRNPD"],
        ["MKI67", "PCNA", "CCNB1", "CDK1", "UBE2C"]),
    "Ciliated epithelium": (
        ["FAM183A", "CFAP100"],
        ["FOXJ1", "PIFO", "TPPP3", "CAPS", "DNAI1"]),
    "Fibroblast / CAF": (
        ["DCN", "LUM", "POSTN", "BGN", "C7", "OGN"],
        ["ASPN", "FMOD", "MGP", "SPARC", "FBLN1"]),
    "Fibroblast, adventitial": (
        ["PI16", "MFAP5", "TIMP3"],
        ["SCARA5", "CD34", "PCOLCE2"]),
    "Fibroblast, matrix-high": (
        ["COL1A1", "COL1A2", "COL4A1", "COL4A2"],
        ["COL3A1", "COL5A1", "COL6A1", "COL6A2", "COL6A3", "FN1"]),
    "Smooth muscle": (
        ["MYH11", "MYL9", "TAGLN", "C11orf96"],
        ["ACTA2", "CNN1", "DES", "ACTG2"]),
    "Endothelial": (
        ["FLT1", "EPAS1", "AQP1", "PECAM1", "SOCS3", "TFPI"],
        ["VWF", "CDH5", "CLDN5", "RAMP2", "EGFL7"]),
    "T cell": (
        ["TRAC", "TRBC1", "CD52", "CXCR4"],
        ["CD3D", "CD3E", "CD2", "IL7R", "CD8A"]),
    "Myeloid / macrophage": (
        ["C1QC", "MS4A6A", "AIF1", "FCGR3A", "FCGBP"],
        ["CD68", "CD14", "C1QA", "C1QB", "LYZ", "TYROBP"]),
}


def build_signature(var_names, primary, alternates=()):
    """Genes on the panel, topped up from alternates if too few survive."""
    lower = {g.lower(): g for g in var_names}
    hit = lambda n: n if n in var_names else lower.get(n.lower())
    found = [g for g in (hit(n) for n in primary) if g]
    missing = [n for n in primary if hit(n) is None]
    used_alt = []
    if len(found) < 3:
        for n in alternates:
            g = hit(n)
            if g and g not in found:
                found.append(g)
                used_alt.append(g)
            if len(found) >= 4:
                break
    return found, missing, used_alt


# QC thresholds — keep in sync with notebook 02
MIN_COUNTS = 10
MIN_GENES = 5
MAX_CONTROL_FRAC = 0.05   # on true negative controls this is very lax
LEIDEN_RES = 1.0


def main() -> int:
    src = DATA / "ovarian_subset.h5ad"
    if not src.exists():
        sys.exit(f"{src} not found — run scripts/prepare_workshop_data.py first")

    print("loading", src)
    adata = sc.read_h5ad(src)
    adata.layers["counts"] = adata.X.copy()

    # ---- controls ------------------------------------------------------- #
    # recompute rather than trusting a stored column: feature_types is
    # authoritative and older prep runs flagged by name only
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from prepare_workshop_data import flag_controls

    adata.var["control"] = flag_controls(adata.var, verbose=True)
    is_ctrl = adata.var["control"].to_numpy()
    if not is_ctrl.any():
        print("WARNING: no control features found; control_frac will be 0 for all cells")

    # Background must come from the TRUE negative controls only. Deprecated
    # codewords are retired probes that still detect real transcripts, so
    # including them turns control_frac into a measure of biology — and any
    # threshold on it then removes cells for expressing the wrong genes.
    BACKGROUND_CLASSES = ["Negative Control Probe", "Negative Control Codeword"]
    ftypes = (adata.var["feature_types"].astype(str)
              if "feature_types" in adata.var.columns
              else pd.Series("unknown", index=adata.var_names))
    is_background = ftypes.isin(BACKGROUND_CLASSES).to_numpy() & is_ctrl
    if not is_background.any():
        print("WARNING: no designated negative controls; using all non-gene features")
        is_background = is_ctrl
    print(f"background from {int(is_background.sum())} negative-control features")

    adata.obs["control_counts"] = np.asarray(
        adata.layers["counts"][:, is_background].sum(axis=1)
    ).ravel()
    adata = adata[:, ~is_ctrl].copy()

    # ---- QC -------------------------------------------------------------- #
    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True, log1p=False)
    adata.obs["control_frac"] = adata.obs["control_counts"] / (
        adata.obs["total_counts"] + adata.obs["control_counts"]
    )
    lo, hi = np.percentile(adata.obs["cell_area"], [1, 99])
    keep = (
        (adata.obs["total_counts"] >= MIN_COUNTS)
        & (adata.obs["n_genes_by_counts"] >= MIN_GENES)
        & (adata.obs["control_frac"] <= MAX_CONTROL_FRAC)
        & (adata.obs["cell_area"].between(lo, hi))
    ).to_numpy()
    print(f"QC keeps {keep.sum():,} / {adata.n_obs:,} cells ({100*keep.mean():.1f}%)")
    adata = adata[keep].copy()

    # ---- normalise, cluster ---------------------------------------------- #
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=50, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, key_added="leiden",
                 flavor="igraph", n_iterations=2)
    adata.X = adata.layers["lognorm"].copy()
    print(f"{adata.obs['leiden'].nunique()} Leiden clusters")

    # ---- annotate (draft) ------------------------------------------------ #
    present = {}
    for name, (primary, alternates) in MARKERS.items():
        found, missing, used_alt = build_signature(adata.var_names, primary, alternates)
        note = ""
        if missing:
            note += f"   missing: {missing}"
        if used_alt:
            note += f"   alternates: {used_alt}"
        print(f"  {name:<26} {len(found)} genes{note}")
        if found:
            present[name] = found
    for name, genes in present.items():
        sc.tl.score_genes(adata, genes, score_name=f"score_{name}")
    cols = [c for c in adata.obs.columns if c.startswith("score_")]
    means = adata.obs.groupby("leiden", observed=True)[cols].mean()
    auto = means.idxmax(axis=1).str.replace("score_", "", regex=False)

    # ------------------------------------------------------------------ #
    # HAND CORRECTIONS — edit this dict after inspecting the result.
    # e.g. OVERRIDES = {"7": "Tumour (proliferating)", "13": "Segmentation artefact"}
    OVERRIDES: dict[str, str] = {}
    # ------------------------------------------------------------------ #
    mapping = {**auto.to_dict(), **OVERRIDES}
    adata.obs["cell_type"] = adata.obs["leiden"].map(mapping).astype("category")

    print("\ndraft annotation:")
    print(adata.obs["cell_type"].value_counts().to_string())
    if not OVERRIDES:
        print("\n*** OVERRIDES is empty — inspect and correct before shipping ***")

    dest = DATA / "ovarian_annotated.h5ad"
    adata.write_h5ad(dest, compression="gzip")
    print(f"\nwrote {dest}  ({dest.stat().st_size/1e6:.0f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
