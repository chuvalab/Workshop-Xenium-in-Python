#!/usr/bin/env python
"""Verify the `spatialtx` environment before the workshop starts.

Run from the repository root:

    conda activate spatialtx
    python scripts/check_install.py

It prints a version table, runs a tiny end-to-end smoke test (build an AnnData,
compute a spatial graph, run one spatial statistic, draw a figure) and reports
whether the workshop data is present. Anything that fails here should be fixed
BEFORE the first morning — not during it.
"""

from __future__ import annotations

import importlib
import platform
import sys
from pathlib import Path

REQUIRED = [
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
    "seaborn",
    "sklearn",
    "skimage",
    "h5py",
    "pyarrow",
    "shapely",
    "geopandas",
    "igraph",
    "leidenalg",
    "anndata",
    "scanpy",
    "squidpy",
]

OPTIONAL = [
    "spatialdata",
    "spatialdata_io",
    "spatialdata_plot",
    "tifffile",
    "datashader",
    "omnipath",
]

DATA_FILES = [
    "data/ovarian_subset.h5ad",
    "data/ovarian_annotated.h5ad",
    "data/transcripts_crop.parquet",
    "data/cell_boundaries_crop.parquet",
    "data/nucleus_boundaries_crop.parquet",
    "data/crop_metadata.json",
]


def _version(name: str) -> str:
    mod = importlib.import_module(name)
    return str(getattr(mod, "__version__", "installed"))


def check_imports() -> bool:
    ok = True
    print("\nRequired packages")
    print("-" * 46)
    for name in REQUIRED:
        try:
            print(f"  {name:<18} {_version(name)}")
        except Exception as exc:  # noqa: BLE001
            print(f"  {name:<18} MISSING  ({exc})")
            ok = False

    print("\nOptional packages")
    print("-" * 46)
    for name in OPTIONAL:
        try:
            print(f"  {name:<18} {_version(name)}")
        except Exception:  # noqa: BLE001
            print(f"  {name:<18} not installed (some cells will be skipped)")
    return ok


def smoke_test() -> bool:
    """Exercise the exact code path notebooks 02-04 depend on."""
    print("\nSmoke test")
    print("-" * 46)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        import scanpy as sc
        import squidpy as sq
        from anndata import AnnData

        rng = np.random.default_rng(0)
        n_cells, n_genes = 600, 40
        counts = rng.poisson(1.2, size=(n_cells, n_genes)).astype("float32")
        adata = AnnData(counts)
        adata.var_names = [f"gene_{i}" for i in range(n_genes)]
        adata.obsm["spatial"] = rng.uniform(0, 1000, size=(n_cells, 2))

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata, n_comps=10)
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata, key_added="cl", flavor="igraph", n_iterations=2)
        print(f"  scanpy       ok  ({adata.obs['cl'].nunique()} clusters)")

        sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)
        sq.gr.nhood_enrichment(adata, cluster_key="cl", show_progress_bar=False)
        sq.gr.spatial_autocorr(adata, mode="moran", n_perms=20, n_jobs=1)
        print("  squidpy      ok  (graph + neighbourhood enrichment + Moran's I)")

        fig, ax = plt.subplots(figsize=(3, 3))
        ax.scatter(*adata.obsm["spatial"].T, s=2, c=adata.obs["cl"].cat.codes)
        fig.savefig("_smoketest.png", dpi=60)
        plt.close(fig)
        Path("_smoketest.png").unlink(missing_ok=True)
        print("  matplotlib   ok  (figure written and removed)")
        return True
    except Exception as exc:  # noqa: BLE001
        print(f"  FAILED: {type(exc).__name__}: {exc}")
        return False


def check_data() -> bool:
    print("\nWorkshop data")
    print("-" * 46)
    root = Path(__file__).resolve().parents[1]
    ok = True
    for rel in DATA_FILES:
        path = root / rel
        if path.exists():
            size = path.stat().st_size / 1e6
            print(f"  {rel:<38} {size:8.1f} MB")
        else:
            print(f"  {rel:<38} MISSING")
            ok = False
    if not ok:
        print("\n  -> run:  python scripts/copy_data.py")
    return ok


def main() -> int:
    print("=" * 46)
    print("Spatial transcriptomics workshop — environment check")
    print("=" * 46)
    print(f"\nPython   {sys.version.split()[0]}  ({sys.executable})")
    print(f"Platform {platform.platform()}")

    imports_ok = check_imports()
    smoke_ok = smoke_test() if imports_ok else False
    data_ok = check_data()

    print("\n" + "=" * 46)
    if imports_ok and smoke_ok and data_ok:
        print("All good. See you on the first morning.")
        return 0
    if imports_ok and smoke_ok:
        print("Environment fine, data missing — run scripts/copy_data.py")
        return 0
    print("Something is wrong. Please email the organiser BEFORE day 1,")
    print("and paste everything printed above into the message.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
