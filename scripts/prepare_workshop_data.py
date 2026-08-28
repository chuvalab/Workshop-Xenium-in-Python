#!/usr/bin/env python
r"""Carve a laptop-sized teaching subset out of a full Xenium Prime 5K run.

INSTRUCTOR SCRIPT — students never run this. Run it once on the machine that
holds `Xenium_Prime_Ovarian_Cancer_FFPE_XRrun_outs/`, then upload the contents
of `data/` to the department share:
    P:\PI\PI_Chuva_de_Sousa_Lopes\susana\SpatialTranscriptomicsWorkshop\data

    python scripts/prepare_workshop_data.py \
        --xenium-dir /path/to/Xenium_Prime_Ovarian_Cancer_FFPE_XRrun_outs \
        --out-dir data \
        --width 1500 --height 1500 \
        --window-width 400 --window-height 400

What it writes
--------------
data/ovarian_subset.h5ad          raw counts + all QC columns, one spatial crop
data/transcripts_crop.parquet     individual transcripts for a small inner window
data/cell_boundaries_crop.parquet cell polygons for that same inner window
data/nucleus_boundaries_crop.parquet
data/morphology_crop.ome.tif      all morphology channels for the inner window
data/crop_metadata.json           exactly which coordinates were taken, and why

Two nested windows
------------------
--width/--height    the CELL crop: every cell in this rectangle goes into the
                    .h5ad. Make it as big as your students' laptops tolerate.
--window-width/     the INNER window, centred inside the cell crop, written out
--window-height     at transcript and pixel resolution. Keep this small: it
                    carries every molecule and every pixel, so it grows fast.
                    Both windows are independent, and neither has to be square.

Why a crop and not a downsample
-------------------------------
Randomly sampling cells would destroy the thing we are here to teach. Every
neighbourhood statistic in notebooks 04 and 05 assumes you have all the cells in
a region, not a random 10% of them. So we take a contiguous rectangle and keep
every cell inside it.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _read_table(directory: Path, stem: str) -> pd.DataFrame:
    """Xenium ships most tables as both .parquet and .csv.gz. Prefer parquet."""
    parquet = directory / f"{stem}.parquet"
    csv = directory / f"{stem}.csv.gz"
    if parquet.exists():
        return pd.read_parquet(parquet)
    if csv.exists():
        return pd.read_csv(csv)
    raise FileNotFoundError(f"Neither {parquet.name} nor {csv.name} in {directory}")


def _decode(series: pd.Series) -> pd.Series:
    """feature_name / cell_id sometimes arrive as bytes from parquet."""
    if series.dtype == object and len(series) and isinstance(series.iloc[0], bytes):
        return series.str.decode("utf-8")
    return series


def choose_window(
    x: np.ndarray, y: np.ndarray, width: float, height: float, step: float = 100.0
) -> tuple[float, float]:
    """Slide a width x height box over the section, return the densest position.

    Densest is a reasonable default for teaching: it maximises the number of
    cells the students get, and in a tumour section the dense regions are also
    the architecturally interesting ones. Override with --x0/--y0 if you have a
    region you specifically want to teach on.
    """
    best, best_n = (float(x.min()), float(y.min())), -1
    xs = np.arange(x.min(), max(x.min(), x.max() - width) + step, step)
    ys = np.arange(y.min(), max(y.min(), y.max() - height) + step, step)
    for x0 in xs:
        in_x = (x >= x0) & (x < x0 + width)
        if not in_x.any():
            continue
        y_sub = y[in_x]
        for y0 in ys:
            n = int(((y_sub >= y0) & (y_sub < y0 + height)).sum())
            if n > best_n:
                best_n, best = n, (float(x0), float(y0))
    print(f"  densest {width:.0f} x {height:.0f} um window holds {best_n:,} cells")
    return best


def crop_mask(df: pd.DataFrame, xcol: str, ycol: str, box: tuple) -> np.ndarray:
    x0, y0, w, h = box
    return (
        (df[xcol] >= x0) & (df[xcol] < x0 + w) & (df[ycol] >= y0) & (df[ycol] < y0 + h)
    ).to_numpy()


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def flag_controls(var, verbose: bool = False):
    """Boolean Series: which features are controls rather than targeted genes.

    `feature_types` from the 10x h5 is authoritative — anything that is not
    "Gene Expression" is a control of some kind (negative control probe,
    negative control codeword, unassigned codeword, genomic control, deprecated
    codeword). Name matching is only a fallback, because the naming has changed
    between Xenium chemistries and silently flagging nothing is far worse than
    flagging something odd.
    """
    import pandas as pd

    names = pd.Index(var.index).str.lower()
    by_name = (
        names.str.startswith("negcontrol")
        | names.str.startswith("neg_control")
        | names.str.startswith("unassignedcodeword")
        | names.str.startswith("unassigned_codeword")
        | names.str.startswith("deprecatedcodeword")
        | names.str.startswith("deprecated_codeword")
        | names.str.startswith("genomiccontrol")
        | names.str.startswith("genomic_control")
        | names.str.startswith("antisense")
        | names.str.startswith("blank")
        | names.str.contains("codeword")
    )
    by_name = pd.Series(by_name, index=var.index)

    if "feature_types" in var.columns:
        ft = var["feature_types"].astype(str)
        by_type = ~ft.str.strip().str.lower().isin(["gene expression", "gene_expression"])
        control = by_type | by_name
        if verbose:
            print("      feature_types breakdown:")
            for k, v in ft.value_counts().items():
                print(f"        {k:<32} {v:>6,}")
    else:
        control = by_name
        if verbose:
            print("      no feature_types column — falling back to name matching")

    if verbose:
        n = int(control.sum())
        print(f"      {n} control features, {int((~control).sum())} targeted genes")
        if n == 0:
            print("      *** WARNING: no control features found. ***")
            print("      Notebook 02 needs them. First 10 feature names:")
            print("        " + ", ".join(map(str, var.index[:10])))
            print("      If these look like Ensembl IDs, reload the matrix with")
            print("      gene symbols as var_names.")
        else:
            print(f"      examples: {', '.join(map(str, var.index[control][:4]))}")
    return control


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--xenium-dir", required=True, type=Path)
    p.add_argument("--out-dir", default=Path("data"), type=Path)
    p.add_argument("--width", type=float, default=1500.0, help="crop width, um")
    p.add_argument("--height", type=float, default=1500.0, help="crop height, um")
    p.add_argument("--x0", type=float, default=None, help="crop origin x, um (default: auto)")
    p.add_argument("--y0", type=float, default=None, help="crop origin y, um (default: auto)")
    p.add_argument(
        "--window-width",
        type=float,
        default=None,
        help="width of the inner window kept at transcript/image resolution, um "
        "(default: 400, or --transcript-window if given)",
    )
    p.add_argument(
        "--window-height",
        type=float,
        default=None,
        help="height of the inner window, um (default: same as --window-width)",
    )
    p.add_argument(
        "--transcript-window",
        type=float,
        default=None,
        help="shorthand: make the inner window square with this side length, um",
    )
    p.add_argument("--skip-image", action="store_true", help="do not crop the morphology OME-TIFF")
    p.add_argument(
        "--image-downsample",
        type=int,
        default=1,
        help="keep every Nth pixel of the morphology crop (default 1 = full "
        "resolution). Use 2-4 for a large window; the students only ever look "
        "at this image, they do not segment from it",
    )
    p.add_argument(
        "--inspect-image",
        action="store_true",
        help="diagnose the morphology files and exit, without running the pipeline",
    )
    args = p.parse_args()

    if args.inspect_image:
        return inspect_image(args.xenium_dir)

    import scanpy as sc

    xdir: Path = args.xenium_dir
    out: Path = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    if not xdir.exists():
        sys.exit(f"No such directory: {xdir}")

    # ---- 1. counts ------------------------------------------------------- #
    print("[1/6] reading cell_feature_matrix.h5 ...")
    # gex_only=False is ESSENTIAL. The default is True, which keeps only
    # feature_type == "Gene Expression" and silently discards every negative
    # control probe, negative control codeword, unassigned codeword and genomic
    # control — i.e. exactly the features notebook 02 measures background with.
    adata = sc.read_10x_h5(xdir / "cell_feature_matrix.h5", gex_only=False)
    adata.var_names_make_unique()
    print(f"      {adata.n_obs:,} cells x {adata.n_vars:,} features")
    if "feature_types" in adata.var.columns:
        print("      feature types in the file:")
        for k, v in adata.var["feature_types"].astype(str).value_counts().items():
            print(f"        {k:<34} {v:>6,}")
    else:
        print("      warning: no feature_types column in the h5 file")

    # ---- 2. per-cell metadata -------------------------------------------- #
    print("[2/6] reading cells table ...")
    cells = _read_table(xdir, "cells")
    cells["cell_id"] = _decode(cells["cell_id"])
    cells = cells.set_index("cell_id")
    adata.obs_names = adata.obs_names.astype(str)
    common = adata.obs_names.intersection(cells.index)
    if len(common) != adata.n_obs:
        print(f"      warning: matched {len(common):,} of {adata.n_obs:,} cell ids")
    adata = adata[common].copy()
    adata.obs = cells.loc[adata.obs_names].copy()
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy(dtype="float32")

    # ---- 3. pick the crop ------------------------------------------------ #
    print("[3/6] choosing crop window ...")
    x = adata.obs["x_centroid"].to_numpy()
    y = adata.obs["y_centroid"].to_numpy()
    if args.x0 is None or args.y0 is None:
        x0, y0 = choose_window(x, y, args.width, args.height)
    else:
        x0, y0 = args.x0, args.y0
    big_box = (x0, y0, args.width, args.height)

    keep = crop_mask(adata.obs, "x_centroid", "y_centroid", big_box)
    sub = adata[keep].copy()
    print(f"      kept {sub.n_obs:,} cells ({100 * keep.mean():.1f}% of the section)")

    # store the panel / control structure explicitly — students need it in nb02
    if "feature_types" in sub.var.columns:
        sub.var["feature_types"] = sub.var["feature_types"].astype(str)
    sub.var["control"] = flag_controls(sub.var, verbose=True)
    if not sub.var["control"].any():
        sys.exit(
            "\nSTOPPING: no control features in the matrix.\n"
            "Notebook 02 measures the background rate from these, so a dataset\n"
            "without them is not usable for the workshop.\n\n"
            "The usual cause is scanpy's read_10x_h5(gex_only=True) default, which\n"
            "drops every non-'Gene Expression' feature. This script passes\n"
            "gex_only=False, so if you still see this, inspect the file directly:\n\n"
            "    import h5py\n"
            f"    f = h5py.File(r'{xdir / 'cell_feature_matrix.h5'}')\n"
            "    set(f['matrix/features/feature_type'][:])\n"
        )

    sub.uns["crop"] = {"x0": x0, "y0": y0, "width": args.width, "height": args.height}
    sub.write_h5ad(out / "ovarian_subset.h5ad", compression="gzip")
    print(f"      -> {out / 'ovarian_subset.h5ad'}")

    # ---- 4. inner window at transcript resolution ------------------------ #
    # The inner window is independent of the cell crop: it is centred inside it,
    # but you set its size yourself. Keep it small — this is the window that gets
    # written out at transcript and pixel resolution.
    wx = args.window_width or args.transcript_window or 400.0
    wy = args.window_height or args.transcript_window or wx
    wx, wy = min(wx, args.width), min(wy, args.height)   # cannot exceed the cell crop

    cx = x0 + (args.width - wx) / 2
    cy = y0 + (args.height - wy) / 2
    small_box = (cx, cy, wx, wy)
    print(f"[4/6] inner window {wx:.0f} x {wy:.0f} um at ({cx:.0f}, {cy:.0f}) ...")
    if min(wx, wy) < 600:
        print("      note: the optional ovrlpy section in notebook 05 fits a")
        print("      transcriptome embedding on this window and wants more tissue.")
        print("      If you plan to teach it, use --window-width/--window-height 800+")

    tx = _read_table(xdir, "transcripts")
    tx["feature_name"] = _decode(tx["feature_name"])
    tx["cell_id"] = _decode(tx["cell_id"])
    tmask = crop_mask(tx, "x_location", "y_location", small_box)
    tx_crop = tx.loc[tmask].reset_index(drop=True)
    tx_crop.to_parquet(out / "transcripts_crop.parquet", index=False)
    print(f"      {len(tx_crop):,} transcripts -> transcripts_crop.parquet")
    del tx

    # ---- 5. boundaries --------------------------------------------------- #
    print("[5/6] boundaries ...")
    for stem, fname in [
        ("cell_boundaries", "cell_boundaries_crop.parquet"),
        ("nucleus_boundaries", "nucleus_boundaries_crop.parquet"),
    ]:
        try:
            b = _read_table(xdir, stem)
        except FileNotFoundError:
            print(f"      {stem} not found, skipping")
            continue
        b["cell_id"] = _decode(b["cell_id"])
        # keep whole polygons: a cell is in if ANY vertex falls in the window
        inside = crop_mask(b, "vertex_x", "vertex_y", small_box)
        ids = b.loc[inside, "cell_id"].unique()
        b_crop = b[b["cell_id"].isin(ids)].reset_index(drop=True)
        b_crop.to_parquet(out / fname, index=False)
        print(f"      {b_crop['cell_id'].nunique():,} polygons -> {fname}")

    # ---- 6. morphology image --------------------------------------------- #
    if args.skip_image:
        print("[6/6] image crop skipped")
    else:
        print("[6/6] cropping morphology image ...")
        try:
            crop_morphology(xdir, out / "morphology_crop.ome.tif", small_box,
                            downsample=args.image_downsample)
        except Exception:  # noqa: BLE001
            import traceback

            print("\n      IMAGE CROP FAILED — full traceback follows.")
            print("      Everything else above was written correctly; only the")
            print("      .ome.tif is missing, and the notebooks tolerate that.")
            print("      To diagnose, run:")
            print(f"        python {Path(__file__).name} --xenium-dir {xdir} --inspect-image")
            print("      " + "-" * 60)
            traceback.print_exc()
            print("      " + "-" * 60 + "\n")

    meta = {
        "source": str(xdir),
        "cell_window_um": {"x0": x0, "y0": y0, "width": args.width, "height": args.height},
        "inner_window_um": {"x0": cx, "y0": cy, "width": wx, "height": wy},
        "n_cells": int(sub.n_obs),
        "n_features": int(sub.n_vars),
        "n_control_features": int(sub.var["control"].sum()),
        "n_transcripts_in_window": int(len(tx_crop)),
    }
    (out / "crop_metadata.json").write_text(json.dumps(meta, indent=2))
    print("\nNote: notebook 04 places the morphology image under the cells using")
    print("crop_metadata.json — ship that file to students along with the rest.")
    print("\nDone. Contents of", out)
    for f in sorted(out.iterdir()):
        print(f"  {f.name:<38} {f.stat().st_size / 1e6:8.1f} MB")
    print("\nNext: copy these to the department share,")
    print(r"  P:\PI\PI_Chuva_de_Sousa_Lopes\susana\SpatialTranscriptomicsWorkshop\data")
    print("then confirm a participant outside your group can read them.")
    return 0


def _morphology_sources(xdir: Path) -> list:
    """Every morphology file, newest layout first."""
    focus = xdir / "morphology_focus"
    if focus.is_dir():
        found = sorted(focus.glob("*.ome.tif")) + sorted(focus.glob("*.ome.tiff"))
        if found:
            return found
    for name in ("morphology_focus.ome.tif", "morphology_mip.ome.tif", "morphology.ome.tif"):
        if (xdir / name).exists():
            return [xdir / name]
    return []


def _open_level(src: Path):
    """Open one morphology file and return (TiffFile, level, closer).

    `is_ome=False` matters. A Prime run writes one file per channel into
    morphology_focus/, but each file's OME header describes the *whole* set, so
    tifffile will otherwise try to assemble a multi-file series — which reads
    files we did not ask for, can fail outright if any is missing, and makes
    `series[0]` a different shape than the file in front of you. Forcing the
    plain-TIFF reader gives us exactly one channel per file, deterministically.
    """
    import tifffile

    tf = tifffile.TiffFile(src, is_ome=False)
    series = tf.series[0]
    levels = getattr(series, "levels", None) or [series]
    return tf, levels[0]


def _read_window(level, rc: tuple, px: float, verbose: bool = True):
    """Read one rectangular window out of a pyramidal TIFF level.

    Tries the cheap strategies first and falls back rather than giving up: on a
    full Xenium image the difference between slicing and reading everything is
    several GB, but a correct slow answer beats no answer.
    """
    import numpy as np

    r0, r1, c0, c1 = rc
    shape = level.shape

    rr0, rr1 = max(0, r0), min(shape[-2], r1)
    cc0, cc1 = max(0, c0), min(shape[-1], c1)
    if rr1 <= rr0 or cc1 <= cc0:
        raise ValueError(
            f"requested window falls outside the image, which is "
            f"{shape[-1] * px:.0f} x {shape[-2] * px:.0f} um "
            f"({shape[-1]} x {shape[-2]} px). Check --x0/--y0 and the window size."
        )
    if verbose and (rr0, rr1, cc0, cc1) != (r0, r1, c0, c1):
        print(f"      note: window clipped to the image edge "
              f"({(cc1 - cc0) * px:.0f} x {(rr1 - rr0) * px:.0f} um)")

    errors = []

    # 1) slice through zarr — only the covering tiles get decoded
    try:
        import zarr

        # A pyramidal TIFF opens as a zarr *group*, one array per level, not as a
        # single array. Asking for level 0 up front gives a plain array store;
        # older tifffile builds ignore the argument, so we also unwrap below.
        try:
            store = level.aszarr(level=0)
        except TypeError:
            store = level.aszarr()

        try:
            z = zarr.open(store, mode="r")

            if not hasattr(z, "shape"):
                # still a group: pick the full-resolution array, which tifffile
                # and the OME-NGFF convention both name "0"
                keys = sorted(z.array_keys()) if hasattr(z, "array_keys") else sorted(z)
                if not keys:
                    raise RuntimeError("zarr group contains no arrays")
                z = z["0"] if "0" in keys else z[keys[0]]

            # guard against having grabbed a downsampled level: the pixel
            # coordinates below are only valid at full resolution
            if tuple(z.shape[-2:]) != tuple(shape[-2:]):
                raise RuntimeError(
                    f"zarr array is {z.shape[-2:]} but the level is {shape[-2:]}; "
                    "refusing to slice at the wrong resolution"
                )

            return np.asarray(z[..., rr0:rr1, cc0:cc1])
        finally:
            store.close()
    except Exception as exc:  # noqa: BLE001
        errors.append(f"zarr slice: {type(exc).__name__}: {exc}")

    # 2) read the whole level, then slice. Correct, but memory-hungry.
    try:
        if verbose:
            n_px = shape[-1] * shape[-2]
            print(f"      note: falling back to a full read of {shape[-1]} x "
                  f"{shape[-2]} px (~{n_px * 2 / 1e9:.1f} GB at 16-bit)")
            print(f"            reason: {errors[0]}")
        return np.asarray(level.asarray()[..., rr0:rr1, cc0:cc1])
    except Exception as exc:  # noqa: BLE001
        errors.append(f"full read: {type(exc).__name__}: {exc}")

    raise RuntimeError("could not read the image window. Tried:\n  - " + "\n  - ".join(errors))


def crop_morphology(xdir: Path, dest: Path, box: tuple, downsample: int = 1) -> None:
    """Crop the morphology image to the inner window, all channels.

    Xenium images are indexed in pixels; the physical pixel size lives in
    experiment.xenium as `pixel_size` (um/px), so micron coordinates are divided
    by it. A Prime run stores DAPI, 18S, ATP1A1, aCD45 and E-cadherin as separate
    files in morphology_focus/; we read each and stack to (C, Y, X).
    """
    import numpy as np
    import tifffile

    sources = _morphology_sources(xdir)
    if not sources:
        raise FileNotFoundError(
            f"no morphology image found under {xdir}. Looked in morphology_focus/ "
            "and for morphology_focus/mip/morphology .ome.tif"
        )

    exp_path = xdir / "experiment.xenium"
    px = 0.2125
    if exp_path.exists():
        try:
            px = float(json.loads(exp_path.read_text()).get("pixel_size", px))
        except Exception:  # noqa: BLE001
            print(f"      note: could not read pixel_size, assuming {px} um/px")
    else:
        print(f"      note: experiment.xenium missing, assuming {px} um/px")

    x0, y0, w, h = box
    c0, r0 = int(x0 / px), int(y0 / px)
    c1, r1 = int(round((x0 + w) / px)), int(round((y0 + h) / px))

    d = max(1, int(downsample))
    est = (r1 - r0) * (c1 - c0) * 2 * len(sources) / d**2
    print(f"      {len(sources)} channel(s), window {c1 - c0} x {r1 - r0} px"
          + (f", downsampled {d}x" if d > 1 else "")
          + f"  ~{est / 1e6:.1f} MB uncompressed")
    if est > 200e6:
        print(f"      -> that is large for a teaching file. Consider "
              f"--image-downsample {min(4, max(2, int((est / 150e6) ** 0.5) + 1))}, "
              "or a smaller --window-width/--window-height.")

    planes, names = [], []
    for n, src in enumerate(sources):
        tf, level = _open_level(src)
        try:
            block = _read_window(level, (r0, r1, c0, c1), px, verbose=(n == 0))
        finally:
            tf.close()
        block = np.asarray(block)
        if d > 1:
            block = block[..., ::d, ::d]
        block = block.reshape(-1, block.shape[-2], block.shape[-1])
        planes.append(block)
        names.extend([src.stem] * block.shape[0])

    crop = np.concatenate(planes, axis=0)
    if crop.shape[0] == 1:
        crop = crop[0]

    dest.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        dest,
        crop,
        photometric="minisblack",
        metadata={
            "axes": "YX" if crop.ndim == 2 else "CYX",
            "PhysicalSizeX": px * d,
            "PhysicalSizeY": px * d,
            "PhysicalSizeXUnit": "\u00b5m",
            "PhysicalSizeYUnit": "\u00b5m",
            "Channels": names,
        },
    )
    ny, nx = crop.shape[-2], crop.shape[-1]
    n_ch = 1 if crop.ndim == 2 else crop.shape[0]
    print(f"      wrote {n_ch} channel(s), {nx} x {ny} px @ {px * d:.4g} um/px "
          f"= {nx * px * d:.0f} x {ny * px * d:.0f} um -> {dest.name}")


def inspect_image(xdir: Path) -> int:
    """Print everything we can see about the morphology files, and stop.

    Run this when the image crop fails. It does not need the rest of the
    pipeline, so it works even if scanpy is unhappy.
    """
    print("=" * 68)
    print("Morphology image inspection")
    print("=" * 68)

    try:
        import tifffile

        print(f"tifffile  {tifffile.__version__}")
    except ImportError:
        print("tifffile  NOT INSTALLED  <- this alone would break the image crop")
        return 1
    try:
        import zarr

        print(f"zarr      {zarr.__version__}")
    except ImportError:
        print("zarr      not installed (crop still works, but reads whole levels)")

    exp = xdir / "experiment.xenium"
    print(f"\nexperiment.xenium: {'found' if exp.exists() else 'MISSING'}")
    if exp.exists():
        try:
            meta = json.loads(exp.read_text())
            print(f"  pixel_size = {meta.get('pixel_size', 'absent')}")
        except Exception as exc:  # noqa: BLE001
            print(f"  could not parse: {exc}")

    focus = xdir / "morphology_focus"
    print(f"\nmorphology_focus/: {'directory' if focus.is_dir() else 'not a directory'}")
    if focus.is_dir():
        for f in sorted(focus.iterdir()):
            print(f"  {f.name:<44} {f.stat().st_size / 1e6:9.1f} MB")

    sources = _morphology_sources(xdir)
    print(f"\nfiles the script would use: {[s.name for s in sources] or 'NONE'}")
    if not sources:
        print("\n-> nothing to read. Is --xenium-dir pointing at the outs/ folder itself?")
        return 1

    for src in sources:
        print(f"\n--- {src.name}")
        for label, kwargs in [("is_ome=False (what we use)", {"is_ome": False}),
                              ("default OME handling", {})]:
            try:
                with tifffile.TiffFile(src, **kwargs) as tf:
                    s = tf.series[0]
                    levels = getattr(s, "levels", None) or [s]
                    print(f"  {label}: {len(tf.series)} series, {len(levels)} level(s)")
                    for i, lv in enumerate(levels[:4]):
                        print(f"      level {i}: shape {lv.shape} dtype {lv.dtype}")
            except Exception as exc:  # noqa: BLE001
                print(f"  {label}: FAILED {type(exc).__name__}: {exc}")

        try:
            tf, level = _open_level(src)
            try:
                probe = _read_window(level, (0, min(64, level.shape[-2]),
                                             0, min(64, level.shape[-1])), 0.2125)
                print(f"  64x64 px test read: OK, shape {probe.shape}, dtype {probe.dtype}")
            finally:
                tf.close()
        except Exception as exc:  # noqa: BLE001
            print(f"  64x64 px test read: FAILED {type(exc).__name__}: {exc}")

    print("\nIf every test read says OK, the image itself is fine and the problem")
    print("is the window: check that --x0/--y0 and the window size fall inside the")
    print("pixel dimensions printed above (multiply px by pixel_size for microns).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
