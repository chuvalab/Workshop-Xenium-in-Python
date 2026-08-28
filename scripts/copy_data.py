#!/usr/bin/env python
"""Copy the workshop data from the department share into your own folder.

    conda activate spatialtx
    python scripts/copy_data.py

The data lives on the P: drive:

    P:\\PI\\PI_Chuva_de_Sousa_Lopes\\susana\\SpatialTranscriptomicsWorkshop\\data

You need your own copy because the notebooks write files back into `data/`
(the QC'd, annotated and niche objects), and twenty-five people cannot write to
the same folder at once. Copy once, then work locally — it is also much faster
than reading a 500 MB h5ad over the network in every cell.

If P: is not mapped, or you are on macOS or Linux, point the script at wherever
the share is mounted:

    python scripts/copy_data.py --source "/Volumes/PI_Chuva_de_Sousa_Lopes/susana/SpatialTranscriptomicsWorkshop/data"

Re-running is safe: files already present with the right size are skipped.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "data"

# Where the data lives, in order of preference. The first that exists wins.
CANDIDATE_SOURCES = [
    r"P:\PI\PI_Chuva_de_Sousa_Lopes\susana\SpatialTranscriptomicsWorkshop\data",
    # UNC path, for when the drive letter is not mapped but the server is reachable
    r"\\vf-d3-home\divg\PI\PI_Chuva_de_Sousa_Lopes\susana\SpatialTranscriptomicsWorkshop\data",
    # typical macOS / Linux mount points
    "/Volumes/PI_Chuva_de_Sousa_Lopes/susana/SpatialTranscriptomicsWorkshop/data",
    "/Volumes/divg/PI/PI_Chuva_de_Sousa_Lopes/susana/SpatialTranscriptomicsWorkshop/data",
    "/mnt/p/PI/PI_Chuva_de_Sousa_Lopes/susana/SpatialTranscriptomicsWorkshop/data",
]

# Files the notebooks need. Anything else in the share is copied too, but these
# are the ones whose absence is reported as a problem.
REQUIRED = [
    "ovarian_subset.h5ad",
    "ovarian_annotated.h5ad",
    "transcripts_crop.parquet",
    "cell_boundaries_crop.parquet",
    "nucleus_boundaries_crop.parquet",
    "crop_metadata.json",
]
OPTIONAL = ["morphology_crop.ome.tif"]


def find_source(explicit: str | None) -> Path:
    if explicit:
        p = Path(explicit)
        if not p.is_dir():
            sys.exit(
                f"\n--source path does not exist or is not a folder:\n    {p}\n\n"
                "Check the spelling, and that you are connected to the LUMC network "
                "(VPN if you are working from home)."
            )
        return p

    for candidate in CANDIDATE_SOURCES:
        try:
            p = Path(candidate)
            if p.is_dir():
                return p
        except OSError:
            continue

    sys.exit(
        "\nCould not find the data share. Tried:\n  "
        + "\n  ".join(CANDIDATE_SOURCES)
        + "\n\nIf you are on Windows, open File Explorer and check that the P: drive\n"
        "is connected. If you are on macOS or Linux, mount the share and pass the\n"
        "path yourself:\n\n"
        '    python scripts/copy_data.py --source "/path/to/the/share/data"\n'
    )


def copy_one(src: Path, dst: Path) -> tuple[str, float]:
    """Copy with a progress readout. Returns (status, megabytes)."""
    size = src.stat().st_size
    mb = size / 1e6

    if dst.exists() and dst.stat().st_size == size:
        return "already present", mb

    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(dst.suffix + ".part")
    copied, chunk = 0, 4 << 20
    t0 = time.time()
    with src.open("rb") as fin, tmp.open("wb") as fout:
        while True:
            block = fin.read(chunk)
            if not block:
                break
            fout.write(block)
            copied += len(block)
            pct = 100 * copied / size if size else 100
            rate = copied / 1e6 / max(time.time() - t0, 1e-6)
            print(f"\r    {dst.name}: {pct:5.1f}%  ({copied/1e6:6.0f}/{mb:.0f} MB, "
                  f"{rate:.0f} MB/s)", end="", flush=True)
    print()
    tmp.replace(dst)

    if dst.stat().st_size != size:
        return "SIZE MISMATCH", mb
    return "copied", mb


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", default=os.environ.get("SPATIALTX_DATA_SOURCE"),
                    help="folder to copy from (default: look for the P: drive)")
    ap.add_argument("--dest", default=str(DEST), help="where to put it (default: ./data)")
    ap.add_argument("--force", action="store_true", help="re-copy files that already exist")
    args = ap.parse_args()

    src = find_source(args.source)
    dst = Path(args.dest)
    print(f"copying from : {src}")
    print(f"          to : {dst.resolve()}\n")

    available = {p.name: p for p in src.iterdir() if p.is_file()}
    if not available:
        sys.exit(f"The share exists but contains no files:\n    {src}")

    total_mb = sum(p.stat().st_size for p in available.values()) / 1e6
    print(f"{len(available)} file(s), {total_mb:.0f} MB total\n")

    problems = []
    for name in sorted(available):
        target = dst / name
        if args.force and target.exists():
            target.unlink()
        try:
            status, mb = copy_one(available[name], target)
            if status == "copied":
                print(f"  {name:<38} {mb:8.1f} MB  copied")
            elif status == "already present":
                print(f"  {name:<38} {mb:8.1f} MB  already present, skipped")
            else:
                print(f"  {name:<38} {status}")
                problems.append(name)
        except Exception as exc:  # noqa: BLE001
            print(f"  {name:<38} FAILED: {type(exc).__name__}: {exc}")
            problems.append(name)

    print("\nchecking the files the notebooks need")
    missing = [n for n in REQUIRED if not (dst / n).exists()]
    for name in REQUIRED:
        print(f"  {'OK  ' if (dst / name).exists() else 'MISSING'} {name}")
    for name in OPTIONAL:
        state = "OK  " if (dst / name).exists() else "absent (optional)"
        print(f"  {state} {name}")

    print()
    if problems or missing:
        if missing:
            print(f"Missing required file(s): {missing}")
        print("Something went wrong. Re-run the script; if it persists, tell the organiser.")
        return 1

    free = shutil.disk_usage(dst).free / 1e9
    print(f"All data copied. {free:.1f} GB free on this disk.")
    print("Next: python scripts/check_install.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
