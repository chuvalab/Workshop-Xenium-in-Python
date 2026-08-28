#!/usr/bin/env python
"""Download the workshop data files listed in data/data_manifest.yml.

    conda activate spatialtx
    python scripts/download_data.py

Re-running is safe: files that already exist with the right checksum are
skipped. If you are on the LUMC network and the download is slow, ask the
organiser for the USB stick instead — just copy its contents into data/.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import requests
import yaml

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "data_manifest.yml"
CHUNK = 1 << 20


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        done = 0
        with tmp.open("wb") as fh:
            for chunk in r.iter_content(CHUNK):
                fh.write(chunk)
                done += len(chunk)
                if total:
                    pct = 100 * done / total
                    print(f"\r    {dest.name}: {pct:5.1f}%  ({done/1e6:.0f}/{total/1e6:.0f} MB)",
                          end="", flush=True)
        print()
    tmp.replace(dest)


def main() -> int:
    if not MANIFEST.exists():
        sys.exit(f"Manifest not found: {MANIFEST}")
    manifest = yaml.safe_load(MANIFEST.read_text())

    files = manifest.get("files", [])
    if not files:
        sys.exit("Manifest lists no files.")

    print(f"{len(files)} file(s) to fetch into {ROOT / 'data'}\n")
    failures = 0
    for entry in files:
        name, url = entry["name"], entry["url"]
        want = entry.get("sha256")
        dest = ROOT / "data" / name

        if url.startswith("REPLACE_ME"):
            print(f"  {name}: manifest still has a placeholder URL — ask the organiser")
            failures += 1
            continue

        if dest.exists() and (want is None or sha256(dest) == want):
            print(f"  {name}: already present, skipping")
            continue

        print(f"  {name}: downloading")
        try:
            download(url, dest)
        except Exception as exc:  # noqa: BLE001
            print(f"    failed: {type(exc).__name__}: {exc}")
            failures += 1
            continue

        if want:
            got = sha256(dest)
            if got != want:
                print(f"    checksum mismatch!\n      expected {want}\n      got      {got}")
                failures += 1
            else:
                print("    checksum ok")

    print()
    if failures:
        print(f"{failures} file(s) had problems. Re-run, or ask for the USB stick.")
        return 1
    print("All data present. Now run: python scripts/check_install.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
