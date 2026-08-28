# Installation

Twenty minutes, most of it waiting. Do it at least a few days before the workshop so
there is time to fix problems.

## What you need

- A laptop with **at least 8 GB RAM** (16 GB is comfortable) and **~5 GB free disk**
- conda — [Miniforge](https://github.com/conda-forge/miniforge) is the recommended
  installer. It defaults to conda-forge and avoids the Anaconda commercial licence
  question entirely, which matters for institutional use.
- Windows users: run everything in **Miniforge Prompt**, not PowerShell or CMD.

## Steps

```bash
git clone https://github.com/<your-org>/spatialtx-workshop.git
cd spatialtx-workshop

conda env create -f environment.yml
conda activate spatialtx

python -m ipykernel install --user --name spatialtx --display-name "Python (spatialtx)"

python scripts/copy_data.py
python scripts/check_install.py

jupyter lab
```

`mamba env create -f environment.yml` is several times faster and solves identically.
Recent conda versions use the same solver, so plain conda is fine too.

---

## Troubleshooting

### The solve takes forever, or fails with a long conflict message

Use mamba:

```bash
conda install -n base -c conda-forge mamba
mamba env create -f environment.yml
```

If it still conflicts, the usual culprit is an old conda. Check with `conda --version`
(want ≥ 23.10) and `conda update -n base conda`.

### `PackagesNotFoundError`

Almost always a channel problem. This environment is conda-forge only. Check:

```bash
conda config --show channels
```

If `defaults` appears first, that is likely the cause. Either install Miniforge, or:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### macOS with Apple Silicon (M1/M2/M3/M4)

The environment installs natively on arm64 and is faster than under Rosetta. If you
have an old x86 conda installation, remove it and install the arm64 Miniforge.

If a package genuinely has no arm64 build, force an x86 environment:

```bash
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
conda activate spatialtx
conda config --env --set subdir osx-64
```

### Windows: `Failed building wheel for ...` during the pip step

The pip section installs the spatialdata family. If a pure-Python package tries to
build from source, install the Microsoft C++ Build Tools, or skip it — spatialdata is
optional and only used for a demonstration cell in notebook 01. Everything else runs
without it.

### The `spatialtx` kernel does not appear in Jupyter

You skipped the `ipykernel install` line, or ran it in the wrong environment:

```bash
conda activate spatialtx
python -m ipykernel install --user --name spatialtx --display-name "Python (spatialtx)"
```

Restart JupyterLab afterwards. Verify inside a notebook with:

```python
import sys; print(sys.executable)   # must contain /envs/spatialtx/
```

### `sc.tl.leiden` warns about the flavour, or errors

Recent scanpy deprecated the old `leidenalg` path. The notebooks use
`flavor="igraph", n_iterations=2` throughout, which is the current recommendation. If
your scanpy is older than 1.10 that argument does not exist — update the environment:

```bash
conda env update -f environment.yml --prune
```

### Kernel dies while loading data

Memory. Close other applications and Chrome tabs. The teaching crop is deliberately
sized for 8 GB, but a browser with forty tabs will beat you to it.

### `ovrlpy` is missing

It is optional — section 4 of notebook 05 detects it and skips itself if absent.
To install it into the workshop environment:

```bash
conda activate spatialtx
pip install ovrlpy
```

It is also on bioconda (`conda install bioconda::ovrlpy`) if you prefer.

### Notebook 02: "Mean of empty slice" / no control features

`adata.var["control"]` is all False, so there is nothing to compute a background
rate from. Notebook 02 now recomputes the flag itself using `feature_types`, which
is authoritative, so re-running the cell should fix it. If it still reports zero,
the diagnostic it prints tells you what is in `var` — send that to the organiser.

Instructors: regenerate the data with the current
`scripts/prepare_workshop_data.py`, which prints a `feature_types` breakdown and
warns loudly if it finds no controls.

### The image crop fails, or `morphology_crop.ome.tif` is missing

Run the built-in diagnostic — it does not need scanpy, so it works even when the
rest of the pipeline does not:

```bash
python scripts/prepare_workshop_data.py --xenium-dir /path/to/outs --inspect-image
```

It prints which files were found, what tifffile makes of each, the pixel dimensions,
and the result of a 64x64 px test read. If every test read says OK, the image is fine
and your window is outside it: compare `--x0/--y0` plus the window size against the
pixel dimensions multiplied by `pixel_size`.

Note that the pipeline no longer hides this failure. If the crop fails you get a full
traceback, and every other output file is still written correctly.

### `tifffile` cannot read the morphology image

You are missing `imagecodecs` (the OME-TIFFs are JPEG-2000 compressed):

```bash
conda install -c conda-forge imagecodecs
```

Nothing downstream depends on the image, so you can also just skip that cell.

### The data copy fails, or P: is not there

The data is on the department share:

```
P:\PI\PI_Chuva_de_Sousa_Lopes\susana\SpatialTranscriptomicsWorkshop\data
```

- **"Could not find the data share"** — on Windows, open File Explorer and check the
  P: drive is connected; reconnect it if not. Working from home means you need the
  LUMC VPN.
- **macOS or Linux** — mount the share, then pass the path:
  `python scripts/copy_data.py --source "/Volumes/.../SpatialTranscriptomicsWorkshop/data"`.
  You can also set `SPATIALTX_DATA_SOURCE` once instead of typing it every time.
- **"Permission denied"** — you may not have access to that share yet. Ask the
  organiser; this can take a day or two to arrange, so do not leave it to day 1.
- **The copy stops partway** — just re-run it. Files already present with the right
  size are skipped, so it picks up where it left off.
- **"No space left"** — you need roughly 2 GB free, plus room for the derived objects
  the notebooks write.

Do not run the notebooks directly against P:. They write results back into `data/`,
and a shared network folder is both slow and unwritable for a room of people.

---

## Nuclear option

Start clean:

```bash
conda deactivate
conda env remove -n spatialtx
conda clean --all
conda env create -f environment.yml
```

## Still stuck

Email the organiser **before day 1** with:

```bash
conda --version
conda info
python scripts/check_install.py
```

and paste all of it. Screenshots of error messages are much harder to help with than
pasted text.
