# Spatial transcriptomics: a hands-on introduction

**Department of Anatomy & Embryology, LUMC**
Two mornings, 09:30–12:30. For PhD students and postdocs.

We use a real Xenium Prime 5K run — human ovarian cancer, FFPE — and work from raw
counts to tissue architecture. The through-line of the whole course is one question:

> **What can you learn from a tissue section that you fundamentally cannot learn from
> dissociated cells?**

---

## Before day 1 — do this, it takes 20 minutes

You must arrive with a working environment. There is no time to install software on
the morning; if you turn up unprepared you will spend the first session watching
conda solve.

### 1. Get the code

```bash
git clone https://github.com/<your-org>/spatialtx-workshop.git
cd spatialtx-workshop
```

No git? Download the ZIP from the repository page and unzip it.

### 2. Build the environment

You need conda. If you have none, install [Miniforge](https://github.com/conda-forge/miniforge)
(free, no licence issues, works on Windows/macOS/Linux).

```bash
conda env create -f environment.yml     # or: mamba env create -f environment.yml
conda activate spatialtx
python -m ipykernel install --user --name spatialtx --display-name "Python (spatialtx)"
```

This takes 5–15 minutes. Full details and platform-specific fixes: [`docs/INSTALL.md`](docs/INSTALL.md).

### 3. Copy the data

The dataset lives on the department share:

```
P:\PI\PI_Chuva_de_Sousa_Lopes\susana\SpatialTranscriptomicsWorkshop\data
```

Copy it into your own clone with:

```bash
python scripts/copy_data.py
```

The script finds the P: drive automatically. If P: is not mapped, or you are on a
Mac, point it at wherever the share is mounted:

```bash
python scripts/copy_data.py --source "/Volumes/.../SpatialTranscriptomicsWorkshop/data"
```

**Copy it — do not work directly off P:.** The notebooks write their results back
into `data/`, twenty-five people cannot write to the same network folder at once,
and reading a large `.h5ad` over the network in every cell is painfully slow. You
need roughly 2 GB free.

### 4. Check everything works

```bash
python scripts/check_install.py
```

Then open `notebooks/00_setup_check.ipynb` and run every cell. **If anything fails,
email the organiser with the complete output before day 1.**

---

## Programme

### Day 1 — the data, and whether you can trust it

| Time | | |
|---|---|---|
| 09:30 | Lecture | Why spatial, how Xenium works, what the technologies trade off |
| 10:00 | `01_what_is_in_the_box` | The four data layers; counts, transcripts, polygons, image |
| 10:45 | *break* | |
| 11:00 | `02_quality_control` | Negative controls, segmentation failure, QC in space |
| 12:00 | `03_from_counts_to_cell_types` (start) | Normalisation and clustering for a targeted panel |
| 12:30 | end | |

### Day 2 — using the coordinates

| Time | | |
|---|---|---|
| 09:30 | `03_from_counts_to_cell_types` (finish) | Annotation, and the clusters that are artefacts |
| 10:15 | `04_spatial_statistics` | Spatial graphs, neighbourhood enrichment, Moran's I, niches |
| 11:15 | *break* | |
| 11:30 | `05_beyond_single_cell` | Distance fields, your own axis, contact, segmentation-free |
| 12:20 | Wrap-up | What you would do with your own tissue; questions |
| 12:30 | end | |

---

## What you will be able to do afterwards

- Explain what a Xenium run measures, and what it infers
- Assess data quality using negative controls and segmentation diagnostics, and
  recognise when a "quality" filter is silently deleting biology
- Normalise, cluster and annotate a targeted panel without importing scRNA-seq
  assumptions that do not hold
- Produce a DEG table per cluster, export it, and know which columns to trust
- Build a spatial neighbourhood graph and defend its parameters
- Test cell-type co-localisation against a null that respects tissue architecture
- Identify tissue niches and relate them to structures you can name
- Measure gene expression as a function of distance to an anatomical structure, or to
  an axis you place yourself
- Judge whether a result depends on the segmentation, by re-testing it without cells
- Judge whether a published spatial result is real or an artefact of segmentation or
  of a too-easy null
- Design a spatial experiment, including the hard question of what *n* is
  (see [`docs/DESIGNING_YOUR_STUDY.md`](docs/DESIGNING_YOUR_STUDY.md))

---

## Repository layout

```
notebooks/    the workshop, in order (plus A1 and A2, optional appendices)
scripts/      environment check, notebook check, data copy and preparation
data/         downloaded data lands here (not in git)
solutions/    worked answers — look after trying
docs/         install guide, Python basics, study design, instructor notes
slides/       the lecture deck and the pre-workshop install guide
```

## New to Python?

You will be fine — every cell runs as written, and the workshop is designed for a
very mixed room. When you want to change something, [`docs/PYTHON_BASICS.md`](docs/PYTHON_BASICS.md)
covers the five things you will actually edit and what the common errors mean.

Look for the **Try it yourself** boxes in each notebook. They have one line marked
`<-- CHANGE THIS`, run in a couple of seconds, and print the result. Start there.

## Data

Distributed on the department share at `P:\PI\PI_Chuva_de_Sousa_Lopes\susana\SpatialTranscriptomicsWorkshop\data`.
A spatial crop of the public 10x Genomics demonstration dataset
*Xenium Prime 5K — Human Ovarian Cancer, FFPE*
(`Xenium_Prime_Ovarian_Cancer_FFPE_XRrun_outs`), redistributed for teaching.
Consult the 10x dataset licence terms before publishing anything derived from it.
See `data/data_manifest.yml` for provenance and `scripts/prepare_workshop_data.py`
for exactly how the crop was made.

## Getting help

- During the workshop: raise a hand, or put a red sticky note on your laptop lid
- Afterwards: open an issue on this repository
- Install problems before day 1: email the organiser, and **paste the full output of
  `python scripts/check_install.py`**
- No access to the P: drive: tell the organiser well before day 1 — getting share
  permissions can take a day or two

## Licence

Teaching materials (notebooks, slides, documentation): CC BY 4.0.
Code (`scripts/`): MIT.
The dataset is 10x Genomics' and carries its own terms.
