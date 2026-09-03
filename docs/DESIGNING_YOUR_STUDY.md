# Designing your own spatial experiment

Not covered in the workshop — there was not time — but this is the part you need when
you plan your own study. Read it before you book instrument time.

## The thinking checklist

Work through these for **your** tissue, before you spend a euro.

### 1. Is the question spatial?
Write your question in one sentence. If you can answer it from dissociated cells,
do that instead — it is cheaper, deeper, and better powered. Spatial is worth it when
the answer contains one of: *next to*, *distance from*, *within the region of*,
*organised as*, *gradient*, *interface*, *architecture*.

### 2. What is n?
The commonest fatal design error in this field. Cells are not independent replicates.
Two cells 10 um apart share everything: patient, section, staining batch, focus. Your
n is the number of **sections**, and really the number of **donors**.

A study with 400,000 cells from two patients is an n of 2. Plan the statistics
accordingly: aggregate to a per-sample value, then test across samples. If that
leaves you with n=3 per group, say so at the design stage rather than discovering it
at review.

### 3. Panel or whole transcriptome?
- Imaging (Xenium, MERFISH, CosMx): a fixed panel, single-cell resolution,
  sub-cellular detail.
- Sequencing (Visium HD, Stereo-seq, Slide-seq): unbiased transcriptome, but
  resolution is a bin, not a cell, and cell boundaries must be inferred.

If you cannot list the genes that answer your question, you are not ready for a
panel — do a pilot with an unbiased method or use an existing atlas first.

### 4. Does the panel contain your answer?
Take your marker list to `gene_panel.json` and check every gene. Missing markers are
the single most common regret. Add-on custom panels exist; budget for them early.

### 5. Section geometry
A tissue section is a 2D slice through a 3D structure. A tubular gland cut
transversely looks like a ring; cut longitudinally, a stripe. Any "distance to
structure" measurement is a distance *within the plane*. Decide the orientation
deliberately and record it — as an anatomy department this is the thing you are best
placed in the building to get right.

### 6. Controls
- A tissue region you know the answer for, on the same slide.
- Where possible, more than one section per block, so you can separate section
  effects from biology.
- Adjacent H&E or IF, for orientation and for arguing with reviewers.

### 7. Storage and compute
One Xenium Prime run is roughly 20–100 GB. Plan storage before you book instrument
time. You will not be doing this on a laptop — arrange server or HPC access at the
same time as the experiment.

## Where to go next

**In this repository**
- `A1_controls_and_detection_limits.ipynb` — the control classes in full, and how to
  turn negative controls into a gene-level detection threshold. Worth an hour before
  you design your own panel.
- `A2_normalisation_and_depth.ipynb` — four ways to normalise a targeted panel, a
  diagnostic for whether a principal component is measuring depth, and an explorer
  that runs six remedies. This dataset has a real depth artefact, so it is worked
  through on live data rather than a toy example.

**Read**
- Moses & Pachter (2022), *Nature Methods* — museum of spatial transcriptomics; the
  best single orientation to the field and its history.
- Palla, Fischer, Regev & Theis (2022), *Nature Biotechnology* — spatial components
  of molecular tissue biology; the conceptual framing used in this workshop.
- The squidpy and SpatialData documentation tutorials. They are unusually good.
- Bhatt & Bhatt, and the CellCharter / UTAG / Banksy papers, for niche detection
  beyond the k-means recipe in notebook 04.

**Tools worth knowing beyond this workshop**
- **Xenium Explorer** (free, 10x) — browse the full run interactively with the images.
  Use it constantly; it will catch artefacts no plot will.
- **Baysor**, **Proseg**, **Cellpose/Segger** — alternative segmentation. If your
  finding changes when you re-segment, it was a segmentation finding.
- **CellCharter**, **Banksy**, **UTAG** — niche and domain detection.
- **SpatialData** — the container, especially once you have multiple sections or
  modalities to align.
- **Voyager** (R) — spatial statistics with a proper geospatial foundation, if you
  prefer R or want variograms and local Moran's I.

**Local**
- The LUMC Sequencing Analysis Support Core (SASC) and the Cell Observatory for
  instrument access and study design consultation. Talk to them *before* the
  experiment, not after.

**In this repository**
- `A1_controls_and_detection_limits.ipynb` — control classes, and how to turn negative
  controls into a gene-level detection threshold. Worth an hour before you commit to a
  panel.
- `A2_normalisation_and_depth.ipynb` — normalisation options for a targeted panel, and
  a diagnostic for whether a principal component is measuring depth rather than
  biology.

---

---

The mistake to avoid is treating a spatial dataset as scRNA-seq with a bonus scatter
plot at the end. The coordinates are not decoration — they are the measurement that
justified the cost. If your figure 4 would look the same without them, redesign the
analysis.