# Instructor notes

Everything you need to run this twice and have it go better the second time.

---

## Before anything else: prepare the data

Students never touch the full run. You do, once.

```bash
conda activate spatialtx
python scripts/prepare_workshop_data.py \
    --xenium-dir /path/to/Xenium_Prime_Ovarian_Cancer_FFPE_XRrun_outs \
    --out-dir data \
    --width 1500 --height 1500 \
    --window-width 400 --window-height 400
```

Then **look at the crop** before you commit to it. Open the full run in Xenium
Explorer, find the region the script chose (coordinates are printed and stored in
`data/crop_metadata.json`), and check it contains what you want to teach:

- a tumour nest with a clear edge — needed for notebook 05's distance analysis
- some stroma, some vessels, some immune infiltrate
- ideally one region of visibly poorer quality — a fold, an edge, a low-count patch.
  This is *pedagogically valuable*. A pristine crop makes notebook 02 boring and
  teaches students that QC is a formality.

If the auto-selected densest window is not it, pass `--x0` and `--y0` explicitly.

The script prints the estimated image size before reading and suggests a
downsample factor if it is too large. Students only ever *look* at this image —
they never segment from it — so 2x or 3x costs nothing pedagogically.

If the image step fails, run `--inspect-image` (see `docs/INSTALL.md`) before
anything else; it reports what tifffile sees in `morphology_focus/` and whether a
test read succeeds.

There are **two nested windows** and they are set independently. `--width/--height`
is the cell crop that becomes the `.h5ad`; `--window-width/--window-height` is the
smaller inner rectangle written out at transcript and pixel resolution. Neither has
to be square — if your region of interest is a long invasive front, a 900 x 250 µm
inner window is a better teaching object than a square one. `--transcript-window N`
remains as shorthand for a square N x N inner window.

### If you want to teach the ovrlpy section

Section 4 of notebook 05 is **optional** and off the critical path — it is skipped
automatically if `ovrlpy` is not installed. Two things to decide before you commit:

- **Window size.** ovrlpy fits a transcriptome embedding on the transcript window, so
  400 x 400 um is on the small side. Use `--window-width 800 --window-height 800` or
  larger if you plan to run it. That makes `transcripts_crop.parquet` roughly four
  times bigger, so check the file sizes afterwards.
- **Time.** Budget 15 minutes including the discussion, and run it yourself first to
  time `ovr.analyse()` on your window. If day 2 is running late, this is the first
  thing to cut — the segmentation-free section that follows makes an overlapping
  point without an extra dependency.

Worth it if your group works on dense tissue (lymphoid, brain, developing embryo)
where vertical stacking is common, or on anything where tissue folds are a known
problem. Less worth it for sparse, well-spread tissue.

Then generate the day-2 starting point:

```bash
python scripts/make_annotated_checkpoint.py
```

This runs the notebook 02–03 pipeline headlessly and writes `ovarian_annotated.h5ad`.
Students overwrite that file when they finish notebook 03; the first overwrite copies
your version to `ovarian_annotated_reference.h5ad` so they can always get back to it.
If you want the whole room on an identical annotation for a demonstration, tell them
to load the reference file.
**Open it and fix the labels by hand before shipping it.** The automatic
signature-argmax annotation is a reasonable draft and no more; students will take
whatever labels you ship as ground truth, so make them right.

Finally, copy the finished `data/` folder to the department share:

```
P:\PI\PI_Chuva_de_Sousa_Lopes\susana\SpatialTranscriptomicsWorkshop\data
```

Then check the permissions: every participant needs **read** access, and they need
it before day 1, because arranging it can take a day or two. Ask one person outside
your group to run `python scripts/copy_data.py` as a test — you will not discover a
permissions problem from your own account.

If anyone in the group is on macOS, find out how the share mounts on their machine
and add that path to `CANDIDATE_SOURCES` in `scripts/copy_data.py` so it is found
automatically next time.

### Check the control features before you ship anything

The script now prints the `feature_types` breakdown from the h5 and **exits** if it
finds no controls. Compare that breakdown against the run's
`analysis_summary.html` — the numbers per class should match.

**Panel and control counts to expect.** The Prime 5K pre-designed panel is 5,001
genes, plus up to 100 custom add-on genes — so 5,101 features is a full add-on, and
worth finding out which 100 genes were chosen. Controls: 609 negative control
codewords + 40 negative control targets = 649 designated negative controls, plus 21
genomic control probes and a variable number of unassigned and deprecated codewords.
If your numbers do not decompose that way, something is being mislabelled.

**Calibrate the background rate correctly.** A good Prime run puts the background
rate around 1e-4, not the 1e-2 that scRNA-seq intuition suggests. On our crop, ~650
negative controls carry a few dozen counts in total. The consequence for teaching:
the FDR-based detection threshold lands at two or three counts and excludes almost
nothing, so the *binding* constraint on which genes are usable is statistical power,
not background. Notebook 02 computes both thresholds and makes that point explicitly —
it is a better lesson than a clean run that simply passes every check.

Second trap, discovered on this dataset: **deprecated codewords are not controls.**
10x defines them as codewords the onboard pipeline does not use — retired from the
active panel but still in the codebook — and they can carry hundreds of thousands of
counts. Including them in a background estimate inflated it roughly ninefold in our
crop. Notebook 02 now computes background from Negative Control Probe and Negative
Control Codeword only, and reports the rest separately. Expect a question about this;
the top-ten table makes it concrete.

The trap that caused this check to exist: `sc.read_10x_h5()` defaults to
`gex_only=True`, which keeps only `Gene Expression` features and silently discards
every control. The script passes `gex_only=False`. If you ever load a Xenium matrix
yourself, remember to do the same.

### Sizing sanity check

| File | Target |
|---|---|
| `ovarian_subset.h5ad` | 150–400 MB |
| `ovarian_annotated.h5ad` | similar |
| `transcripts_crop.parquet` | < 100 MB |
| boundaries | < 50 MB each |
| `morphology_crop.ome.tif` | < 100 MB — all channels at full resolution, so it grows fast; use `--image-downsample 2` or `3` if the script warns |

Much bigger than this and someone's 8 GB laptop will die in front of the room.
Shrink `--width/--height` and regenerate.

---

## Two weeks before

- Send the repository link with a hard instruction to complete `00_setup_check.ipynb`
  before day 1, and a deadline for reporting problems.
- Ask each participant for **one sentence about their own tissue and question**.
  Collect these. Use them as examples throughout — "for Anna's kidney sections, the
  distance-to-glomerulus version of this would be..." is worth ten generic examples.
- Book a room with power at every seat. Ask participants to copy the data **before**
  they arrive; twenty-five simultaneous copies off the same share will crawl.

### After any notebook edit

```bash
python scripts/check_notebooks.py
```

Static check, under a second, no data needed. It catches syntax errors and — the one
that actually bites — names used before assignment, which is what happens when you
edit or delete a cell that defined something a later cell relies on. A student
running top to bottom in a fresh kernel hits these; you usually do not, because your
kernel still holds the old variable.

It is not a substitute for running the notebooks for real, which is the next item.

## The day before

- Prepare a USB stick with `data/` on it anyway. Someone will turn up without share
  access, or without VPN, and it saves the first half hour.
- Run every notebook top to bottom yourself, in a fresh environment, and record the
  timings. Cell timings drift with package versions.
- Have Xenium Explorer installed on the projected machine with the **full** run
  loaded. Ten minutes of live browsing during the day-1 lecture does more than any
  slide.

---

## Run of show

### Day 1

| Time | Item | Notes |
|---|---|---|
| 09:20 | Doors open | Catch install casualties now; pair them with a working neighbour immediately, do not try to fix conda live |
| 09:30 | Lecture, slides 1–12 | 30 min. Include 5 min of live Xenium Explorer |
| 10:00 | Notebook 01 | Let them read the markdown themselves; talk over the two plotting cells |
| 10:45 | Break | |
| 11:00 | Notebook 02 | The core of day 1. Do not rush the "QC in space" section |
| 12:00 | Notebook 03, sections 1–2 |  Get them to a clustered object; annotation waits for day 2. Slides 13–20 (normalisation, then PCA/UMAP/Leiden, PCs, resolution, order) belong here rather than in the opening lecture — teach them against the notebook, not cold. **This dataset has a real depth-driven PC1**, so the fix-explorer in notebook 03 is a live demonstration rather than a hypothetical — budget 15 minutes for it |
| 12:25 | Wrap | One question each for tomorrow, written on a sticky note |

### Day 2

| Time | Item | Notes |
|---|---|---|
| 09:30 | Recap + slides 21–24 | 15 min. Answer yesterday's sticky notes |
| 09:45 | Notebook 03, sections 3–5 | Notebook 03 writes `ovarian_annotated.h5ad`, the same name day 2 loads — so anyone who finishes uses their own annotation, and anyone who does not still has the shipped file. Nobody is blocked either way |
| 10:15 | Notebook 04 | Neighbourhood enrichment and niches are the priority; cut Ripley if behind |
| 11:15 | Break | |
| 11:30 | Notebook 05 | The point of the whole workshop. Protect this slot ruthlessly. Section 4 (ovrlpy) is optional — cut it first if late |
| 12:15 | Notebook 06 | Challenge 6 (their own tissue) if the group is engaged; challenge 1 if they want more code |
| 12:30 | End | |

**The schedule will slip.** Planned sacrifices, in order: the ovrlpy section (05.4),
Ripley (04), the segmentation-free section (05.5), the manual contact-null
section (05.2). Never cut
the distance-field section (05.1) or QC-in-space (02.4).

---

## Teaching a very mixed room

Assume a third have never used pandas and a third could teach the scanpy parts.

- **Pair deliberately.** Ask on arrival: "who has clustered scRNA-seq data before?"
  and seat mixed pairs. Peer explanation carries the room better than you can.
- **Every cell runs unmodified.** Beginners execute and read; advanced participants
  do the exercises. Never require someone to write code to reach the next cell.
- **Red sticky note on the laptop lid = stuck.** Silent, non-embarrassing, and lets
  you scan the room in two seconds. Green when unstuck.
- **Give the advanced ones the exercises up front**, on day 1. Exercises 2.3, 3.1,
  4.4 and 5.1 are substantial and will absorb a strong participant for the full
  session.
- **Checkpoints are not optional.** At 10:10 on day 2, everyone loads the shipped
  annotated object. Say it out loud: falling behind is expected and costs nothing.
- **Use the Try it yourself boxes for pacing.** Each notebook has two, with a single
  line marked `<-- CHANGE THIS` and a printed result. They take a beginner two to
  three minutes and are the natural thing to say from the front: "spend three minutes
  on the box, then we move on." They also give you a way to keep a slow half of the
  room busy while you help someone individually.
- **Slides 15–16 (PCA/UMAP/Leiden, and what a UMAP will not tell you) are for the
  beginners.** If the room is stronger than expected, run them fast; if a third have
  never clustered anything, they are the most valuable eight minutes of day 1.
- **Point novices at `docs/PYTHON_BASICS.md` on day 1, not at the coffee break.**
  It is one page and covers the errors they will hit.
- **Talk about anatomy, not code.** This room's advantage is that they already know
  what an invasive front is. Every time a plot appears, ask *them* what structure
  they are seeing before you say anything.

---

## Where people reliably get confused

| Symptom | Actual cause | Fix |
|---|---|---|
| "Why is my plot upside down?" | image y-axis convention | `ax.invert_yaxis()`; explain once, early |
| "My cluster numbers differ from the screen" | Leiden is not deterministic across versions | Say so at the start of 03; use the shipped object for anything shared |
| "Is this cell a doublet?" | conflating scRNA-seq doublets with segmentation merges | Notebook 03 §4; go back to the polygon picture in 01 |
| "Everything is significantly co-localised" | global-shuffle null is too easy | Notebook 05 §2 is the answer; flag it in 04 when the first z-score of 300 appears |
| "Should I regress out total counts?" | scRNA-seq reflex | Discuss: here it partly encodes cell size, which is biology |
| "Why no highly variable genes?" | the panel already selected | Notebook 03 §1 callout |
| Silence during exercises | question too open for a beginner | Have a scaffolded version ready: give them the first two lines |

---

## Adapting this to your own dataset

The notebooks are written against Xenium column names, but the structure transfers.

- **Different tissue, same platform:** change the marker dictionary in notebook 03
  and the `TUMOUR`/`CELLTYPE` picks in notebook 05. Everything else runs.
- **MERFISH / CosMx:** the loading in notebook 01 changes; QC concepts are identical,
  though the control-probe structure differs. Check what your platform's blanks are.
- **Visium HD or other binned data:** notebooks 04–05 apply directly to bins.
  Notebook 02's segmentation content becomes a discussion of binning artefacts instead.
- **Adding a second sample** is the single highest-value extension, because it lets
  you teach the *n* problem concretely rather than as a warning. If you have two
  sections, use them.

---

## Feedback

Collect at the end of day 2, on paper, three questions:

1. What is one thing you will do differently in your own project?
2. Which section was too fast, and which too slow?
3. What did you expect to learn and did not?

Question 3 is the one that improves the course.
