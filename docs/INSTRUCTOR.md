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

- Send `slides/install_guide.pptx` with the repository link and a hard deadline for
  reporting problems — about three days before day 1. Edit the repo URL and add your
  email to the closing slide first.
- Ask everyone to confirm they have run `00_setup_check.ipynb` successfully. A reply
  saying "done" is worth chasing; silence usually means not started.
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
| 11:00 | Notebook 02 | The core of day 1. Section 4 is a 15-min choose-your-own-thresholds exercise — time-box it, then compare answers out loud |
| 12:00 | Notebook 03, sections 1–2 |  Get them to a clustered object; annotation waits for day 2. Slides 13–20 (normalisation, then PCA/UMAP/Leiden, PCs, resolution, order) belong here rather than in the opening lecture — teach them against the notebook, not cold. **This dataset has a real depth-driven PC1**, so the fix-explorer in notebook 03 is a live demonstration rather than a hypothetical — budget 15 minutes for it |
| 12:25 | Wrap | One question each for tomorrow, written on a sticky note |

### Day 2

| Time | Item | Notes |
|---|---|---|
| 09:30 | Recap + slides 21–24 | 15 min. Answer yesterday's sticky notes |
| 09:45 | Notebook 03, sections 3–5 | Notebook 03 writes `ovarian_annotated.h5ad`, the same name day 2 loads — so anyone who finishes uses their own annotation, and anyone who does not still has the shipped file. Nobody is blocked either way |
| 10:15 | Notebook 04 | Neighbourhood enrichment and niches are the priority; cut Ripley if behind |
| 11:15 | Break | |
| 11:30 | Notebook 05 | The point of the whole workshop. Protect this slot ruthlessly |
| 12:20 | Wrap-up | Go round the room: one sentence each on what they would do with their own tissue. Point them at `docs/DESIGNING_YOUR_STUDY.md` for the design checklist we did not have time for |
| 12:30 | End | |

**The schedule will slip.** Planned sacrifices, in order: Ripley (04), the
segmentation-free section (05.3), the draw-your-own-axis section (05.1b), the manual
contact-null section (05.2). Never cut the distance-field section (05.1) or QC-in-space
(02.3) — they carry the argument of the whole workshop. Never cut
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

## The DEG export

Notebook 03 section 3 writes three files into `results/`: the full DEG table for every
cluster, a top-25-per-cluster version, and an Excel workbook with one sheet per
cluster. People do ask for this, and it is the artefact they take back to their PI.

Two things to say while it runs, both of which apply to every scRNA-seq paper as well:

- **The p-values are anti-conservative.** Clusters were defined from the same data the
  test uses, so the test is not independent of the grouping — every cluster gets
  "highly significant" DEGs, including clusters split out of a homogeneous population.
  Use the ranking and effect sizes; the honest test of whether a cluster is real is
  whether it reappears in another sample.
- **Sort by `above_background`, not by p-value.** The flag comes from notebook 02, and
  a DEG that never beat the negative-control null is not a finding.

The `pct_in` versus `pct_out` columns are the practical ones. A gene in 12% of the
cluster and 8% of everything else can be "significant" with enough cells and is not a
marker.

## Permutation counts in notebook 04

Every method in notebook 04 builds its null by shuffling, and the defaults in the
squidpy docs are too slow for a laptop in a session. We use `n_perms=100` for
neighbourhood enrichment, `n_simulations=20` for Ripley and `N_PERMS=10` for Moran's I.

There is a boxed note explaining what this costs, and it is worth reading aloud rather
than skipping: with 10 permutations the smallest possible p-value is 0.09, so nothing
can be significant. Across 5,000 genes, Bonferroni would need roughly 100,000
permutations before any gene could survive in principle.

The point to land: **the statistic does not depend on n_perms — only the p-value
does.** So the ranking of genes is trustworthy today and the p-values are not, which
is exactly the right tool for choosing what to look at. Someone always asks whether
the results are "real"; this is the honest answer.

If a student's machine is fast and they finish early, `N_PERMS = 1000` on Moran's I is
a reasonable thing for them to set running while they do the exercises.

## Running the choose-your-own-thresholds exercise

Section 4 of notebook 02 is the one place students make a real decision, and it is the
most valuable fifteen minutes of day 1. Run it deliberately.

**Time-box it.** Announce fifteen minutes for steps 1–4 and hold to it. Without a
limit, careful people will still be adjusting numbers at the coffee break.

**Insist on step 3.** The instinct is to pick a number, read the cell count, and move
on. The learning is in plotting what was removed *in space*. Walk the room and ask to
see that plot rather than their thresholds.

**Then compare, out loud.** The last cell prints each person's choice as a single
line. Ask four or five to read theirs out and write them on the board. On our crop the
realistic spread is roughly:

| Choice | kept overall | kept in the tumour nest |
|---|---|---|
| lenient (counts ≥ 5) | ~96% | ~90% |
| middle (counts ≥ 10) | ~91% | ~77% |
| strict (counts ≥ 30, area p5–p95) | ~56% | ~24% |

That third row is the teaching moment: a threshold that sounds merely cautious removes
three quarters of the tumour. Nobody picks it intending that.

**Do not adjudicate.** There is no right answer, and saying so is the point. Push
instead on the methods sentence — if someone cannot write one they would defend in
review, the threshold is not yet a decision.

**Day 2 is safe regardless.** Notebook 04 loads the annotation, and anyone who has not
finished notebook 03 gets the shipped file. Tell them this *before* they start, or the
cautious will pick thresholds that remove nothing.

**If you are running late,** tell them to keep the pre-set numbers and skip to step 3.
The spatial plot alone carries most of the lesson.

## The parameter-exploration section in notebook 03

Section 2 is the notebook-03 equivalent of the QC exercise: students pick `n_pcs`,
`n_neighbors` and `resolution` themselves, then commit. It runs on an 8,000-cell
subsample so each attempt takes seconds.

**There is a `for` loop tutorial in the middle of it.** The resolution sweep is the
first loop most beginners meet, so it is explained properly — anatomy, indentation,
the R comparison — followed by a fill-in-the-blanks exercise where they write one for
`n_neighbors`. The placeholders are named `VALUES_TO_TRY` and `LOOP_VARIABLE`, so an
unfilled blank raises a `NameError` that says which one. Give this three minutes and
walk the room; it is the single most transferable thing in the notebook for someone
who has never coded.

**Budget ten minutes and time-box it**, same as the QC section. The four boxed
try-it cells are meant to be run once each, not optimised.

**The most valuable one is the last.** `min_dist` and `spread` change the UMAP and
nothing else — the cluster count printed above does not move. That is direct evidence
for the claim on slide 16 that a UMAP is a drawing rather than an analysis, and it
lands far better here, on their own screen, than it does from the front.

**On resolution, resist giving a number.** The honest answer is that it depends on the
question: immune subtypes need finer than tissue compartments. Push them towards the
test that does work — can you name every cluster with markers?

**Chosen values are recorded** in `adata.uns["clustering_params"]` and carried into the
saved object, so they end up in a methods section rather than in someone's memory.

**If you are running late,** tell them to keep the defaults in the commit cell and skip
the four try-it boxes. Nothing downstream depends on having explored.

## Two hands-on additions worth protecting

**Colour choice (notebook 03, section 5).** Students set their own palette, stored in
`adata.uns["cell_type_colors"]` so scanpy and the tissue plot both use it, and it
travels into day 2. The exercise that makes the point is the last one: everything grey
except one lineage. Ask the room to compare that against the ten-colour default and
say which they would put in a paper. Five minutes, and it changes how people make
figures afterwards.

**Draw your own axis (notebook 05, section 1b).** Students pick two points off a
coordinate grid and get a signed perpendicular distance for every cell — negative one
side, positive the other — then run the same composition and expression-versus-distance
analyses as section 1.

This generalises the whole of section 1: distance-to-a-nest is one special case, and
an axis you place by eye works for structures no marker labels (a capsule, a lumen, a
fold, a lobe boundary). It is the version they will actually use on their own tissue.

Two things to say while they do it. **Fix the line before looking at the expression
plots** — you can otherwise slide it until a gene looks interesting, and the notebook
says so. And **a straight line is a model**: if the front curves, the ends of the line
mean something different from the middle, and the niche-based distance from section 1
is the better tool.

Ten minutes. If the session is running late this is a reasonable cut, because
section 1 already teaches the concept.

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
