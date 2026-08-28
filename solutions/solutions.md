# Worked solutions

Try the exercise first. Most of the value is in the attempt; several of these have
more than one defensible answer, and the reasoning matters more than the code.

Assumes the notebook's variables are already in memory.

---

## 1.1 — Which control feature has the highest count?

```python
tot = np.asarray(adata.layers["counts"].sum(axis=0)).ravel() if "counts" in adata.layers \
      else np.asarray(adata.X.sum(axis=0)).ravel()
ctrl = pd.Series(tot[adata.var["control"].to_numpy()],
                 index=adata.var_names[adata.var["control"]]).sort_values(ascending=False)
print(ctrl.head())
```

**Interpretation.** One control standing far above the rest is usually a specific
probe cross-hybridising, or a codeword that is one bit away from a highly expressed
gene. That is informative rather than alarming. What *would* be alarming is the whole
control distribution shifting upwards, which indicates systematic background.

Compare against the gene distribution before deciding: if the top control sits below
the 10th percentile of gene counts, the run is clean.

---

## 1.2 — Fraction of transcripts assigned to a cell

```python
tx = pd.read_parquet(DATA / "transcripts_crop.parquet")
unassigned = tx["cell_id"].astype(str).isin(["UNASSIGNED", "-1", "0"])
print(f"unassigned: {unassigned.mean():.1%}")

by_gene = (tx.assign(un=unassigned)
             .groupby("feature_name")
             .agg(n=("un", "size"), frac_unassigned=("un", "mean"))
             .query("n >= 200")
             .sort_values("frac_unassigned"))
print(by_gene.head(10))   # most reliably captured
print(by_gene.tail(10))   # most often lost
```

**Interpretation.** 10–40% unassigned is normal. Genes with a *high* unassigned
fraction tend to be secreted or extracellular-matrix transcripts, which genuinely sit
outside cell boundaries, and transcripts from cells with thin cytoplasm. The point
for students: unassigned is not uniform across genes, so segmentation loss is a
**gene-specific bias**, not just a loss of sensitivity.

---

## 2.1 — Where do your genes of interest sit?

```python
genes = ["TOP2A", "TRAC", "CXCR4"]
tot = np.asarray(adata.layers["counts"].sum(axis=0)).ravel()
s = pd.Series(tot, index=adata.var_names)
ctrl_q = np.percentile(tot[adata.var["control"].to_numpy()], [50, 95, 99]) \
         if "control" in adata.var else None

for g in genes:
    if g in s.index:
        print(f"{g:<8} {s[g]:>8,.0f} counts   control p50/p95/p99 = {ctrl_q}")
```

**Rule of thumb.** Below the control 95th percentile: not measurable, do not use.
Within a few-fold of it: usable as a cluster-level mean, never per cell. An order of
magnitude above: fine per cell. `CXCR4` and other low-abundance immune genes are very often in the "cluster mean only"
band, which is exactly why single-cell exhaustion claims from spatial panels need care.

---

## 2.2 — A region to exclude wholesale

```python
x, y = adata.obsm["spatial"].T
# adjust to whatever your crop shows
box = dict(x0=x.min(), x1=x.min() + 200, y0=y.min(), y1=y.max())
m = (x >= box["x0"]) & (x < box["x1"]) & (y >= box["y0"]) & (y < box["y1"])
print(f"{m.sum():,} cells ({100*m.mean():.1f}%) in the excluded strip")
print(adata.obs.loc[m, ["total_counts", "control_frac"]].median())
```

**Methods sentence.** "Cells within 200 µm of the section edge were excluded prior to
analysis because median transcript counts in this band were <X> compared with <Y>
elsewhere, consistent with tissue-edge damage."

The teaching point is that this is a *stated, spatially defined* exclusion rather
than an invisible consequence of a count threshold.

---

## 2.3 — Stricter filtering

```python
for thr in (10, 25):
    k = (adata.obs["total_counts"] >= thr).to_numpy()
    print(f"min_counts={thr}: keep {k.sum():,} ({100*k.mean():.1f}%)")
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x[k], y[k], s=0.6, c="0.85", linewidths=0, rasterized=True)
    ax.scatter(x[~k], y[~k], s=1.4, c="#F26B43", linewidths=0, rasterized=True)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_title(f"removed at {thr}")
    plt.show()
```

**Expected outcome.** The cells lost between 10 and 25 are disproportionately small,
low-cytoplasm cells — lymphocytes above all. Raising the threshold therefore depletes
the immune compartment preferentially and will change every immune co-localisation
result downstream. This is a real, published source of disagreement between studies.

---

## 3.1 — Does normalisation change the answer?

```python
def cluster_from(adata, layer=None, use_rep=None, key="leiden_alt"):
    b = adata.copy()
    if layer is not None:
        b.X = b.layers[layer].copy()
        sc.pp.scale(b, max_value=10)
        sc.pp.pca(b, n_comps=50)
        sc.pp.neighbors(b, n_neighbors=15, n_pcs=30)
    else:
        sc.pp.neighbors(b, n_neighbors=15, use_rep=use_rep)
    sc.tl.leiden(b, resolution=1.0, key_added=key, flavor="igraph", n_iterations=2)
    return b.obs[key].values

adata.obs["cl_area"] = cluster_from(adata, layer="lognorm_area", key="cl_area")
adata.obs["cl_pearson"] = cluster_from(adata, use_rep="X_pca_pearson", key="cl_pearson")

ct = pd.crosstab(adata.obs["cell_type"], adata.obs["cl_area"], normalize="index")
sns.heatmap(ct, cmap="Blues"); plt.xlabel("area-normalised clusters")
plt.ylabel("count-normalised annotation"); plt.show()

from sklearn.metrics import adjusted_rand_score as ari
print("median vs area   :", round(ari(adata.obs["leiden"], adata.obs["cl_area"]), 3))
print("median vs pearson:", round(ari(adata.obs["leiden"], adata.obs["cl_pearson"]), 3))
```

**The rule.** Normalisation matters most where cell **size** differs systematically
between the populations you are comparing. Lymphocytes (tiny) versus tumour cells
(large) are the worst case: count-normalisation partly encodes size, so the two
schemes disagree there. Populations of similar size — different fibroblast states,
say — are stable either way. State which you used, and check your main claim under
both.

**Do the schemes disagree about the same cells?** Usually yes, and that is the useful
part. Pull out the cells where the labelings differ and look at where they sit:

```python
moved = (adata.obs["leiden"].astype(str) + "|" + adata.obs["cl_area"].astype(str))
unstable = moved.map(moved.value_counts()) < 50          # rare label combinations
x, y = adata.obsm["spatial"].T
plt.scatter(x, y, s=0.5, c="0.88"); plt.scatter(x[unstable], y[unstable], s=2, c="#F26B43")
plt.gca().set_aspect("equal"); plt.gca().invert_yaxis(); plt.show()
```

If the unstable cells cluster at compartment boundaries, the instability is
segmentation-driven — those are the cells whose polygons are least trustworthy, so
their size factors are least trustworthy too.

## 3.2 — Is your smallest cluster real?

Three lines of evidence, at least one spatial:

1. **Markers.** Coherent and lineage-consistent, or a grab-bag from several lineages?
2. **QC.** Does it sit at an extreme of `total_counts` or `cell_area`? Technical.
3. **Space.** Scattered at random → likely artefact or doublet. Concentrated in a
   structure → likely real, and the structure names it.

```python
small = adata.obs["cell_type"].value_counts().idxmin()
m = (adata.obs["cell_type"] == small).to_numpy()
print(adata.obs.loc[m, ["total_counts", "cell_area", "control_frac"]].median())
print(adata.obs.loc[m, "niche"].value_counts(normalize=True).head() if "niche" in adata.obs else "")
```

Concentration in one niche is strong evidence for a real population.

---

## 3.3 — Profile a marker gene instead

```python
genes = [g for g in ["TRAC", "C1QC", "COL1A1"] if g in adata.var_names]
fig, axes = plt.subplots(1, len(genes), figsize=(4.6 * len(genes), 3.8))
for k, (g, ax) in enumerate(zip(genes, np.atleast_1d(axes))):
    for layer, label, colour in schemes:
        centres, prof = area_profile(layer, genes)
        v = prof[:, k]
        ax.plot(centres, v / max(v.mean(), 1e-9), "o-", color=colour, label=label, ms=4)
    ax.set_title(g); ax.set_xlabel("cell area (µm²)")
np.atleast_1d(axes)[-1].legend(frameon=False, fontsize=8)
plt.tight_layout(); plt.show()
```

**Expected outcome.** The curves separate far more than they did for an abundant
housekeeping-type gene. `TRAC` and `CD52` mark lymphocytes, which are small and
have little cytoplasm, so they sit in the bottom area deciles. Under area
normalisation their expression is divided by a small number and goes up; under
median-scaling it goes up less.

**So for an immune-focused study**, the choice is not cosmetic — it changes how
strongly immune cells separate from everything else, and therefore how many immune
clusters you get. The defensible move is to report the immune fraction under both and
say so. The indefensible move is to try both silently and keep whichever gave more
clusters.

## 4.1 — Clustered, dispersed, or random?

```python
sq.gr.ripley(adata, cluster_key="cell_type", mode="L", n_simulations=200, seed=0)
res = adata.uns["cell_type_ripley_L"]["L_stat"]
t = "T cell"
sub = res[res["cell_type"] == t]
i = (sub["stats"] - sub["bins"]).abs().idxmax()
print(f"{t}: largest deviation at r = {sub.loc[i, 'bins']:.0f} um")
```

**Interpretation.** T cells in tumours are typically clustered, with the deviation
peaking somewhere in the 50–200 µm range — the size of perivascular cuffs and small
aggregates. If it peaks near 300–500 µm you are picking up tertiary lymphoid
structures. The radius is the finding; report it in microns.

---

## 4.2 — High Moran's I, low between-cluster variance

```python
sc.tl.rank_genes_groups(adata, "cell_type", method="wilcoxon")
means = pd.DataFrame(
    {ct: np.asarray(adata[adata.obs["cell_type"] == ct].X.mean(axis=0)).ravel()
     for ct in adata.obs["cell_type"].cat.categories}, index=adata.var_names)
spread = means.std(axis=1) / (means.mean(axis=1) + 1e-9)

mi = adata.uns["moranI"]
df = pd.DataFrame({"moran": mi["I"], "between_ct_cv": spread.reindex(mi.index)})
print(df.query("moran > 0.15 and between_ct_cv < 0.35").sort_values("moran", ascending=False).head(15))
```

**What produces this.** Genes responding to something continuous that cuts across
cell types: hypoxia (`VEGFA`, `CA9`) near necrosis, interferon response radiating
from an inflammatory focus, or a proliferation gradient. Also, less excitingly,
technical gradients — a focus or permeabilisation artefact will do this too, so check
the pattern is not a stripe.

These genes are precisely the ones a cluster-first analysis will miss.

---

## 4.3 — Choosing K

There is no statistically correct K. Defensible procedure:

1. Run K = 4, 6, 8, 10.
2. For each, print the composition heatmap and the spatial map.
3. Keep the largest K at which **every niche is still nameable** as a structure you
   would point to on a slide.
4. Check stability: re-run with three random seeds and confirm the niches persist.

```python
from sklearn.metrics import adjusted_rand_score
labs = {}
for seed in (0, 1, 2):
    labs[seed] = KMeans(n_clusters=6, n_init=10, random_state=seed).fit_predict(comp)
print("ARI 0 vs 1:", adjusted_rand_score(labs[0], labs[1]))
print("ARI 0 vs 2:", adjusted_rand_score(labs[0], labs[2]))
```

ARI below ~0.7 across seeds means your K is too high for the data to support.

---

## 4.4 — Same cell type, different niche

```python
CT = "Fibroblast / CAF"   # the most abundant fibroblast state in this section
sub = adata[adata.obs["cell_type"] == CT].copy()
counts = sub.obs["niche"].value_counts()
keep = counts[counts >= 100].index.tolist()
sub = sub[sub.obs["niche"].isin(keep)].copy()
sub.obs["niche"] = sub.obs["niche"].cat.remove_unused_categories()

sc.tl.rank_genes_groups(sub, "niche", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(sub, n_genes=5, standard_scale="var")
```

**Caveat to raise with the group.** Some of the difference will be spillover: a
fibroblast surrounded by tumour cells picks up tumour transcripts through
segmentation error, so `EPCAM` appearing as a "niche marker" for tumour-adjacent
fibroblasts is an artefact, not a state.

The check: is the difference in genes the *neighbours* express (spillover) or in
genes the neighbours do **not** express (real state change)? Only the second is
convincing. This distinction is worth ten minutes of discussion.

---

## 5.1 — Two distances at once

```python
from scipy.spatial import cKDTree
ref = {"tumour": "Tumour / epithelial", "vessel": "Endothelial"}
for name, ct in ref.items():
    m = (adata.obs["cell_type"] == ct).to_numpy()
    d, _ = cKDTree(adata.obsm["spatial"][m]).query(adata.obsm["spatial"], k=1)
    adata.obs[f"dist_{name}"] = d

sub = adata[adata.obs["cell_type"] == "T cell"].copy()
print(np.corrcoef(sub.obs["dist_tumour"], sub.obs["dist_vessel"])[0, 1])

import statsmodels.formula.api as smf
g = "TRAC"
sub.obs["expr"] = np.asarray(sub[:, g].X.todense()).ravel()
fit = smf.ols("expr ~ dist_tumour + dist_vessel", data=sub.obs).fit()
print(fit.summary().tables[1])
```

**The point.** Fit both distances together. Vessels sit in stroma, so
`dist_vessel` and `dist_tumour` are correlated, and a univariate association with one
is partly an association with the other. Reporting the marginal effect of a single
distance is the most common analytical error in this style of paper.

---

## 5.2 — Nuclear fraction by niche

```python
tx = pd.read_parquet(DATA / "transcripts_crop.parquet").query("qv >= 20")
niche = adata.obs["niche"].astype(str)
tx["niche"] = tx["cell_id"].astype(str).map(niche.to_dict())
tx = tx.dropna(subset=["niche"])

g = tx["feature_name"].value_counts().idxmax()
res = (tx[tx["feature_name"] == g]
       .groupby("niche")
       .agg(n=("overlaps_nucleus", "size"), nuclear=("overlaps_nucleus", "mean")))
print(g); print(res.round(3))
```

**Before believing it,** check nucleus size does not differ between niches — a larger
nucleus captures a larger share of transcripts by geometry alone. Regress
`nuclear_ratio` out, or restrict to cells in a narrow band of nucleus area.

---

## 5.3 — Grid correlation versus cell correlation

```python
BIN = 10
def raster(gene):
    m = (tx["feature_name"] == gene).to_numpy()
    img = np.zeros((ny, nx)); np.add.at(img, (gy[m], gx[m]), 1); return img.ravel()

g1, g2 = "EPCAM", "DCN"
grid_r = np.corrcoef(raster(g1), raster(g2))[0, 1]
a = np.asarray(adata[:, g1].X.todense()).ravel()
b = np.asarray(adata[:, g2].X.todense()).ravel()
cell_r = np.corrcoef(a, b)[0, 1]
print(f"grid r = {grid_r:.3f}   cell r = {cell_r:.3f}")
```

**How to read the discrepancy.**

- Grid **positive**, cell **negative** → the two genes are in *adjacent but distinct*
  cells. Correct and biologically meaningful: epithelium next to stroma.
- Grid ≈ cell, both positive → genuine co-expression in the same cells.
- Cell positive, grid ≈ 0 → suspicious. Likely spillover creating apparent
  co-expression that the raw molecule positions do not support.

The last case is the one to hunt for in your own data.

---

## Challenge 1 — The invasive front

Sketch of a defensible approach:

1. Define the tumour compartment from **niches**, not cell types, so it is already
   smoothed (notebook 04).
2. Compute a signed distance: negative inside the tumour compartment, positive
   outside. `scipy.ndimage.distance_transform_edt` on a rasterised mask is the
   cleanest way to get both signs.
3. Define the front as |signed distance| < w, and choose w from the data — the width
   over which composition changes, not a round number you liked.
4. Characterise: composition versus signed distance; then, within each cell type,
   expression versus signed distance.
5. Control for the obvious confounder: cells at the front differ in local density,
   which affects segmentation quality and therefore counts. Include `total_counts` as
   a covariate, or match on it.

The deliverable is a width in microns and a list of genes with a front-associated
gradient, per cell type.
