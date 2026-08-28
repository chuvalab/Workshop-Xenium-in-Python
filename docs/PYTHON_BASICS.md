# Python survival guide

For anyone who has not used Python much. You do **not** need to read this to follow
the workshop — every cell runs as written. Read it when you want to change something
and are not sure how.

---

## Running cells

- **Shift+Enter** runs the cell and moves to the next one.
- **Ctrl+Enter** runs it and stays put — useful when editing the same cell repeatedly.
- The `[ ]` to the left shows `[*]` while running and a number when finished.
- Cells share memory and **must be run in order**. If something is "not defined", you
  probably skipped a cell above.
- Broke something? *Kernel → Restart Kernel and Run All Cells* starts clean.

---

## The five things you will actually change

### 1. A number

```python
MIN_COUNTS = 10        # change to 25, re-run
```

### 2. A piece of text (a "string")

Quotes matter, and so does capitalisation. `"EPCAM"` is not `"epcam"`.

```python
GENE = "EPCAM"         # change to "PTPRC"
```

### 3. A list

Square brackets, comma-separated:

```python
MY_GENES = ["EPCAM", "PTPRC", "COL1A1"]
```

### 4. True / False

Capitalised, no quotes:

```python
RUN_PEARSON = True
```

### 5. A colour

Either a name (`"red"`) or a hex code (`"#001158"`).

---

## Reading the objects you will meet

`adata` holds everything. Four parts matter:

| What | Where | Think of it as |
|---|---|---|
| the counts | `adata.X` | a big table, cells × genes |
| per-cell info | `adata.obs` | a spreadsheet, one row per cell |
| per-gene info | `adata.var` | a spreadsheet, one row per gene |
| coordinates | `adata.obsm["spatial"]` | two columns, x and y in microns |

Useful one-liners:

```python
adata                          # a summary of everything
adata.obs.head()               # first five cells
adata.obs.columns              # what can I colour a plot by?
adata.obs["total_counts"].describe()    # min, max, median of one column
list(adata.var_names[:20])     # first twenty gene names
adata.n_obs, adata.n_vars      # how many cells, how many genes
```

Is a gene on the panel?

```python
"EPCAM" in adata.var_names     # True or False
```

Find genes whose name contains something:

```python
[g for g in adata.var_names if "COL" in g]
```

---

## The three algorithms, in one line each

| | What it does |
|---|---|
| **PCA** | squashes 5,000 gene measurements per cell into ~50 numbers, keeping what separates cells |
| **neighbour graph** | links each cell to the cells most similar to it — the object everything else uses |
| **UMAP** | draws that graph on a page. A picture, nothing is calculated from it |
| **Leiden** | cuts the graph into groups. It never looks at the UMAP |

The one thing to remember: **UMAP is for looking, Leiden is for grouping, and they
both read the same graph.** Distances and blob sizes on a UMAP do not mean anything.

---

## Errors, and what they mean

| Error | Usually means |
|---|---|
| `NameError: name 'x' is not defined` | you skipped a cell above — run them in order |
| `KeyError: 'EPCAM'` | that name is not in the table; check spelling and capitals |
| `SyntaxError` | a missing quote, bracket or comma, usually on the line above |
| `IndentationError` | spacing changed; undo with Ctrl+Z |
| `ValueError: ... shape ...` | two things have different lengths; usually a cell was re-run out of order |
| The kernel dies | out of memory — close other apps and *Restart Kernel* |

The single most useful habit: **read the last line of the error first.** Everything
above it is the path Python took to get there; the last line is what went wrong.

---

## If you get stuck

Put the red sticky note on your laptop lid. Then keep going — every notebook has
checkpoints where you can load a pre-computed file and rejoin, so falling behind on
one cell never costs you the rest of the session.
