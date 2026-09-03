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

## Objects, and what the dot does

`adata` is an **object**: one thing holding the counts, the per-cell table, the
per-gene table and the coordinates, so they cannot drift out of sync.

The dot reaches inside it.

```python
adata.obs        # the per-cell table
adata.var        # the per-gene table
adata.obs.head() # a function that lives inside the table
```

**Parentheses are the thing to watch.** Without them you are asking for data; with
them you are asking the object to do something.

```python
adata.obs.shape   # data about the table   -> (41235, 12)
adata.obs.head()  # run this function      -> the first five rows
adata.obs.head    # forgot the ()          -> a description of the function
```

Type `adata.` and press **Tab** to see everything inside it.

### Coming from R?

| Task | R / Seurat | Python / scanpy |
|---|---|---|
| per-cell metadata | `seurat@meta.data` | `adata.obs` |
| first rows | `head(df)` | `df.head()` |
| dimensions | `dim(obj)` | `obj.shape` |
| a column | `df$total_counts` | `df["total_counts"]` |
| number of cells | `ncol(seurat)` | `adata.n_obs` |

In R the dot is just a character in a name (`data.frame`); in Python it always means
"look inside". And a Seurat object is genes × cells, while `AnnData` is cells × genes —
so `adata.obs` has one row per cell. Python also counts from 0, R from 1.

## Repeating something: the `for` loop

Four parts, and the last one is the unfamiliar bit.

```python
for r in [0.3, 0.6, 1.0]:
    explore(resolution=r)
```

- `for` starts the loop
- `r` is a name you invent; it takes each value in turn
- `[0.3, 0.6, 1.0]` is the list to walk through
- the `:` ends the header, and the **indented** line below is the body

Indentation is the syntax, not decoration — it defines what is inside the loop.
Jupyter indents for you after a colon.

```python
for r in [0.3, 1.0]:
    explore(resolution=r)   # inside — runs twice
print("done")               # outside — runs once
```

In R you would write `for (r in c(0.3, 0.6)) { ... }`. Python drops the parentheses
and braces and uses a colon plus indentation instead; `c(...)` becomes `[...]`.

Forgetting the colon gives `SyntaxError`; inconsistent indentation gives
`IndentationError`. Both name the line.

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
