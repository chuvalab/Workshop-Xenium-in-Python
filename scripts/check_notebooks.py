#!/usr/bin/env python
"""Static sanity check on the workshop notebooks.

    python scripts/check_notebooks.py

Catches the errors that only show up when a student runs a notebook top to
bottom in a fresh kernel:

  * syntax errors in any code cell
  * names used before they are assigned (the classic "NameError: name 'x' is
    not defined" caused by editing or deleting an earlier cell)
  * cells that read a file the manifest does not list

It does NOT execute anything, so it needs neither the data nor a GPU and runs
in under a second. Run it after every notebook edit.

Note the limits of static analysis: it assumes cells run in order, and it
cannot see names created dynamically (exec, globals()[...]). A clean report
means "no obvious breakage", not "verified correct" — there is no substitute
for running the notebooks against real data before the workshop.
"""

from __future__ import annotations

import ast
import builtins
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = sorted((ROOT / "notebooks").glob("*.ipynb"))

# Names that exist in a notebook kernel without being assigned in a cell.
KERNEL_GLOBALS = {
    "get_ipython", "display", "In", "Out", "_", "__", "___", "exit", "quit",
}


def collect_bindings(tree: ast.AST) -> set[str]:
    """Every name this cell makes available to later cells (approximate)."""
    out: set[str] = set()
    for node in ast.walk(tree):
        # Function/lambda parameters must be collected unconditionally: putting
        # this in the elif chain below means ast.FunctionDef matches the earlier
        # branch and its arguments are never seen, which reports every default
        # parameter as an undefined name.
        if isinstance(node, (ast.Lambda, ast.FunctionDef, ast.AsyncFunctionDef)):
            args = node.args
            for a in (*getattr(args, "posonlyargs", []), *args.args, *args.kwonlyargs):
                out.add(a.arg)
            if args.vararg:
                out.add(args.vararg.arg)
            if args.kwarg:
                out.add(args.kwarg.arg)

        if isinstance(node, ast.Name) and isinstance(node.ctx, (ast.Store, ast.Del)):
            out.add(node.id)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            out.add(node.name)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                out.add((alias.asname or alias.name).split(".")[0])
        elif isinstance(node, ast.ExceptHandler) and node.name:
            out.add(node.name)
        elif isinstance(node, ast.Global):
            out.update(node.names)
        elif isinstance(node, ast.comprehension):
            for t in ast.walk(node.target):
                if isinstance(t, ast.Name):
                    out.add(t.id)
    return out


def loaded_names(tree: ast.AST) -> set[str]:
    return {
        n.id for n in ast.walk(tree)
        if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)
    }


def clean_source(cell: dict) -> str:
    """Strip IPython magics, which are not valid Python."""
    lines = []
    for line in "".join(cell["source"]).split("\n"):
        stripped = line.lstrip()
        if stripped.startswith(("%", "!", "?")):
            continue
        lines.append(line)
    return "\n".join(lines)


def check_notebook(path: Path) -> list[str]:
    nb = json.loads(path.read_text(encoding="utf-8"))
    problems: list[str] = []

    defined = set(dir(builtins)) | KERNEL_GLOBALS
    for i, cell in enumerate(nb.get("cells", [])):
        if not isinstance(cell, dict):
            problems.append(f"cell {i}: malformed (not an object)")
            continue
        if cell.get("cell_type") != "code":
            continue

        src = clean_source(cell)
        try:
            tree = ast.parse(src)
        except SyntaxError as exc:
            problems.append(f"cell {i}: SyntaxError line {exc.lineno}: {exc.msg}")
            continue

        bound = collect_bindings(tree)
        missing = loaded_names(tree) - defined - bound
        if missing:
            first = src.strip().split("\n")[0][:50]
            problems.append(
                f"cell {i}: uses {sorted(missing)} before assignment  | {first}"
            )
        defined |= bound

    return problems


def check_data_references() -> list[str]:
    """Every DATA / "file" referenced in a notebook should be in the manifest."""
    manifest = ROOT / "data" / "data_manifest.yml"
    if not manifest.exists():
        return []
    listed = set(re.findall(r"- name:\s*(\S+)", manifest.read_text()))
    # files the notebooks create themselves, so they need no manifest entry
    generated = {
        "ovarian_qc.h5ad", "ovarian_niches.h5ad", "ovarian_annotated_reference.h5ad",
    }
    problems = []
    for path in NOTEBOOKS:
        text = path.read_text(encoding="utf-8")
        for ref in set(re.findall(r'DATA / \\"([^"\\]+)\\"', text)):
            if ref not in listed and ref not in generated:
                problems.append(f"{path.name}: reads {ref!r}, not in data_manifest.yml")
    return problems


def main() -> int:
    if not NOTEBOOKS:
        print(f"no notebooks found under {ROOT / 'notebooks'}")
        return 1

    total = 0
    for path in NOTEBOOKS:
        problems = check_notebook(path)
        total += len(problems)
        status = "OK" if not problems else f"{len(problems)} problem(s)"
        print(f"{path.name:<38} {status}")
        for p in problems:
            print(f"    {p}")

    data_problems = check_data_references()
    if data_problems:
        print("\nData references")
        for p in data_problems:
            print(f"    {p}")
        total += len(data_problems)

    print()
    if total:
        print(f"{total} problem(s) found.")
        return 1
    print("All notebooks look consistent when run top to bottom.")
    print("This is a static check — still run them against real data before the workshop.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
