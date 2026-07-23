# LaTeX manuscript (npj Microgravity / Springer Nature style)

This folder assembles the manuscript in the official **Springer Nature LaTeX
template** (`sn-jnl` document class), which is the format npj Microgravity accepts
and typesets from.

```
latex/
├── main.tex          # the assembled manuscript (sn-jnl class, sn-nature refs)
├── references.bib    # 20 references, transcribed from the manuscript
├── figures/          # the 8 figures used in the text (PNG)
└── README.md         # this file
```

## How to compile

The `sn-jnl.cls` class and `sn-nature.bst` bibliography style are **not vendored
here** — they ship with Springer Nature's official template. Two ways to build:

### Option A — Overleaf (recommended, one click)
1. On [Overleaf](https://www.overleaf.com), create a new project from the
   **"Springer Nature Article Template (sn-jnl)"** (Templates → search
   "Springer Nature").
2. Replace that project's `main.tex` with this `main.tex`, and upload
   `references.bib` and the `figures/` folder.
3. Set the compiler to **pdfLaTeX** and compile. The template already contains
   `sn-jnl.cls` and the `.bst` files.

### Option B — local TeX Live / MiKTeX
1. Download the Springer Nature LaTeX template from
   springernature.com (author support → LaTeX) and place `sn-jnl.cls` and
   `sn-nature.bst` in this folder (or your TeX tree).
2. Build:
   ```bash
   pdflatex main
   bibtex   main
   pdflatex main
   pdflatex main
   ```

## Status / TODO before submission

- [ ] **Not yet compile-tested** — authored without a local TeX install; build
      once on Overleaf and fix any stragglers.
- [ ] **Author block** (`main.tex`) — confirm full author list, ORCID(s), and
      affiliation (currently a single-author placeholder).
- [ ] **References** — transcribed verbatim from the manuscript; several lack
      DOIs and some page ranges look like placeholders (e.g. "1--16"). Verify
      each entry and add DOIs (marked `% TODO` in `references.bib`).
- [ ] **Figures** — currently the repo PNGs. For final submission npj prefers
      vector (PDF/EPS) or ≥300 dpi; swap the files in `figures/` and the paths
      resolve unchanged.
- [ ] **Supplementary tables S1–S15** — referenced in-text; the machine-readable
      CSVs live in `../../zenodo_repo/supplementary_tables/`. Bundle per journal
      instructions at submission.

## Source

Ported from `../npj_microgravity_manuscript.docx`. Body text, figure legends,
Table 1, and all references are the author's own content — nothing was invented.
