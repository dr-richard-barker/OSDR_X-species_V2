# Conserved mitochondrial suppression and organelle-specific transcriptomic responses to spaceflight across six eukaryotic species

**A human-anchored, cross-species meta-analysis of spaceflight transcriptomics from the NASA Open Science Data Repository (OSDR).**

> Cross-species meta-analysis · NASA OSDR · OrthoDB v12 · 22 datasets · 6 species

---

## Overview

Spaceflight perturbs gene expression across every organism studied, but how much
of that response is *conserved* at the level of subcellular compartments and
metabolic pathways has been unclear. This project performs a cross-species
meta-analysis of **22 transcriptomic datasets** from NASA OSDR spanning **six
eukaryotic species** — *Homo sapiens*, *Mus musculus*, *Drosophila
melanogaster*, *Caenorhabditis elegans*, *Saccharomyces cerevisiae*, and
*Arabidopsis thaliana* — and asks which organelle-level transcriptional
responses are shared.

Differentially expressed genes from each dataset are mapped onto a **human-anchored
orthology matrix built from OrthoDB v12**, assigned to **16 subcellular
compartments**, and combined with Fisher's combined probability test to detect
direction-consistent, cross-species responses.

## Key result

A **conserved downregulation of mitochondrial genes**: 49 orthologs show
significant, direction-consistent suppression across at least three species. The
signal is reinforced at the pathway level (TCA cycle, oxidative phosphorylation)
and in the organelle-enrichment and cofactor analyses.

| Metric | Value |
|--------|-------|
| Datasets | 22 (NASA OSDR) |
| Species | 6 |
| Differentially expressed genes | 20,333 |
| Orthology backbone | OrthoDB v12 — 18,030 genes (651 with orthologs in all six species) |
| Subcellular compartments tested | 16 |
| Conserved mito-suppressed orthologs (≥3 species) | 49 |

## Repository contents

```
.
├── figures/                     # Publication figures (PNG + SVG)
│   ├── fig1_study_overview       # Datasets / species composition
│   ├── fig2_volcano_plots        # Per-species differential expression
│   ├── fig3_orthology_upset      # Ortholog sharing across species
│   ├── fig4_organelle_enrichment_heatmap
│   ├── fig5_conserved_organelle_schematic
│   ├── fig6_pathway_hsa00020/00190/03010/03050   # TCA, OxPhos, ribosome, proteasome
│   └── fig7_cofactor_summary
├── zenodo_repo/
│   └── supplementary_tables/    # Machine-readable supplementary data
│       ├── Table_S3_GOCC_enrichment_all.csv
│       ├── Table_S4_GOCC_enrichment_classified.csv
│       ├── Table_S7_organelle_conservation_summary.csv
│       ├── Table_S8_organelle_conservation_final.csv
│       └── Table_S11_cofactor_gene_mapping.csv
├── manuscript/
│   └── npj_microgravity_manuscript.docx
└── README.md
```

## Data sources

- **NASA Open Science Data Repository (OSDR / GeneLab):** https://osdr.nasa.gov —
  cite each contributing OSD accession individually (see the manuscript's data
  table).
- **OrthoDB v12** (orthology backbone): https://www.orthodb.org

## Methods & reproduction

The full methodology — dataset selection, per-species differential expression,
orthology mapping, compartment assignment, and Fisher's combined test — is
described in the Methods section of `manuscript/npj_microgravity_manuscript.docx`.
The supplementary tables in `zenodo_repo/supplementary_tables/` are the
machine-readable outputs of those analyses.

## Code availability

The analysis scripts are in [`scripts/`](scripts/); the methodology they
implement is documented in the Methods section of the manuscript. Together they
regenerate the figures and supplementary tables from the NASA OSDR inputs and the
OrthoDB v12 orthology backbone.

## License

**CC BY 4.0** — Creative Commons Attribution 4.0 International (see `LICENSE`).
OSDR-derived data remain subject to the NASA Open Data policy; cite the original
OSD accessions.

## Citation

> **TODO:** add the author list (ORCID), venue, and DOI once available.

Barker, R. et al. (2026). *Conserved mitochondrial suppression and
organelle-specific transcriptomic responses to spaceflight across six eukaryotic
species.* Manuscript in preparation (targeted at *npj Microgravity*).
