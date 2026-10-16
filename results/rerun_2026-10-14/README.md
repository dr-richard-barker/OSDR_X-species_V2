# Re-run with the conservation test as METHODS describes it (2026-10-14)

The `/abai check` of 13–14 Oct found that the code behind the "conserved" results did not do
what `METHODS.md` says. This folder holds the organelle classification and conservation
meta-analysis re-run after three fixes:

| Script | Before | Now |
|---|---|---|
| `02_deg_extraction.py` | OSD-258 silently switched to nominal p < 0.01 when < 10 genes passed, still labelled "relaxed" | relaxed threshold only (padj < 0.10, \|log2FC\| > 0.5), as METHODS states |
| `05_gocc_enrichment.R` (final pass) | "er membrane" etc. matched as substrings, so "inner membrane" and "outer membrane" were filed as ER | the ER abbreviations must start a word |
| `06_conservation_meta_analysis.R` | see below | the METHODS rule |

`06` previously pooled every gene×dataset row (paralogs and several datasets per species) as
independent Fisher inputs; counted species that merely had a p-value toward "≥ 3 species";
called a gene "down" if every DEG call was down, even a single one; dropped p = 1 rows; and
applied BH within each organelle. It now uses one p-value per species (the smallest p among
that species' rows, Šidák-corrected: 1 − (1 − p_min)^k), Fisher's test over those, BH across
all ortholog–organelle tests, and calls a pair **conserved** if Fisher FDR < 0.05, the gene is
a DEG in ≥ 3 species, and ≥ 80% of those species share one direction (a species whose DEG
rows disagree counts as agreeing with neither).

**Nothing else has been regenerated.** `supplementary_tables/`, the figures, the manuscript,
the compiled PDF, the site, the slides and the infographic still describe the earlier analysis.

## Inputs

All tracked in this repository:

- `supplementary_tables/Table_S2_all_DEGs_with_orthologs.csv.gz` (492,431 rows), with the
  OSD-258 rows re-thresholded, because the raw DE tables `02` reads are not tracked:
  `is_deg = padj < 0.10 & |log2fc| > 0.5`, and `deg_direction` set from the sign of log2fc.
  OSD-258 goes from 182 DEGs to 1, and the total from 20,333 to 20,152.
- `supplementary_tables/Table_S3_GOCC_enrichment_all.csv` (the enrichGO results, with gene IDs).
- `results/orthology/unified_orthology_matrix.csv`.

`05` was run from its final pass (block 7), which re-classifies the enrichGO terms in S3;
enrichGO itself was not re-run. R 4.6.1, data.table 1.18.6.1.

**Control:** the original final pass of `05` and the original `06`, run the same way on the
original S2, reproduce S4, S5, S6, S7 and S8 exactly.

## Organelle classification (`gocc_enrichment_classified.csv`, `organelle_enrichment_summary.csv`)

9 of the 200 enriched terms change:

| Term | Species, direction | Before | Now |
|---|---|---|---|
| mitochondrial inner membrane | fly, down | Mitochondrion \| ER | Mitochondrion |
| organelle inner membrane | fly, down | ER | Other |
| organelle inner membrane | Arabidopsis, up | ER | Other |
| outer membrane | Arabidopsis, up | ER | Other |
| organelle outer membrane | Arabidopsis, up | ER | Other |
| plastid inner membrane | Arabidopsis, up | ER \| Chloroplast | Chloroplast |
| chloroplast inner membrane | Arabidopsis, up | ER \| Chloroplast | Chloroplast |
| plastid outer membrane | Arabidopsis, up | ER \| Chloroplast | Chloroplast |
| chloroplast outer membrane | Arabidopsis, up | ER \| Chloroplast | Chloroplast |

ER keeps one genuine term ("ER body"). In the enrichment summary that drives Fig 4, the ER
rows for Arabidopsis up (7 terms, signed score +298.3) and fly down (2 terms, −4.2) are gone,
and "Other" rows appear for Arabidopsis up (3 terms, +248.3) and fly down (1 term, −2.0).
There are 15 categories including "Other" (S4 had 14; METHODS says 16).

## Conservation (`organelle_conservation_*.csv`)

| | Ortholog–organelle tests | Fisher FDR < 0.05 | Conserved |
|---|---|---|---|
| Original (S6–S8) | 1,210 | 1,079 | 815 (any direction), of which 49 mitochondrial "down" |
| This re-run | 1,154 | 1,061 | **4** |

- The 4 conserved pairs are **one gene**, TUBB8B (ENSG00000173213, a β-tubulin), counted in
  four organelle categories (Cell wall/ECM, Cytoskeleton, Nucleus, Plasma membrane). It is a
  down-regulated DEG in 3 of the 5 species it was tested in; Fisher FDR 1.6e-31.
- No ortholog is a DEG in more than 3 species. Over all 1,154 tests the number of DEG species
  is 0 for 172, 1 for 808, 2 for 157 and 3 for 17.
- **Mitochondrion: 65 tests, 61 with Fisher FDR < 0.05, 1 ortholog that is a DEG in ≥ 3
  species, 0 conserved.**
- Fisher FDR alone filters almost nothing: 1,061 of the 1,154 tests pass it, against 17 that
  are a DEG in ≥ 3 species. 172 tests have no DEG species at all, and 108 of those still pass
  Fisher FDR. They are co-orthologs: the organelle gene lists are mapped to human genes
  through the orthology matrix, where one species gene can map to several human genes, while
  S2 assigns each DEG row to a single human ortholog.

`organelle_conservation_fisher.csv` has every test with its per-species counts
(`n_species_tested`, `n_deg_species`, `n_up_species`, `n_down_species`, `n_mixed_species`,
`direction_agreement`), `organelle_conservation_summary.csv` the per-organelle totals, and
`organelle_conservation_final.csv` only the conserved pairs.
