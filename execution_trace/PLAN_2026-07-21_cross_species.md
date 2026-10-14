> **Historical record, not current results.** This is the plan the Biomni agent wrote on 21 Jul 2026
> before running this analysis. It is kept unchanged below for provenance. Some details changed when the
> analysis was run: for example, the plan names OrthoDB v11, but [`METHODS.md`](../METHODS.md) records OrthoDB v12.
> For what was actually done, see `METHODS.md` and the manuscript.

# Plan: Comparative Multi-Species Spaceflight Transcriptomics — Conserved Subcellular Pathway Response

## Summary
End-to-end analysis pipeline + Zenodo-ready repository + npj Microgravity manuscript draft. Pull processed RNA-seq count matrices from the NASA OSDR Biological Data API for 6 species (mouse, human, *Drosophila*, *C. elegans*, yeast, *Arabidopsis*), run per-dataset DESeq2 (spaceflight vs ground, tissue as covariate), map DEGs to Eukaryota-level OrthoDB orthologous groups anchored on human, perform GO Cellular Component (subcellular) enrichment per species, identify conserved organelle-specific responses via cross-species meta-analysis, visualize KEGG pathways with DEG log2FC and cofactor annotations via **ggkegg** (the installed `ggpathway` equivalent), and package everything as a reproducible repo + Word manuscript with figures.

## 1. Data Acquisition (NASA OSDR Biological Data API)
- Query `/v2/query/assays/` filtered by `study.characteristics.organism` × `study.factor value.spaceflight` × RNA-Seq assay technology to enumerate datasets per species.
- For each dataset (OSD-###), pull ISA metadata + GeneLab Reference Consensus Pipeline (RCP) processed **raw count matrices** (CSV) + sample metadata via the OSDR File API.
- Species: *Mus musculus*, *Homo sapiens*, *Drosophila melanogaster*, *Caenorhabditis elegans*, *Saccharomyces cerevisiae*, *Arabidopsis thaliana*.
- **Fallback**: if RNA-seq spaceflight datasets are insufficient for a species (e.g., *C. elegans* is largely microarray), include microarray datasets and process with limma (noted explicitly in methods).
- Multi-tissue: include all available tissues per species; tissue modeled as covariate in DESeq2 design.

## 2. Differential Expression (per dataset)
- DESeq2 per dataset: design `~ tissue + condition` (spaceflight vs ground control) where multi-tissue; `~ condition` where single-tissue.
- Filter: keep genes with ≥10 counts across ≥min(n_samples, 3) samples.
- DEG threshold: padj < 0.05 (BH), |log2FC| > 1.
- Output: per-dataset DEG tables (CSV) in `results/deg/`.
- Within-species meta-analysis: combine per-dataset p-values for each gene via Fisher's method (ortholog-aware), FDR-correct across the union.

## 3. Orthology Mapping (OrthoDB v11, Eukaryota level)
- Download OrthoDB v11 Eukaryota-level orthologous groups (OGs) covering all 6 species.
- Map each species' gene IDs (Ensembl/TAIR/SGD/FlyBase/WormBase) to OGs.
- Anchor on human gene symbols: each OG → representative human ortholog (if any) for downstream annotation.
- Build a unified orthology matrix: rows = OGs, columns = species, cells = gene ID(s).
- Genes without cross-kingdom orthologs retained for species-specific enrichment but excluded from the conserved-response meta-analysis.

## 4. Subcellular Location Enrichment
- GO Cellular Component (GOCC) ORA via clusterProfiler per species DEG list:
  - org.Hs.eg.db (human), org.Mm.eg.db (mouse), org.Dm.eg.db (fly), org.Ce.eg.db (worm), org.Sc.sgd.db (yeast), org.At.tair.db (Arabidopsis).
  - FDR < 0.05 (BH).
- Group GOCC terms into broad organelle categories (nucleus, mitochondrion, ER, Golgi, peroxisome, lysosome, cytoskeleton, ribosome, proteasome, plasma membrane, extracellular) using the GO hierarchy (SubcellulaRVis-style rollup).
- Cross-species conservation: an organelle category is "conserved" if enriched (FDR < 0.05) in the same direction (up/down DEG enrichment) in ≥4 of 6 species; supplemented by Fisher-combined p-value across species on the shared orthologous gene set.
- Output: summary table of enriched subcellular pathways (CSV) — organelle, species enriched, direction, FDR per species, combined FDR, representative GO terms.

## 5. Pathway & Cofactor Visualization (ggkegg)
- Use **ggkegg** (installed; the `ggpathway` equivalent) to render KEGG pathways for the top conserved enriched pathways (e.g., oxidative phosphorylation, TCA cycle, proteasome, ribosome, glutathione metabolism).
- Overlay per-species DEG log2FC onto pathway nodes (ortholog-mapped to KEGG orthology IDs).
- Map required cofactors (NAD+, FAD, CoA, Fe-S, PLP, thiamine-PP, biotin, lipoate, heme, Zn2+) to pathway enzymes via KEGG Compound/Reaction/Module annotations; annotate cofactor demand on the pathway diagrams.
- Consistent species color-blind-friendly palette (viridis discrete) across all figures; species→color mapping fixed globally.

## 6. Figures (all .svg + .png, viridis palette, Liberation Sans)
- Fig 1: Study/dataset overview (species × tissue × dataset count) — bar/tile plot.
- Fig 2: Per-species volcano plots (spaceflight vs ground) — 6-panel.
- Fig 3: Cross-species orthology overlap (UpSet of OGs containing DEGs).
- Fig 4: Subcellular compartment enrichment heatmap (organelle × species, signed −log10 FDR).
- Fig 5: Conserved organelle response schematic (cell diagram with enriched compartments colored by species) — via GenerateImage.
- Fig 6: ggkegg pathway diagram(s) for top conserved pathway with per-species log2FC overlay + cofactor annotations.
- Fig 7: Cofactor requirement summary (cofactor × pathway × species).
- Supplementary figures as needed.

## 7. Zenodo-Ready Repository Structure
```
spaceflight-multispecies-transcriptomics/
├── README.md
├── LICENSE (CC-BY 4.0 code, CC0 data)
├── CITATION.cff
├── requirements.txt / renv.lock
├── data/
│   ├── raw_counts/        # GeneLab RCP downloaded
│   ├── metadata/          # sample metadata per dataset
│   └── orthodb/           # OrthoDB OG mapping
├── analysis/
│   ├── 01_download_osdr.R
│   ├── 02_deseq2_per_dataset.R
│   ├── 03_orthology_mapping.py
│   ├── 04_subcellular_enrichment.R
│   ├── 05_meta_analysis_conservation.R
│   ├── 06_pathway_cofactor_viz.R
│   └── run_all.sh
├── results/
│   ├── deg/               # per-dataset DEG tables
│   ├── enrichment/        # GOCC enrichment tables
│   ├── meta/              # cross-species meta-analysis
│   └── figures/           # all figures svg+png
├── methods/
│   └── technical_notes.md # normalization, orthology, stats
├── supplementary/
│   ├── Supp_Table_S1_datasets.csv
│   ├── Supp_Table_S2_DEGs_per_species.csv
│   ├── Supp_Table_S3_orthology_matrix.csv
│   ├── Supp_Table_S4_subcellular_enrichment.csv
│   └── Supp_Table_S5_cofactor_mapping.csv
└── manuscript/
    ├── npj_microgravity_draft.docx
    ├── figures/           # manuscript-quality figures
    └── references.bib
```

## 8. npj Microgravity Manuscript Draft (Word .docx)
- Format: Word .docx (npj accepts Word; most editable for submission), via the docx-generation skill.
- Structure: Title page, Abstract (~150 words, unreferenced), Introduction (no heading), Results, Discussion, Methods, Data Availability, Author Contributions, Competing Interests, References, Figure Legends, Tables.
- Length: ~5000 words main text, ≤10 display items, ≤70 references.
- Figures embedded inline + separate figure files.
- Methods: full technical detail on OSDR access, DESeq2 design, OrthoDB mapping, GOCC enrichment, conservation meta-analysis, ggkegg visualization.
- Data Availability: point to Zenodo DOI (placeholder) + NASA OSDR accessions.
- AI use disclosure in Methods per npj policy.

## 9. Technical Notes (methods/technical_notes.md)
- **Normalization**: DESeq2 median-of-ratios on GeneLab RCP raw counts; per-dataset to avoid cross-study batch confounding; within-species meta-analysis via Fisher's method.
- **Orthology**: OrthoDB v11 Eukaryota-level OGs; human-anchored; many-to-many handled by 1:1 preferred + in-paralog inclusion; genes lacking cross-kingdom orthologs analyzed species-specifically.
- **Statistics**: BH-FDR < 0.05 throughout; |log2FC| > 1 for DEGs; conservation = enrichment in same direction in ≥4/6 species + Fisher-combined FDR < 0.05 on shared ortholog set.
- **Visualization**: viridis discrete palette, species colors fixed globally; ggkegg for KEGG pathway diagrams; svglite for editable SVG; Liberation Sans font.

## 10. Compute / Resource Estimate
- **Inputs**: ~6 species × ~3-8 datasets each = ~20-40 datasets; processed count matrices ~5-50 MB each → ~1-2 GB total download.
- **DESeq2**: minutes per dataset on 8 CPU / 32 GB sandbox; total ~30-60 min for all datasets.
- **OrthoDB**: ~1-2 GB download, mapping minutes.
- **Enrichment + visualization**: minutes.
- **Total runtime**: ~2-4 hours, sandbox-feasible.
- **Execution target**: single right-sized sandbox machine (8 CPU / 32 GB) via ManageMachine; no HPC needed (processed counts, not raw FASTQ).
- **Storage**: ~5 GB working disk for raw + intermediates; final deliverables to /mnt/results/.

## Assumptions (delegated decisions)
- Manuscript format: Word .docx (most editable for npj submission).
- Reference/anchor species for orthology: human (richest GO/COMPARTMENTS/KEGG annotation).
- Conservation threshold: ≥4/6 species same direction + Fisher-combined FDR < 0.05.
- DEG threshold: padj < 0.05, |log2FC| > 1.
- Cofactor source: KEGG Compound/Reaction/Module annotations.
- If RNA-seq insufficient for a species, include microarray (limma) with explicit methods note.
- Color palette: viridis discrete, species→color fixed globally.

## Execution Order
1. Provision 8 CPU / 32 GB sandbox machine.
2. Query OSDR API, enumerate + download datasets per species; save raw counts + metadata.
3. Run DESeq2 per dataset; output DEG tables.
4. Download OrthoDB v11 Eukaryota OGs; build orthology matrix; map DEGs to OGs (human-anchored).
5. GOCC enrichment per species; roll up to organelle categories; cross-species meta-analysis for conservation.
6. ggkegg pathway diagrams for top conserved pathways + cofactor mapping.
7. Generate all figures (svg+png, viridis, media-output-checked).
8. Assemble Zenodo repo (scripts, data refs, methods, supplementary tables, CITATION.cff, LICENSE).
9. Draft npj Microgravity manuscript (.docx) with embedded figures via docx-generation skill.
10. Final QA: verify repo completeness, figure quality, manuscript sections.
