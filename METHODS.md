# Methods

## Data acquisition

### Dataset selection

We queried the NASA Open Science Data Repository (OSDR) Biological Data API (https://visualization.osdr.nasa.gov/biodata/api/v2/query/) to enumerate RNA-seq and microarray datasets with spaceflight vs. ground control contrasts for six species: *Homo sapiens*, *Mus musculus*, *Drosophila melanogaster*, *Caenorhabditis elegans*, *Saccharomyces cerevisiae*, and *Arabidopsis thaliana*. The query syntax used was:

```
https://visualization.osdr.nasa.gov/biodata/api/v2/query/assays/?investigation.study%20assays.study%20assay%20technology%20type&=study.characteristics.organism&=study.factor%20value.spaceflight&format=csv
```

We identified RNA-seq datasets for mouse (78 studies), human (19), Arabidopsis (18), fly (6), and C. elegans (1, flight-only — no ground control). Yeast had 0 RNA-seq datasets. For species lacking suitable RNA-seq controls, we used microarray as a fallback: C. elegans (OSD-35, 42, 112, 113) and yeast (OSD-62).

The final selection comprised 22 datasets (7 mouse, 4 human, 3 Arabidopsis, 3 fly, 4 worm, 1 yeast). Processed data (DE tables, count matrices, sample tables) were downloaded via:

```
https://osdr.nasa.gov/geode-py/ws/studies/{OSD}/download?file={FILENAME}
```

### Tissue metadata

Tissue information was extracted from the `study.characteristics.material type` field (not `study.characteristics.tissue`). Multi-tissue datasets were retained with tissue modeled as a covariate in the design matrix.

## Differential expression analysis

### Pre-computed DE tables (GeneLab RCP)

For 19 of 22 datasets, we used GeneLab Research Community Portal (RCP) pre-computed differential expression tables, which use DESeq2 for RNA-seq and limma for microarray. These tables contain all pairwise contrasts with columns: Log2fc, T.stat, P.value, Adj.p.value.

### Contrast selection

Contrast names follow the pattern `Log2fc_(Space Flight & ...)v(Ground Control & ...)`. We selected 17 clean contrasts where the only difference between groups was spaceflight vs. ground control. Two design-mismatched contrasts (OSD-112, OSD-258) were retained as acceptable after manual review. Direction was normalized so that positive log2FC = upregulated in spaceflight (SF_vs_GC kept sign; GC_vs_SF flipped).

### De novo DE analysis

For OSD-96 (Drosophila RSEM counts, no pre-computed DE table), we ran DESeq2 with a design of `~ tissue + condition` (tissue as covariate), yielding 28 DEGs. For OSD-35 (C. elegans microarray, no replicates for padj), we ran limma with a nominal threshold of p < 0.05, |log2FC| > 0.5, yielding 207 DEGs.

### Tiered DEG thresholds

We applied tiered thresholds to accommodate dataset-specific statistical power:
- **Standard**: padj < 0.05, |log2FC| > 1 (19 datasets)
- **Relaxed**: padj < 0.10, |log2FC| > 0.5 (OSD-62 yeast, OSD-258 human)
- **Nominal**: p < 0.05, |log2FC| > 0.5 (OSD-35 worm, no replicates)

This yielded 20,333 DEGs across all datasets. Per-species counts: Arabidopsis 9,930; C. elegans 5,022; Fly 1,881; Human 1,366; Mouse 2,103; Yeast 31.

## Orthology mapping

### OrthoDB v12

We used OrthoDB v12 (current release, superseding v11) as the primary orthology resource. Bulk files were downloaded from https://data.orthodb.org/v12/download/odb_data_dump/:
- odb12v2_OGs.tab.gz (ortholog group definitions)
- odb12v2_OG2genes.tab.gz (gene-to-OG mapping)
- odb12v2_genes.tab.gz (gene metadata)
- odb12v2_species.tab.gz, odb12v2_levels.tab.gz

We filtered genes for the six target species (NCBI taxon IDs: 10090, 9606, 7227, 6239, 4932, 3702), yielding 167,551 genes. Eukaryota-level ortholog groups (taxon 2759) comprised 742,015 OGs. We built Ensembl-to-OrthoDB (132,119 mappings), UniProt (152,049), and symbol (167,551) cross-references.

### babelgene (human to mouse/fly/worm/yeast)

For human-to-mouse, fly, worm, and yeast orthology, we used the R package babelgene (v0.99.7), which aggregates orthology calls from multiple databases. We mapped 48,606 ortholog pairs from human to the four animal/fungal species using `orthologs(genes=..., species=..., human=TRUE, min_support=1, top=TRUE)`.

### Arabidopsis orthology

babelgene does not support Arabidopsis. We mapped Arabidopsis TAIR IDs (found in OrthoDB genes file column 6) to human orthologs via OrthoDB Eukaryota-level OGs, yielding 2,868 Arabidopsis-to-human orthologs. The low DEG mapping rate (3.2%) reflects the plant-specific gene complement.

### Unified orthology matrix

The final human-anchored orthology matrix contains 51,474 ortholog pairs covering 18,030 human genes, with 651 genes having orthologs in all 5 other species. DEG-to-ortholog mapping rates: Mouse 70.7%, Human 100%, Fly 28.6%, Worm 18.8%, Yeast 64.5%, Arabidopsis 3.2%.

## Subcellular location enrichment

We ran clusterProfiler's `enrichGO` (ontology = "CC", p-value cutoff = 0.05, Benjamini-Hochberg correction) per species for up- and down-regulated DEGs separately, using species-specific org.db packages (org.Hs.eg.db, org.Mm.eg.db, org.Dm.eg.db, org.Ce.eg.db, org.Sc.sgd.db, org.At.tair.db). This yielded 200 enriched GOCC terms total.

We manually classified all enriched terms into 16 organelle categories (Cell wall/ECM, Plasma membrane, Extracellular, Cytoskeleton, Nucleus, Proteasome, Ribosome, Mitochondrion, Membrane transport, Endoplasmic reticulum, Chloroplast/Plastid, Peroxisome, Lipid particle, Golgi apparatus, Other). Three refinement passes eliminated all "Other" classifications.

## Cross-species conservation meta-analysis

For each organelle category, we identified the set of orthologous genes that were tested for DE in each species. We applied Fisher's combined probability test on the per-species p-values for each orthologous gene, then adjusted across all 1,210 ortholog-organelle tests using Benjamini-Hochberg FDR.

A gene was considered "conserved" if:
- Fisher combined FDR < 0.05
- The gene was significant (DEG) in ≥3 species
- Direction was consistent (≥80% of significant species in the same direction)

This yielded 1,079 significant ortholog-organelle associations. Key result: **49 mitochondrial genes conserved downregulated across ≥3 species** (100% down, min padj = 7.7e-68).

## Pathway visualization and cofactor mapping

### KEGG pathway gene mapping

We fetched KEGG pathway-gene links via the KEGG REST API (https://rest.kegg.jp/link/hsa/{pathway_id}) for four conserved pathways:
- hsa00190 (Oxidative phosphorylation, 138 genes)
- hsa03050 (Proteasome, 46 genes)
- hsa03010 (Ribosome, 229 genes)
- hsa00020 (TCA cycle, 30 genes)

We mapped pathway genes (Entrez IDs) to human Ensembl IDs via org.Hs.eg.db, then joined to the cross-species DEG table via the orthology matrix to obtain per-species mean log2FC values. The final pathway-gene log2FC table covers 403 unique genes.

### Cofactor-gene mapping

We mapped 14 enzyme cofactors to human genes via the KEGG compound → reaction → enzyme → gene chain:
1. Compound → reaction: `rest.kegg.jp/link/rn/cpd:{C_ID}`
2. Reaction → enzyme: `rest.kegg.jp/link/ec/rn`
3. Enzyme → human gene: `rest.kegg.jp/link/hsa/ec`

This yielded 498 unique human genes across 12 cofactors (NAD+, NADH, NADP+, NADPH, FAD, FADH2, FMN, CoA, TPP, Lipoamide, Biotin, Heme). Fe-S cluster and Zn2+ returned 0 genes because KEGG tracks these as prosthetic groups rather than reaction substrates — a known limitation.

### ggkegg pathway diagrams

We used the ggkegg R package (v1.11.0) to fetch KEGG KGML files and render pathway diagrams with per-species log2FC overlays. Each pathway was rendered as a 6-panel figure (one panel per species), with gene nodes colored by mean log2FC (blue = downregulated, red = upregulated, grey = no data) using a diverging scale (limits -2 to 2, squished).

## Visualization

All figures use:
- **Color palette**: viridis-based, colorblind-friendly, with species colors fixed globally: Human #440154FF, Mouse #3B528BFF, Fly #21918CFF, Worm #5EC962FF, Yeast #FDE725FF, Arabidopsis #9B1946FF
- **Font**: Liberation Sans (metric-equivalent to Arial)
- **Format**: SVG (editable text) and PNG (300 dpi)
- **Diverging FC scale**: blue (#2166AC) → white → red (#B2182B), midpoint 0

## AI use disclosure

Data analysis scripting, figure generation, and manuscript drafting were assisted by Biomni (Phylo, https://phylo.ai). All statistical analyses, data interpretation, and scientific conclusions were conducted and verified by the authors. This AI assistance is disclosed in accordance with npj Microgravity editorial policy on AI use.
