# Supplementary Tables

S1 and S2 are stored gzipped (about 15 and 16 MB; 57 and 61 MB uncompressed). Run `gunzip -k` on them
to get the plain CSVs.

| Table | File | Description |
|-------|------|-------------|
| S1 | Table_S1_all_DEGs_combined.csv.gz | All genes tested across the 22 datasets (492,431 rows), with 20,333 flagged as DEGs in `is_deg` (gene_id, symbol, log2fc, pvalue, padj, is_deg, deg_direction, osd, species, threshold) |
| S2 | Table_S2_all_DEGs_with_orthologs.csv.gz | DEGs with mapped human orthologs (492,431 rows including non-DEG genes with FC values) |
| S3 | Table_S3_GOCC_enrichment_all.csv | All 200 enriched GOCC terms across species and directions |
| S4 | Table_S4_GOCC_enrichment_classified.csv | GOCC terms classified into organelle categories (14 categories occur in this table) |
| S5 | Table_S5_organelle_enrichment_summary.csv | Per-organelle, per-species enrichment summary with signed scores |
| S6 | Table_S6_organelle_conservation_fisher.csv | Fisher's combined probability test results (1,210 ortholog-organelle tests) |
| S7 | Table_S7_organelle_conservation_summary.csv | Organelle-level conservation summary (n_orthologs, n_sig, n_conserved, n_up, n_down) |
| S8 | Table_S8_organelle_conservation_final.csv | Final conserved ortholog-organelle associations (FDR < 0.05) |
| S9 | Table_S9_cross_species_DEGs.csv | DEGs significant in ≥2 species with ortholog mapping |
| S10 | Table_S10_conserved_organelles.csv | Organelles conserved across ≥3 species with direction |
| S11 | Table_S11_cofactor_gene_mapping.csv | 498 human genes mapped to 12 enzyme cofactors via KEGG |
| S12 | Table_S12_pathway_gene_log2fc.csv | Per-species mean log2FC for 403 genes across 4 KEGG pathways |
| S13 | Table_S13_cofactor_pathway_summary.csv | Cofactor × pathway × species mean FC summary |
| S14 | Table_S14_pathway_gene_log2fc_wide.csv | Pathway gene log2FC in wide format (gene × species) |
| S15 | Table_S15_pathway_DEGs_only.csv | Pathway genes that are DEGs (110 rows, 97 unique genes) |
