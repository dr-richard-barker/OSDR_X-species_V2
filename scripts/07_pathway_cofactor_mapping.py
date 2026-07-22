# ======================================================================
# Step 7: KEGG Pathway and Cofactor Mapping
# 
# Description: Fetch KEGG pathway-gene and cofactor-gene links via REST API; build log2FC tables
# 
# Inputs: KEGG REST API (rest.kegg.jp), results/deg/all_degs_with_orthologs.csv
# Outputs: results/meta/cofactor_gene_mapping.csv, results/meta/pathway_gene_log2fc.csv, results/kegg/*.csv
# 
# Language: PYTHON
# Extracted from analysis notebook (worker-0.ipynb)
# See METHODS.md for full parameter details.
# ======================================================================

# --- Python code block 1 ---
import pandas as pd
import os

CACHE = "/workspace/kegg_cache"
RESULTS = "/mnt/shared-workspace/shared/results/meta"
KEGG_OUT = "/mnt/shared-workspace/shared/results/kegg"
os.makedirs(KEGG_OUT, exist_ok=True)

# Cofactor metadata
COFACTORS = {
    "C00003": ("NAD+", "Oxidoreductase (hydride transfer)"),
    "C00004": ("NADH", "Oxidoreductase (hydride transfer)"),
    "C00006": ("NADP+", "Anabolic oxidoreductase"),
    "C00005": ("NADPH", "Anabolic oxidoreductase"),
    "C00016": ("FAD", "Oxidoreductase (2e- transfer)"),
    "C00018": ("FADH2", "Oxidoreductase (2e- transfer)"),
    "C00010": ("CoA", "Acyl-group transfer"),
    "C20124": ("Fe-S cluster", "1e- transfer, redox sensing"),
    "C00061": ("FMN", "Oxidoreductase (1-2e- transfer)"),
    "C00068": ("TPP (thiamine PP)", "Decarboxylation, C-C cleavage"),
    "C00120": ("Biotin", "CO2 fixation (carboxylation)"),
    "C00252": ("Lipoamide", "Acyl transfer, reductive acylation"),
    "C00032": ("Heme", "O2 transport, e- transfer, peroxidation"),
    "C00038": ("Zn2+", "Hydrolase, structural, catalytic"),
}

# Pathways
PATHWAYS = {
    "hsa00190": "Oxidative phosphorylation",
    "hsa03050": "Proteasome",
    "hsa03010": "Ribosome",
    "hsa00020": "TCA cycle",
}

# 1. Build cofactor -> reaction -> ec -> hsa gene mapping
# Load reaction->ec and ec->hsa
rn_ec = pd.read_csv(f"{CACHE}/all_rn_ec.txt", sep="\t", header=None, names=["rn","ec"])
ec_hsa = pd.read_csv(f"{CACHE}/all_ec_hsa.txt", sep="\t", header=None, names=["ec","hsa_gene"])
print(f"rn->ec: {len(rn_ec)}, ec->hsa: {len(ec_hsa)}")

cofactor_gene_rows = []
for cid, (name, function) in COFACTORS.items():
    fn = f"{CACHE}/cof_{cid}_rn.txt"
    if not os.path.exists(fn):
        continue
    cof_rn = pd.read_csv(fn, sep="\t", header=None, names=["cpd","rn"])
    # cofactor -> reaction -> ec
    merged = cof_rn.merge(rn_ec, on="rn", how="inner")
    # ec -> hsa gene
    merged = merged.merge(ec_hsa, on="ec", how="inner")
    merged["cofactor_id"] = cid
    merged["cofactor_name"] = name
    merged["cofactor_function"] = function
    cofactor_gene_rows.append(merged[["cofactor_id","cofactor_name","cofactor_function","rn","ec","hsa_gene"]])
    print(f"{name} ({cid}): {len(cof_rn)} rxns -> {merged['hsa_gene'].nunique()} human genes")

cofactor_gene = pd.concat(cofactor_gene_rows, ignore_index=True).drop_duplicates()
# Clean hsa gene id (strip 'hsa:')
cofactor_gene["entrez"] = cofactor_gene["hsa_gene"].str.replace("hsa:","",regex=False)
print(f"\nTotal cofactor-gene mappings: {len(cofactor_gene)} rows, {cofactor_gene['entrez'].nunique()} unique genes")
cofactor_gene.to_csv(f"{RESULTS}/cofactor_gene_mapping.csv", index=False)
print(f"Saved: {RESULTS}/cofactor_gene_mapping.csv")
print(cofactor_gene.head())

# --- Python code block 2 ---
import pandas as pd

CACHE = "/workspace/kegg_cache"
RESULTS = "/mnt/shared-workspace/shared/results/meta"
KEGG_OUT = "/mnt/shared-workspace/shared/results/kegg"

PATHWAYS = {
    "hsa00190": "Oxidative phosphorylation",
    "hsa03050": "Proteasome",
    "hsa03010": "Ribosome",
    "hsa00020": "TCA cycle",
}

# 1. Load pathway genes
pathway_gene_rows = []
for pid, pname in PATHWAYS.items():
    fn = f"{CACHE}/{pid}_genes.txt"
    df = pd.read_csv(fn, sep="\t", header=None, names=["pathway_ref","hsa_gene"])
    df["pathway_id"] = pid
    df["pathway"] = pname
    df["entrez"] = df["hsa_gene"].str.replace("hsa:","",regex=False)
    pathway_gene_rows.append(df[["pathway_id","pathway","entrez"]])
    print(f"{pname} ({pid}): {len(df)} genes")

pathway_genes = pd.concat(pathway_gene_rows, ignore_index=True).drop_duplicates()
print(f"\nTotal pathway-gene rows: {len(pathway_genes)}, {pathway_genes['entrez'].nunique()} unique genes")

# 2. Load cross-species DEG data with orthologs
print("\nLoading all_degs_with_orthologs.csv (61MB)...")
degs = pd.read_csv("/mnt/shared-workspace/shared/results/deg/all_degs_with_orthologs.csv")
print(f"DEG table: {degs.shape}")
print(f"Columns: {list(degs.columns)}")
print(degs.head(3))

# --- Python code block 3 ---
import pandas as pd

# Load orthology matrix wide to get human ensembl -> entrez/symbol mapping
ortho_wide = pd.read_csv("/mnt/shared-workspace/shared/orthodb/orthology_matrix_wide.csv")
print(f"Orthology matrix wide: {ortho_wide.shape}")
print(f"Columns: {list(ortho_wide.columns)}")
print(ortho_wide.head(3))

# --- Python code block 4 ---
import pandas as pd
import numpy as np

CACHE = "/workspace/kegg_cache"
RESULTS = "/mnt/shared-workspace/shared/results/meta"
KEGG_OUT = "/mnt/shared-workspace/shared/results/kegg"

PATHWAYS = {
    "hsa00190": "Oxidative phosphorylation",
    "hsa03050": "Proteasome",
    "hsa03010": "Ribosome",
    "hsa00020": "TCA cycle",
}

# Load pathway genes
pathway_gene_rows = []
for pid, pname in PATHWAYS.items():
    df = pd.read_csv(f"{CACHE}/{pid}_genes.txt", sep="\t", header=None, names=["pathway_ref","hsa_gene"])
    df["pathway_id"] = pid
    df["pathway"] = pname
    df["entrez"] = df["hsa_gene"].str.replace("hsa:","",regex=False)
    pathway_gene_rows.append(df[["pathway_id","pathway","entrez"]])
pathway_genes = pd.concat(pathway_gene_rows, ignore_index=True).drop_duplicates()

# Load human ensembl->entrez->symbol
ens_map = pd.read_csv(f"{CACHE}/human_ensembl_entrez_symbol.csv")
ens_map["ENTREZID"] = ens_map["ENTREZID"].astype(str)
pathway_genes["entrez"] = pathway_genes["entrez"].astype(str)

# Join pathway genes -> human ensembl
pathway_genes = pathway_genes.merge(ens_map, left_on="entrez", right_on="ENTREZID", how="left")
pathway_genes = pathway_genes.rename(columns={"ENSEMBL":"human_ensembl","SYMBOL":"SYMBOL"})
print(f"Pathway genes with human Ensembl: {pathway_genes['human_ensembl'].notna().sum()}/{len(pathway_genes)}")
print(pathway_genes.head())

# Load DEG data with orthologs
degs = pd.read_csv("/mnt/shared-workspace/shared/results/deg/all_degs_with_orthologs.csv")
# Keep only DEGs to reduce size, but we also want non-DEG log2FC for heatmap
# Actually for pathway viz we want all genes with FC, but let's keep DEGs + significant
degs_sig = degs[degs["is_deg"]==True].copy()
print(f"\nDEG rows (is_deg=True): {len(degs_sig)}")

# Join pathway genes to DEG data via human_ortholog
# degs has 'human_ortholog' column = human ensembl
pathway_fc = pathway_genes.merge(
    degs_sig[["gene_id","symbol","log2fc","pvalue","padj","osd","species","human_ortholog"]],
    left_on="human_ensembl", right_on="human_ortholog", how="left"
)
print(f"\nPathway-DEG join: {len(pathway_fc)} rows, {pathway_fc['log2fc'].notna().sum()} with FC data")
print(f"Pathway genes with DEG in >=1 species: {pathway_fc[pathway_fc['log2fc'].notna()]['human_ensembl'].nunique()}")

# Aggregate per species: mean log2FC across datasets for that species
agg = pathway_fc.dropna(subset=["log2fc"]).groupby(
    ["pathway_id","pathway","human_ensembl","SYMBOL","entrez","species"], as_index=False
).agg(
    mean_log2fc=("log2fc","mean"),
    n_datasets=("osd","nunique"),
    min_padj=("padj","min"),
)
print(f"\nAggregated pathway FC: {len(agg)} rows")
print(agg.head())

# Save
agg.to_csv(f"{RESULTS}/pathway_gene_log2fc.csv", index=False)
print(f"\nSaved: {RESULTS}/pathway_gene_log2fc.csv")

# Also make a wide version: gene x species
wide = agg.pivot_table(index=["pathway_id","pathway","human_ensembl","SYMBOL","entrez"],
                       columns="species", values="mean_log2fc", aggfunc="first").reset_index()
wide.to_csv(f"{KEGG_OUT}/pathway_gene_log2fc_wide.csv", index=False)
print(f"Saved wide: {KEGG_OUT}/pathway_gene_log2fc_wide.csv")
print(f"\nPathway gene coverage (DEGs only):")
print(wide.groupby("pathway")[["Mus musculus","Homo sapiens","Drosophila melanogaster",
    "Caenorhabditis elegans","Saccharomyces cerevisiae","Arabidopsis thaliana"]].apply(lambda x: x.notna().sum()))

# --- Python code block 5 ---
import pandas as pd
import numpy as np

CACHE = "/workspace/kegg_cache"
RESULTS = "/mnt/shared-workspace/shared/results/meta"
KEGG_OUT = "/mnt/shared-workspace/shared/results/kegg"

PATHWAYS = {
    "hsa00190": "Oxidative phosphorylation",
    "hsa03050": "Proteasome",
    "hsa03010": "Ribosome",
    "hsa00020": "TCA cycle",
}

# Rebuild pathway genes with human ensembl
pathway_gene_rows = []
for pid, pname in PATHWAYS.items():
    df = pd.read_csv(f"{CACHE}/{pid}_genes.txt", sep="\t", header=None, names=["pathway_ref","hsa_gene"])
    df["pathway_id"] = pid
    df["pathway"] = pname
    df["entrez"] = df["hsa_gene"].str.replace("hsa:","",regex=False)
    pathway_gene_rows.append(df[["pathway_id","pathway","entrez"]])
pathway_genes = pd.concat(pathway_gene_rows, ignore_index=True).drop_duplicates()

ens_map = pd.read_csv(f"{CACHE}/human_ensembl_entrez_symbol.csv")
ens_map["ENTREZID"] = ens_map["ENTREZID"].astype(str)
pathway_genes["entrez"] = pathway_genes["entrez"].astype(str)
pathway_genes = pathway_genes.merge(ens_map, left_on="entrez", right_on="ENTREZID", how="left")
pathway_genes = pathway_genes.rename(columns={"ENSEMBL":"human_ensembl","SYMBOL":"SYMBOL"})

# Load ALL DEG table rows (not just is_deg) to get log2FC for all genes
degs = pd.read_csv("/mnt/shared-workspace/shared/results/deg/all_degs_with_orthologs.csv")
print(f"Full DEG table: {degs.shape}")

# Join pathway genes to ALL gene FC data via human_ortholog
pathway_fc = pathway_genes.merge(
    degs[["gene_id","symbol","log2fc","pvalue","padj","is_deg","osd","species","human_ortholog"]],
    left_on="human_ensembl", right_on="human_ortholog", how="left"
)
print(f"Pathway-FC join (all genes): {len(pathway_fc)} rows, {pathway_fc['log2fc'].notna().sum()} with FC data")

# Aggregate per species: mean log2FC across datasets
agg = pathway_fc.dropna(subset=["log2fc"]).groupby(
    ["pathway_id","pathway","human_ensembl","SYMBOL","entrez","species"], as_index=False
).agg(
    mean_log2fc=("log2fc","mean"),
    n_datasets=("osd","nunique"),
    n_deg=("is_deg","sum"),
    min_padj=("padj","min"),
)
print(f"\nAggregated pathway FC (all genes): {len(agg)} rows, {agg['human_ensembl'].nunique()} unique genes")

# Save long format
agg.to_csv(f"{RESULTS}/pathway_gene_log2fc.csv", index=False)
print(f"Saved: {RESULTS}/pathway_gene_log2fc.csv")

# Wide format: gene x species
wide = agg.pivot_table(index=["pathway_id","pathway","human_ensembl","SYMBOL","entrez"],
                       columns="species", values="mean_log2fc", aggfunc="first").reset_index()
wide.to_csv(f"{KEGG_OUT}/pathway_gene_log2fc_wide.csv", index=False)
print(f"Saved wide: {KEGG_OUT}/pathway_gene_log2fc_wide.csv")

# Coverage stats
species_cols = ["Mus musculus","Homo sapiens","Drosophila melanogaster",
    "Caenorhabditis elegans","Saccharomyces cerevisiae","Arabidopsis thaliana"]
print(f"\nPathway gene coverage (all genes with FC):")
print(wide.groupby("pathway")[species_cols].apply(lambda x: x.notna().sum()))

# DEG-only subset for highlighting
deg_only = agg[agg["n_deg"]>0]
print(f"\nDEG-only pathway genes: {len(deg_only)} rows, {deg_only['human_ensembl'].nunique()} unique genes")
deg_only.to_csv(f"{KEGG_OUT}/pathway_degs_only.csv", index=False)
print(f"Saved: {KEGG_OUT}/pathway_degs_only.csv")

# --- Python code block 6 ---
import pandas as pd

RESULTS = "/mnt/shared-workspace/shared/results/meta"
KEGG_OUT = "/mnt/shared-workspace/shared/results/kegg"

# Load cofactor-gene mapping
cofactor_gene = pd.read_csv(f"{RESULTS}/cofactor_gene_mapping.csv")
cofactor_gene["entrez"] = cofactor_gene["entrez"].astype(str)

# Load pathway FC long
pathway_fc = pd.read_csv(f"{RESULTS}/pathway_gene_log2fc.csv")
pathway_fc["entrez"] = pathway_fc["entrez"].astype(str)

# Join cofactors to pathway FC
cofactor_fc = cofactor_gene[["cofactor_id","cofactor_name","cofactor_function","entrez"]].drop_duplicates().merge(
    pathway_fc, on="entrez", how="left"
)
print(f"Cofactor-pathway-FC: {len(cofactor_fc)} rows")
print(f"Cofactor genes with FC data: {cofactor_fc.dropna(subset=['mean_log2fc'])['entrez'].nunique()}")

# Summary: cofactor x pathway x species mean FC
cofactor_summary = cofactor_fc.dropna(subset=["mean_log2fc"]).groupby(
    ["cofactor_name","pathway","species"], as_index=False
).agg(
    n_genes=("entrez","nunique"),
    mean_fc=("mean_log2fc","mean"),
    n_down=("mean_log2fc", lambda x: (x<0).sum()),
    n_up=("mean_log2fc", lambda x: (x>0).sum()),
)
print(f"\nCofactor-pathway-species summary: {len(cofactor_summary)} rows")
print(cofactor_summary.head(10))

cofactor_summary.to_csv(f"{KEGG_OUT}/cofactor_pathway_summary.csv", index=False)
print(f"\nSaved: {KEGG_OUT}/cofactor_pathway_summary.csv")

# Also save cofactor-gene-level FC for the figure
cofactor_fc.to_csv(f"{KEGG_OUT}/cofactor_gene_log2fc.csv", index=False)
print(f"Saved: {KEGG_OUT}/cofactor_gene_log2fc.csv")

# Quick stats: which cofactors are most affected in OxPhos?
oxphos_cof = cofactor_summary[cofactor_summary["pathway"]=="Oxidative phosphorylation"]
print(f"\n=== Cofactor usage in OxPhos (mean FC across genes) ===")
print(oxphos_cof.pivot_table(index="cofactor_name", columns="species", values="mean_fc"))

