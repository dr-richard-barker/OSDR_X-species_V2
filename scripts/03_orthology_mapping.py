# ======================================================================
# Step 3: Orthology Mapping (OrthoDB v12 + babelgene)
# 
# Description: Build human-anchored orthology matrix using OrthoDB v12 and babelgene
# 
# Inputs: OrthoDB v12 bulk files, babelgene R package
# Outputs: results/orthology/orthology_matrix_wide.csv, unified_orthology_matrix.csv
# 
# Language: PYTHON
# Extracted from analysis notebook (worker-0.ipynb)
# See METHODS.md for full parameter details.
# ======================================================================

# --- Python code block 1 ---
import pandas as pd
import os

ortho_dir = '/mnt/shared-workspace/shared/orthodb'
deg = pd.read_csv('/mnt/shared-workspace/shared/results/deg/all_degs_combined.csv', low_memory=False)

# Load mappings
ensembl_map = pd.read_csv(os.path.join(ortho_dir, 'ensembl_to_orthodb.tab'), sep='\t', header=None, names=['ensembl_id','orthodb_gene_id','org_id'])
uniprot_map = pd.read_csv(os.path.join(ortho_dir, 'uniprot_to_orthodb.tab'), sep='\t', header=None, names=['uniprot_id','orthodb_gene_id','org_id'])
symbol_map = pd.read_csv(os.path.join(ortho_dir, 'symbol_to_orthodb.tab'), sep='\t', header=None, names=['symbol','orthodb_gene_id','org_id'])

# Load OG2genes mapping
og2genes = pd.read_csv(os.path.join(ortho_dir, 'our_genes_OG_mapping.tab'), sep='\t', header=None, names=['og_id','orthodb_gene_id'])

# Load Eukaryota OGs
euk_ogs = pd.read_csv(os.path.join(ortho_dir, 'eukaryota_OGs.tab'), sep='\t', header=None, names=['og_id','level','description'])
euk_og_set = set(euk_ogs['og_id'])

# Filter OG2genes to Eukaryota-level OGs only
og2genes_euk = og2genes[og2genes['og_id'].isin(euk_og_set)]
print(f"OG2genes (Eukaryota level): {len(og2genes_euk)} mappings, {og2genes_euk['og_id'].nunique()} OGs")

# Build orthodb_gene_id -> og_id mapping
gene_to_og = dict(zip(og2genes_euk['orthodb_gene_id'], og2genes_euk['og_id']))

# Add org_id to mappings
ensembl_map['species_taxid'] = ensembl_map['org_id'].str.split('_').str[0]
uniprot_map['species_taxid'] = uniprot_map['org_id'].str.split('_').str[0]
symbol_map['species_taxid'] = symbol_map['org_id'].str.split('_').str[0]

# Map species names
taxid_to_species = {
    '10090': 'Mus musculus', '9606': 'Homo sapiens', '7227': 'Drosophila melanogaster',
    '6239': 'Caenorhabditis elegans', '4932': 'Saccharomyces cerevisiae', '3702': 'Arabidopsis thaliana'
}

# Check what gene ID types our DEGs use per species
print("\n=== DEG gene ID formats per species ===")
for sp in deg['species'].unique():
    sub = deg[deg['species']==sp]
    sample_ids = sub['gene_id'].dropna().unique()[:5]
    print(f"  {sp}: {sample_ids}")

# Map DEGs to OrthoDB gene IDs
# Strategy: 
# 1. Try Ensembl ID (strip version) for mouse/human/fly/worm
# 2. Try UniProt ID for yeast/Arabidopsis (via symbol bridge)
# 3. Try gene symbol as fallback

def map_gene_to_orthodb(gene_id, symbol, species):
    if pd.isna(gene_id):
        return None
    gene_id = str(gene_id).strip()
    
    # Try Ensembl (strip version number)
    ens_base = gene_id.split('.')[0]
    match = ensembl_map[ensembl_map['ensembl_id'] == ens_base]
    if len(match) > 0:
        return match.iloc[0]['orthodb_gene_id']
    
    # Try UniProt (via symbol)
    if pd.notna(symbol):
        sym = str(symbol).strip()
        match = uniprot_map[uniprot_map['uniprot_id'] == sym]
        if len(match) > 0:
            return match.iloc[0]['orthodb_gene_id']
    
    # Try symbol
    if pd.notna(symbol):
        sym = str(symbol).strip()
        match = symbol_map[(symbol_map['symbol'] == sym) & (symbol_map['species_taxid'] == str(species))]
        if len(match) > 0:
            return match.iloc[0]['orthodb_gene_id']
    
    return None

# Test mapping on a sample
print("\n=== Test mapping ===")
test_samples = deg.groupby('species').first().reset_index()[['species','gene_id','symbol']]
for _, row in test_samples.iterrows():
    taxid = [k for k,v in taxid_to_species.items() if v == row['species']][0]
    odb_id = map_gene_to_orthodb(row['gene_id'], row['symbol'], taxid)
    print(f"  {row['species']}: gene_id={row['gene_id']}, symbol={row['symbol']} -> orthodb={odb_id}")

# --- Python code block 2 ---
import pandas as pd
import urllib.request
import json
import time
import os

# Load Arabidopsis DEGs
deg = pd.read_csv('/mnt/shared-workspace/shared/results/deg/all_degs_combined.csv', low_memory=False)
ath_degs = deg[deg['species'] == 'Arabidopsis thaliana']
ath_gene_ids = ath_degs['gene_id'].unique()
print(f"Arabidopsis genes to map: {len(ath_gene_ids)}")

# Use OrthoDB API: search for Arabidopsis genes by name (TAIR ID)
# The API endpoint: https://data.orthodb.org/v12/search?query=AT1G01010&species=3702
# But this searches text, not gene IDs. Let me try the genesearch endpoint

# Actually, let me use a different approach: the OrthoDB /tab endpoint
# First, get all OGs at Eukaryota level that contain Arabidopsis genes
# Then map those OGs to human genes

# We already have the OG2genes mapping filtered for our species
# Let me use the genes file to find Arabidopsis gene IDs in OrthoDB
# The issue is that Arabidopsis TAIR IDs aren't in the genes file directly

# Let me try the OrthoDB API with UniProt IDs as a bridge
# Arabidopsis genes in our DEG data use TAIR IDs (AT1G01010)
# OrthoDB genes file has UniProt IDs for Arabidopsis

# Better approach: use the symbol mapping we already built
symbol_map = pd.read_csv('/mnt/shared-workspace/shared/orthodb/symbol_to_orthodb.tab', sep='\t', header=None, names=['symbol','orthodb_gene_id','org_id'])
ath_symbol_map = symbol_map[symbol_map['org_id'].str.startswith('3702')]
print(f"Arabidopsis symbol mappings in OrthoDB: {len(ath_symbol_map)}")

# Map Arabidopsis DEGs by symbol
ath_degs_with_symbol = ath_degs[ath_degs['symbol'].notna()].copy()
ath_degs_with_symbol['symbol_clean'] = ath_degs_with_symbol['symbol'].astype(str).str.split('|').str[0]

# Merge with symbol map
ath_mapped = ath_degs_with_symbol.merge(ath_symbol_map[['symbol','orthodb_gene_id']], 
                                        left_on='symbol_clean', right_on='symbol', how='left',
                                        suffixes=('','_ortho'))
n_mapped = ath_mapped['orthodb_gene_id'].notna().sum()
print(f"Arabidopsis DEGs mapped by symbol: {n_mapped}/{len(ath_degs_with_symbol)} ({n_mapped/len(ath_degs_with_symbol)*100:.1f}%)")

# Also try mapping by TAIR ID directly (some may be in the genes file as locus tags)
# The genes file column 4 has gene symbols, column 6 has Ensembl IDs
# For Arabidopsis, the "Ensembl" column might have TAIR-like IDs
genes_ath = pd.read_csv('/mnt/shared-workspace/shared/orthodb/our_species_genes.tab', sep='\t', header=None,
                        names=['orthodb_gene_id','org_id','protein_id','symbol','uniprot_id','ensembl_id','col7','description','location','scaffold','quality'])
genes_ath = genes_ath[genes_ath['org_id'].str.startswith('3702')]
print(f"\nArabidopsis genes in OrthoDB: {len(genes_ath)}")
print("Sample ensembl_id values:", genes_ath['ensembl_id'].dropna().unique()[:5])

# Check if TAIR IDs are in the protein_id or symbol columns
print("Sample protein_ids:", genes_ath['protein_id'].iloc[:5].tolist())
print("Sample symbols:", genes_ath['symbol'].dropna().unique()[:10])

# --- Python code block 3 ---
import pandas as pd
import os

ortho_dir = '/mnt/shared-workspace/shared/orthodb'

# Load Arabidopsis genes from OrthoDB with TAIR IDs
genes_ath = pd.read_csv(os.path.join(ortho_dir, 'our_species_genes.tab'), sep='\t', header=None,
                        names=['orthodb_gene_id','org_id','protein_id','symbol','uniprot_id','ensembl_id','col7','description','location','scaffold','quality'])
genes_ath = genes_ath[genes_ath['org_id'].str.startswith('3702')]
genes_ath_with_tair = genes_ath[genes_ath['ensembl_id'].notna() & (genes_ath['ensembl_id'] != '')]
print(f"Arabidopsis genes with TAIR IDs: {len(genes_ath_with_tair)}")

# Build TAIR -> orthodb_gene_id mapping
tair_to_orthodb = dict(zip(genes_ath_with_tair['ensembl_id'], genes_ath_with_tair['orthodb_gene_id']))

# Load OG2genes mapping
og2genes = pd.read_csv(os.path.join(ortho_dir, 'our_genes_OG_mapping.tab'), sep='\t', header=None, names=['og_id','orthodb_gene_id'])
euk_ogs = pd.read_csv(os.path.join(ortho_dir, 'eukaryota_OGs.tab'), sep='\t', header=None, names=['og_id','level','description'])
euk_og_set = set(euk_ogs['og_id'])
og2genes_euk = og2genes[og2genes['og_id'].isin(euk_og_set)]
gene_to_og = dict(zip(og2genes_euk['orthodb_gene_id'], og2genes_euk['og_id']))

# Map Arabidopsis DEGs to OGs via TAIR IDs
deg = pd.read_csv('/mnt/shared-workspace/shared/results/deg/all_degs_combined.csv', low_memory=False)
ath_degs = deg[deg['species'] == 'Arabidopsis thaliana'].copy()
ath_degs['orthodb_gene_id'] = ath_degs['gene_id'].map(tair_to_orthodb)
n_mapped = ath_degs['orthodb_gene_id'].notna().sum()
print(f"Arabidopsis DEGs mapped to OrthoDB: {n_mapped}/{len(ath_degs)} ({n_mapped/len(ath_degs)*100:.1f}%)")

# Map to OGs
ath_degs['og_id'] = ath_degs['orthodb_gene_id'].map(gene_to_og)
n_og = ath_degs['og_id'].notna().sum()
print(f"Arabidopsis DEGs mapped to OGs: {n_og}/{len(ath_degs)} ({n_og/len(ath_degs)*100:.1f}%)")

# Now find human genes in the same OGs
# Get all human genes from OrthoDB
genes_all = pd.read_csv(os.path.join(ortho_dir, 'our_species_genes.tab'), sep='\t', header=None,
                        names=['orthodb_gene_id','org_id','protein_id','symbol','uniprot_id','ensembl_id','col7','description','location','scaffold','quality'])
genes_human = genes_all[genes_all['org_id'].str.startswith('9606')]
human_ensembl_to_orthodb = dict(zip(genes_human['ensembl_id'], genes_human['orthodb_gene_id']))
human_orthodb_to_ensembl = dict(zip(genes_human['orthodb_gene_id'], genes_human['ensembl_id']))

# Build OG -> human gene mapping
og2genes_human = og2genes_euk.merge(genes_human[['orthodb_gene_id','ensembl_id']], on='orthodb_gene_id')
og_to_human = {}
for _, row in og2genes_human.iterrows():
    og = row['og_id']
    ens = row['ensembl_id']
    if pd.notna(ens) and ens != '':
        if og not in og_to_human:
            og_to_human[og] = []
        og_to_human[og].append(ens)

# Map Arabidopsis DEGs to human orthologs via OGs
ath_degs['human_ortholog'] = ath_degs['og_id'].map(lambda og: '|'.join(og_to_human.get(og, [])) if pd.notna(og) else None)
n_human = ath_degs['human_ortholog'].notna().sum()
print(f"Arabidopsis DEGs with human ortholog via OG: {n_human}/{len(ath_degs)} ({n_human/len(ath_degs)*100:.1f}%)")

# Save Arabidopsis orthology mapping
ath_ortho = ath_degs[ath_degs['human_ortholog'].notna()][['gene_id','symbol','og_id','human_ortholog']].drop_duplicates()
ath_ortho.to_csv(os.path.join(ortho_dir, 'arabidopsis_to_human_orthologs.csv'), index=False)
print(f"\nSaved {len(ath_ortho)} Arabidopsis-to-human ortholog mappings")
print("Sample:")
print(ath_ortho.head())

