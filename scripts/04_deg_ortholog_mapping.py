# ======================================================================
# Step 4: Map DEGs to Human Orthologs
# 
# Description: Map all DEGs to human orthologs and find cross-species overlaps
# 
# Inputs: results/deg/all_degs_combined.csv, results/orthology/orthology_matrix_wide.csv
# Outputs: results/deg/all_degs_with_orthologs.csv, results/meta/cross_species_degs.csv
# 
# Language: PYTHON
# Extracted from analysis notebook (worker-0.ipynb)
# See METHODS.md for full parameter details.
# ======================================================================

# --- Python code block 1 ---
import pandas as pd

# Load orthology matrix
ortho_wide = pd.read_csv('/mnt/shared-workspace/shared/orthodb/orthology_matrix_wide.csv')
ortho_long = pd.read_csv('/mnt/shared-workspace/shared/orthodb/unified_orthology_matrix.csv')

# Load DEGs
deg = pd.read_csv('/mnt/shared-workspace/shared/results/deg/all_degs_combined.csv', low_memory=False)
deg['gene_id'] = deg['gene_id'].astype(str)

# Map each species' DEGs to human orthologs
species_cols = {
    'Mus musculus': 'Mus musculus',
    'Homo sapiens': 'human_ensembl',  # human maps to itself
    'Drosophila melanogaster': 'Drosophila melanogaster',
    'Caenorhabditis elegans': 'Caenorhabditis elegans',
    'Saccharomyces cerevisiae': 'Saccharomyces cerevisiae',
    'Arabidopsis thaliana': 'Arabidopsis thaliana'
}

# Build reverse mapping: species_gene_id -> human_ensembl
deg_with_ortholog = deg.copy()
deg_with_ortholog['human_ortholog'] = None

for species, col in species_cols.items():
    if species == 'Homo sapiens':
        # Human genes map to themselves (strip version)
        mask = deg_with_ortholog['species'] == 'Homo sapiens'
        deg_with_ortholog.loc[mask, 'human_ortholog'] = deg_with_ortholog.loc[mask, 'gene_id'].str.split('.').str[0]
    else:
        # Map via orthology matrix
        sp_map = ortho_long[ortho_long['species'] == species].set_index('species_gene_id')['human_ensembl'].to_dict()
        mask = deg_with_ortholog['species'] == species
        deg_with_ortholog.loc[mask, 'human_ortholog'] = deg_with_ortholog.loc[mask, 'gene_id'].map(sp_map)

# Summary
print("=== DEG to human ortholog mapping ===")
for species in deg_with_ortholog['species'].unique():
    sub = deg_with_ortholog[deg_with_ortholog['species'] == species]
    n_total = len(sub)
    n_mapped = sub['human_ortholog'].notna().sum()
    n_deg_mapped = sub[sub['is_deg']]['human_ortholog'].notna().sum()
    n_deg = sub['is_deg'].sum()
    print(f"  {species}: {n_mapped}/{n_total} genes mapped ({n_mapped/n_total*100:.1f}%), {n_deg_mapped}/{n_deg} DEGs mapped ({n_deg_mapped/n_deg*100:.1f}%)")

# Save
deg_with_ortholog.to_csv('/mnt/shared-workspace/shared/results/deg/all_degs_with_orthologs.csv', index=False)
print(f"\nSaved {len(deg_with_ortholog)} rows with ortholog mappings")

# Cross-species DEG overlap (via human ortholog)
deg_only = deg_with_ortholog[deg_with_ortholog['is_deg'] & deg_with_ortholog['human_ortholog'].notna()].copy()
print(f"\nDEGs with human ortholog: {len(deg_only)}")
cross_species = deg_only.groupby('human_ortholog').agg(
    n_species=('species', 'nunique'),
    species_list=('species', lambda x: '|'.join(sorted(set(x)))),
    directions=('deg_direction', lambda x: '|'.join(x))
).reset_index()

# Genes that are DEGs in multiple species
multi_species = cross_species[cross_species['n_species'] >= 2]
print(f"Human orthologs that are DEGs in >=2 species: {len(multi_species)}")
print(f"Human orthologs that are DEGs in >=3 species: {len(cross_species[cross_species['n_species']>=3])}")
print(f"Human orthologs that are DEGs in >=4 species: {len(cross_species[cross_species['n_species']>=4])}")

# Conserved direction (same direction in all species where DEG)
def check_conserved_direction(directions):
    dirs = directions.split('|')
    return 'all_up' if all(d == 'up' for d in dirs) else ('all_down' if all(d == 'down' for d in dirs) else 'mixed')

multi_species['conservation'] = multi_species['directions'].apply(check_conserved_direction)
print(f"\nConserved direction in multi-species DEGs:")
print(multi_species['conservation'].value_counts())

multi_species.to_csv('/mnt/shared-workspace/shared/results/meta/cross_species_degs.csv', index=False)
print(f"\nSaved cross-species DEG table")

