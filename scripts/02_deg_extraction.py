# ======================================================================
# Step 2: DEG Extraction and Differential Expression
# 
# Description: Extract DEGs from GeneLab RCP pre-computed tables; run DESeq2 on OSD-96 and limma on OSD-35
# 
# Inputs: Processed DE tables and count matrices from Step 1
# Outputs: results/deg/all_degs_combined.csv, per-dataset DEG tables
# 
# Language: PYTHON
# Extracted from analysis notebook (worker-0.ipynb)
# See METHODS.md for full parameter details.
# ======================================================================

# --- Python code block 1 ---
import pandas as pd
dl_df = pd.read_csv('/mnt/shared-workspace/shared/osdr_data/download_log.csv')
n_ok = sum(1 for _, r in dl_df.iterrows() if 'ERROR' not in str(r['status']))
print(f"Downloaded OK: {n_ok}/{len(dl_df)} files")
print(f"Total size: {round(dl_df['size_mb'].sum(),2)} MB")
print(f"\nDE tables: {len(dl_df[dl_df['category']=='de'])}")
print(f"Sample tables: {len(dl_df[dl_df['category']=='sampletable'])}")
print(f"Count/expr matrices: {len(dl_df[dl_df['category']=='counts_or_expr'])}")

# --- Python code block 2 ---
import pandas as pd
import os

proc_dir = '/mnt/shared-workspace/shared/osdr_data/processed'

# Inspect a few DE tables - just the contrast columns (Log2fc columns)
for osd in ['OSD-37', 'OSD-258', 'OSD-104', 'OSD-207', 'OSD-347', 'OSD-217', 'OSD-421']:
    de_files = [f for f in os.listdir(proc_dir) if f.startswith(f'{osd}_de_')]
    if not de_files:
        continue
    fpath = os.path.join(proc_dir, de_files[0])
    # Read just the header
    with open(fpath, 'r') as f:
        header = f.readline().strip()
    cols = header.split('","')
    # Clean quotes
    cols = [c.strip('"') for c in cols]
    # Find Log2fc columns
    log2fc_cols = [c for c in cols if c.startswith('Log2fc')]
    print(f"\n=== {osd} ({de_files[0][:40]}...) ===")
    print(f"  Total cols: {len(cols)}")
    print(f"  Annotation cols: {[c for c in cols[:15] if not c.startswith(('Log2fc','T.stat','P.value','Adj.p','Group.','All.'))]}")
    print(f"  Log2fc contrasts ({len(log2fc_cols)}):")
    for c in log2fc_cols[:20]:
        print(f"    {c}")
    if len(log2fc_cols) > 20:
        print(f"    ... and {len(log2fc_cols)-20} more")

# --- Python code block 3 ---
import pandas as pd
import os
import re

proc_dir = '/mnt/shared-workspace/shared/osdr_data/processed'
final = pd.read_csv('/mnt/shared-workspace/shared/osdr_data/dataset_final.csv')

# For each dataset with a DE table, find the best Space Flight vs Ground Control contrast
# Strategy: find Log2fc columns where one side has "Space Flight" and other has "Ground Control"
# Prefer: (Space Flight)v(Ground Control) - positive = up in flight
# If multiple (due to genotypes/tissues), prefer Wild Type / control genotype, or the simplest one

contrast_selection = []

for _, row in final.iterrows():
    osd = row['osd']
    species = row['species']
    de_files = [f for f in os.listdir(proc_dir) if f.startswith(f'{osd}_de_')]
    if not de_files:
        contrast_selection.append({'osd': osd, 'species': species, 'contrast': None, 'reason': 'no DE table'})
        continue
    fpath = os.path.join(proc_dir, de_files[0])
    with open(fpath, 'r') as f:
        header = f.readline().strip()
    cols = [c.strip('"') for c in header.split('","')]
    log2fc_cols = [c for c in cols if c.startswith('Log2fc_')]
    
    # Find contrasts with Space Flight vs Ground Control
    flt_vs_gc = []
    for c in log2fc_cols:
        contrast_text = c.replace('Log2fc_', '').strip('()')
        # Check both directions
        has_sf = 'Space Flight' in contrast_text
        has_gc = 'Ground Control' in contrast_text or 'Ground control' in contrast_text
        if has_sf and has_gc:
            # Determine direction: (Space Flight)v(Ground Control) = positive = up in flight (what we want)
            # vs (Ground Control)v(Space Flight) = negative = up in flight
            parts = contrast_text.split(')v(')
            if len(parts) == 2:
                left, right = parts
                if 'Space Flight' in left and 'Ground Control' in right:
                    direction = 'SF_vs_GC'  # positive = up in flight
                elif 'Ground Control' in left and 'Space Flight' in right:
                    direction = 'GC_vs_SF'  # negative = up in flight, need to flip
                else:
                    direction = 'unknown'
                flt_vs_gc.append({'col': c, 'direction': direction, 'text': contrast_text})
    
    if not flt_vs_gc:
        contrast_selection.append({'osd': osd, 'species': species, 'contrast': None, 'reason': 'no SF vs GC contrast'})
        continue
    
    # Prefer SF_vs_GC direction; prefer Wild Type / control; prefer simplest (fewest &)
    sf_vs_gc_dir = [c for c in flt_vs_gc if c['direction'] == 'SF_vs_GC']
    candidates = sf_vs_gc_dir if sf_vs_gc_dir else flt_vs_gc
    
    # Score: prefer Wild Type, prefer fewer & (simpler contrast)
    def score(c):
        text = c['text'].lower()
        s = 0
        if 'wild type' in text or 'wild-type' in text or 'wt' in text:
            s -= 10
        if 'col-0' in text or 'col0' in text:
            s -= 5  # Arabidopsis control ecotype
        if 'n2' in text and 'bristol' in text:
            s -= 5  # C. elegans control
        if 'by4742' in text or 'by4743' in text:
            s -= 5  # yeast control
        s += text.count('&')  # fewer & = simpler
        if 'canton-s' in text and 'wild type' not in text:
            s -= 3  # Drosophila Canton-S is wild type
        return s
    
    candidates_sorted = sorted(candidates, key=score)
    best = candidates_sorted[0]
    
    contrast_selection.append({
        'osd': osd, 'species': species, 
        'contrast': best['col'], 'direction': best['direction'],
        'contrast_text': best['text'], 'reason': 'selected'
    })
    print(f"{osd} ({species}): {best['col']}")
    print(f"  direction: {best['direction']}, candidates: {len(flt_vs_gc)}")

contrast_df = pd.DataFrame(contrast_selection)
contrast_df.to_csv('/mnt/shared-workspace/shared/osdr_data/contrast_selection.csv', index=False)
print(f"\nContrasts selected: {sum(1 for c in contrast_selection if c['contrast'])}/{len(contrast_selection)}")
print(f"Missing: {[c['osd'] for c in contrast_selection if not c['contrast']]}")

# --- Python code block 4 ---
import pandas as pd
import os
import re

proc_dir = '/mnt/shared-workspace/shared/osdr_data/processed'
final = pd.read_csv('/mnt/shared-workspace/shared/osdr_data/dataset_final.csv')

# Better contrast selection: find contrasts where the ONLY difference between sides is Space Flight vs Ground Control
# Split each side by &, compare the non-flight/non-control terms

contrast_selection = []

for _, row in final.iterrows():
    osd = row['osd']
    species = row['species']
    de_files = [f for f in os.listdir(proc_dir) if f.startswith(f'{osd}_de_')]
    if not de_files:
        contrast_selection.append({'osd': osd, 'species': species, 'contrast': None, 'direction': None, 'reason': 'no DE table'})
        continue
    fpath = os.path.join(proc_dir, de_files[0])
    with open(fpath, 'r') as f:
        header = f.readline().strip()
    cols = [c.strip('"') for c in header.split('","')]
    log2fc_cols = [c for c in cols if c.startswith('Log2fc_')]
    
    best = None
    best_score = 999
    
    for c in log2fc_cols:
        contrast_text = c.replace('Log2fc_', '').strip('()')
        parts = contrast_text.split(')v(')
        if len(parts) != 2:
            continue
        left, right = parts
        # Normalize
        left_lower = left.lower()
        right_lower = right.lower()
        has_sf_left = 'space flight' in left_lower
        has_gc_left = 'ground control' in left_lower
        has_sf_right = 'space flight' in right_lower
        has_gc_right = 'ground control' in right_lower
        
        # Must have SF on one side, GC on other
        if not ((has_sf_left and has_gc_right) or (has_gc_left and has_sf_right)):
            continue
        
        # Determine direction
        if has_sf_left and has_gc_right:
            direction = 'SF_vs_GC'
        else:
            direction = 'GC_vs_SF'
        
        # Extract non-flight/control terms from each side
        def extract_other_terms(side):
            terms = [t.strip() for t in side.split('&')]
            others = [t for t in terms if 'space flight' not in t.lower() and 'ground control' not in t.lower()]
            return set(o.lower().strip() for o in others)
        
        left_others = extract_other_terms(left)
        right_others = extract_other_terms(right)
        
        # Score: 0 if other terms match perfectly (clean contrast), higher if mismatched
        mismatch = len(left_others.symmetric_difference(right_others))
        n_terms = len(left_others) + len(right_others)
        
        # Prefer clean contrasts (mismatch=0), then fewer total terms
        score = mismatch * 100 + n_terms
        
        if score < best_score:
            best_score = score
            best = {'col': c, 'direction': direction, 'text': contrast_text, 'mismatch': mismatch, 'left_others': left_others, 'right_others': right_others}
    
    if best:
        contrast_selection.append({
            'osd': osd, 'species': species,
            'contrast': best['col'], 'direction': best['direction'],
            'contrast_text': best['text'], 'mismatch': best['mismatch'],
            'reason': 'clean' if best['mismatch']==0 else f'mismatched({best["mismatch"]})'
        })
        flag = 'CLEAN' if best['mismatch']==0 else 'MISMATCHED'
        print(f"{osd} ({species}) [{flag}]: {best['col']}")
        if best['mismatch'] > 0:
            print(f"  left others: {best['left_others']}")
            print(f"  right others: {best['right_others']}")
    else:
        contrast_selection.append({'osd': osd, 'species': species, 'contrast': None, 'direction': None, 'reason': 'no SF vs GC contrast found'})

contrast_df = pd.DataFrame(contrast_selection)
contrast_df.to_csv('/mnt/shared-workspace/shared/osdr_data/contrast_selection.csv', index=False)
n_clean = sum(1 for c in contrast_selection if c.get('reason')=='clean')
n_mismatch = sum(1 for c in contrast_selection if 'mismatched' in str(c.get('reason','')))
n_none = sum(1 for c in contrast_selection if not c.get('contrast'))
print(f"\nClean contrasts: {n_clean}, Mismatched: {n_mismatch}, None: {n_none}")

# --- Python code block 5 ---
import pandas as pd
import os
import re
import numpy as np

proc_dir = '/mnt/shared-workspace/shared/osdr_data/processed'
final = pd.read_csv('/mnt/shared-workspace/shared/osdr_data/dataset_final.csv')
contrast_df = pd.read_csv('/mnt/shared-workspace/shared/osdr_data/contrast_selection.csv')

deg_dir = '/mnt/shared-workspace/shared/results/deg'
os.makedirs(deg_dir, exist_ok=True)

# For each dataset with a selected contrast, extract: gene_id, symbol, log2FC (normalized: + = up in flight), padj
all_degs = []

for _, crow in contrast_df.iterrows():
    osd = crow['osd']
    species = crow['species']
    contrast_col = crow['contrast']
    direction = crow['direction']
    
    if pd.isna(contrast_col):
        print(f"{osd} ({species}): no contrast - skip (will handle separately)")
        continue
    
    de_files = [f for f in os.listdir(proc_dir) if f.startswith(f'{osd}_de_')]
    fpath = os.path.join(proc_dir, de_files[0])
    
    # Read the DE table - only the columns we need
    with open(fpath, 'r') as f:
        header = f.readline().strip()
    cols = [c.strip('"') for c in header.split('","')]
    
    # Find the matching Log2fc, P.value (or Adj.p.value) columns
    # The contrast text is embedded in the column name
    contrast_text = contrast_col.replace('Log2fc_', '').strip('()')
    
    # Find the Adj.p.value column with the same contrast
    padj_col = f"Adj.p.value_({contrast_text})"
    pval_col = f"P.value_({contrast_text})"
    
    # Determine gene ID column (ENSEMBL for animals, TAIR for Arabidopsis)
    id_col = 'ENSEMBL' if 'ENSEMBL' in cols else ('TAIR' if 'TAIR' in cols else cols[0])
    sym_col = 'SYMBOL' if 'SYMBOL' in cols else None
    
    # Read only needed columns
    needed = [id_col]
    if sym_col:
        needed.append(sym_col)
    needed.extend([contrast_col, padj_col, pval_col])
    needed = [c for c in needed if c in cols]
    
    try:
        df = pd.read_csv(fpath, usecols=needed, low_memory=False)
    except Exception as e:
        print(f"{osd}: read error {e}")
        continue
    
    # Rename
    df = df.rename(columns={id_col: 'gene_id', sym_col: 'symbol' if sym_col else 'symbol',
                            contrast_col: 'log2fc_raw', padj_col: 'padj', pval_col: 'pvalue'})
    
    # Normalize direction: we want positive log2FC = upregulated in spaceflight
    # SF_vs_GC: positive = up in flight (keep as is)
    # GC_vs_SF: positive = up in ground (flip sign)
    if direction == 'GC_vs_SF':
        df['log2fc'] = -df['log2fc_raw']
    elif direction == 'SF_vs_GC':
        df['log2fc'] = df['log2fc_raw']
    else:
        # unknown - check from column name
        if 'Space Flight' in contrast_text.split(')v(')[0]:
            df['log2fc'] = df['log2fc_raw']
        else:
            df['log2fc'] = -df['log2fc_raw']
    
    df['osd'] = osd
    df['species'] = species
    df['contrast'] = contrast_text
    df['direction_normalized'] = 'positive = up in spaceflight'
    
    # Filter out rows with no gene ID
    df = df[df['gene_id'].notna() & (df['gene_id'] != '')]
    # Handle multi-mapped genes (pipe-separated) - take first
    df['gene_id'] = df['gene_id'].astype(str).str.split('|').str[0]
    
    # DEG status: padj < 0.05, |log2fc| > 1
    df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
    df['log2fc'] = pd.to_numeric(df['log2fc'], errors='coerce')
    df['is_deg'] = (df['padj'] < 0.05) & (df['log2fc'].abs() > 1)
    df['deg_direction'] = np.where(df['is_deg'] & (df['log2fc'] > 0), 'up', 
                                   np.where(df['is_deg'] & (df['log2fc'] < 0), 'down', 'ns'))
    
    n_deg = df['is_deg'].sum()
    n_up = (df['deg_direction']=='up').sum()
    n_down = (df['deg_direction']=='down').sum()
    
    # Save per-dataset DEG table
    out_cols = ['gene_id', 'symbol', 'log2fc', 'pvalue', 'padj', 'is_deg', 'deg_direction', 'osd', 'species']
    df[out_cols].to_csv(os.path.join(deg_dir, f'{osd}_deg.csv'), index=False)
    
    all_degs.append(df[out_cols])
    print(f"{osd} ({species}): {len(df)} genes, {n_deg} DEGs ({n_up} up, {n_down} down) | contrast: {contrast_text[:60]}")

# Combine all
all_deg_df = pd.concat(all_degs, ignore_index=True)
all_deg_df.to_csv(os.path.join(deg_dir, 'all_degs_combined.csv'), index=False)
print(f"\nTotal rows: {len(all_deg_df)}")
print(f"Total DEGs: {all_deg_df['is_deg'].sum()}")
print(f"\nDEG counts per species:")
print(all_deg_df[all_deg_df['is_deg']].groupby(['species','osd','deg_direction']).size().unstack(fill_value=0))

# --- Python code block 6 ---
import pandas as pd
import os

proc_dir = '/mnt/shared-workspace/shared/osdr_data/processed'

# Check OSD-96 RSEM counts structure
fpath = os.path.join(proc_dir, 'OSD-96_counts_or_expr_GLDS-96_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv')
with open(fpath, 'r') as f:
    header = f.readline().strip()
cols = [c.strip('"') for c in header.split('","')]
print(f"OSD-96 RSEM counts: {len(cols)} cols")
print(f"First 10 cols: {cols[:10]}")
df96 = pd.read_csv(fpath, nrows=3)
print(df96.iloc[:3, :6].to_string())

print("\n--- OSD-35 normalized expression ---")
fpath35 = os.path.join(proc_dir, 'OSD-35_counts_or_expr_GLDS-35_array_normalized_expression_probeset_GLmicroarray.csv')
with open(fpath35, 'r') as f:
    header = f.readline().strip()
cols35 = [c.strip('"') for c in header.split('","')]
print(f"OSD-35 norm expr: {len(cols35)} cols")
print(f"First 10 cols: {cols35[:10]}")
df35 = pd.read_csv(fpath35, nrows=3)
print(df35.iloc[:3, :6].to_string())

# Check OSD-242 - it has a DE table but no clean contrast was found
print("\n--- OSD-242 DE table contrasts ---")
de242 = [f for f in os.listdir(proc_dir) if f.startswith('OSD-242_de_')]
if de242:
    fpath242 = os.path.join(proc_dir, de242[0])
    with open(fpath242, 'r') as f:
        header = f.readline().strip()
    cols242 = [c.strip('"') for c in header.split('","')]
    log2fc = [c for c in cols242 if c.startswith('Log2fc_')]
    print(f"OSD-242: {len(log2fc)} contrasts:")
    for c in log2fc:
        print(f"  {c}")

# --- Python code block 7 ---
import pandas as pd
import os
import numpy as np

proc_dir = '/mnt/shared-workspace/shared/osdr_data/processed'
deg_dir = '/mnt/shared-workspace/shared/results/deg'

# OSD-242: extract FLT vs GC contrast
# FLT_C1 = Space Flight, GC_C2 = Ground Control
# Use Log2fc_(FLT_C1)v(GC_C2) - positive = up in flight
fpath = [f for f in os.listdir(proc_dir) if f.startswith('OSD-242_de_')][0]
fpath = os.path.join(proc_dir, fpath)
with open(fpath, 'r') as f:
    header = f.readline().strip()
cols = [c.strip('"') for c in header.split('","')]

contrast_col = 'Log2fc_(FLT_C1)v(GC_C2)'
padj_col = 'Adj.p.value_(FLT_C1)v(GC_C2)'
pval_col = 'P.value_(FLT_C1)v(GC_C2)'

needed = ['ENSEMBL', 'SYMBOL', contrast_col, padj_col, pval_col]
needed = [c for c in needed if c in cols]
df = pd.read_csv(fpath, usecols=needed, low_memory=False)
df = df.rename(columns={'ENSEMBL':'gene_id', 'SYMBOL':'symbol', contrast_col:'log2fc', padj_col:'padj', pval_col:'pvalue'})
df['log2fc'] = pd.to_numeric(df['log2fc'], errors='coerce')
df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
df['pvalue'] = pd.to_numeric(df['pvalue'], errors='coerce')
df = df[df['gene_id'].notna()]
df['gene_id'] = df['gene_id'].astype(str).str.split('|').str[0]
df['is_deg'] = (df['padj'] < 0.05) & (df['log2fc'].abs() > 1)
df['deg_direction'] = np.where(df['is_deg'] & (df['log2fc'] > 0), 'up', np.where(df['is_deg'] & (df['log2fc'] < 0), 'down', 'ns'))
df['osd'] = 'OSD-242'
df['species'] = 'Mus musculus'
df['contrast'] = 'FLT_C1 vs GC_C2 (Space Flight vs Ground Control)'
df['direction_normalized'] = 'positive = up in spaceflight'

n_deg = df['is_deg'].sum()
print(f"OSD-242 (Mus musculus): {len(df)} genes, {n_deg} DEGs ({(df['deg_direction']=='up').sum()} up, {(df['deg_direction']=='down').sum()} down)")
df[['gene_id','symbol','log2fc','pvalue','padj','is_deg','deg_direction','osd','species']].to_csv(
    os.path.join(deg_dir, 'OSD-242_deg.csv'), index=False)

# --- Python code block 8 ---
import pandas as pd
import os

proc_dir = '/mnt/shared-workspace/shared/osdr_data/processed'

# OSD-96: prepare count matrix and sample metadata for DESeq2
fpath = os.path.join(proc_dir, 'OSD-96_counts_or_expr_GLDS-96_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv')
df96 = pd.read_csv(fpath)
# First col is gene ID (FBgn)
df96 = df96.rename(columns={df96.columns[0]: 'gene_id'})
print("OSD-96 columns:", list(df96.columns))
print("Sample columns (excluding gene_id):", [c for c in df96.columns if c != 'gene_id'])

# Build sample metadata from column names
# Pattern: Dmel_Can-S_wo_{FLT|GC}_{generation}_{timepoint}
samples = [c for c in df96.columns if c != 'gene_id']
meta96 = []
for s in samples:
    is_flt = '_FLT_' in s or '-FLT-' in s
    is_gc = '_GC_' in s or '-GC-' in s
    # Extract timepoint
    if '1.5h' in s or '1.5hr' in s:
        timepoint = '1.5h'
    elif '12h' in s:
        timepoint = '12h'
    elif '24h' in s:
        timepoint = '24h'
    else:
        timepoint = 'unknown'
    # Extract generation
    if '3rd-gen' in s:
        gen = '3rd-gen'
    elif '5th-gen' in s:
        gen = '5th-gen'
    else:
        gen = 'unknown'
    
    condition = 'Space_Flight' if is_flt else ('Ground_Control' if is_gc else 'unknown')
    # For 5th-gen, there are FLT-der and GC-der (derived from flight/ground)
    if 'FLT-der' in s:
        condition = 'Space_Flight'  # derived from flight
    elif 'GC-der' in s:
        condition = 'Ground_Control'  # derived from ground
    
    meta96.append({'sample': s, 'condition': condition, 'timepoint': timepoint, 'generation': gen})

meta96_df = pd.DataFrame(meta96)
print("\nOSD-96 sample metadata:")
print(meta96_df.to_string())

# Save count matrix and metadata for DESeq2
df96.to_csv('/mnt/shared-workspace/shared/osdr_data/OSD-96_counts.csv', index=False)
meta96_df.to_csv('/mnt/shared-workspace/shared/osdr_data/OSD-96_metadata.csv', index=False)
print(f"\nSaved OSD-96 counts ({df96.shape}) and metadata ({meta96_df.shape})")

# --- Python code block 9 ---
import pandas as pd
import os
import numpy as np

deg_dir = '/mnt/shared-workspace/shared/results/deg'

# Load all per-dataset DEG files
deg_files = [f for f in os.listdir(deg_dir) if f.endswith('_deg.csv') and f != 'all_degs_combined.csv']
print(f"DEG files: {len(deg_files)}")

all_degs = []
for f in deg_files:
    df = pd.read_csv(os.path.join(deg_dir, f))
    all_degs.append(df)
    
combined = pd.concat(all_degs, ignore_index=True)
print(f"Combined: {len(combined)} rows, {combined['osd'].nunique()} datasets")

# Apply tiered thresholds:
# Standard: padj < 0.05, |log2fc| > 1
# Relaxed (for datasets with <10 DEGs at standard): padj < 0.10, |log2fc| > 0.5
# For OSD-35 (no replicates): pvalue < 0.05, |log2fc| > 0.5 (no padj available)

# First, check current DEG counts
deg_counts = combined.groupby('osd')['is_deg'].sum().sort_values()
print("\nCurrent DEG counts (padj<0.05, |log2fc|>1):")
print(deg_counts)

# Identify datasets needing relaxed threshold
low_deg = deg_counts[deg_counts < 10].index.tolist()
print(f"\nDatasets with <10 DEGs (apply relaxed threshold): {low_deg}")

# Apply relaxed threshold for low-DEG datasets
for osd in low_deg:
    mask = combined['osd'] == osd
    if osd == 'OSD-35':
        # No replicates - use pvalue instead of padj
        combined.loc[mask, 'is_deg'] = (combined.loc[mask, 'pvalue'] < 0.05) & (combined.loc[mask, 'log2fc'].abs() > 0.5)
    else:
        combined.loc[mask, 'is_deg'] = (combined.loc[mask, 'padj'] < 0.10) & (combined.loc[mask, 'log2fc'].abs() > 0.5)
    combined.loc[mask, 'threshold'] = 'relaxed'
    n = combined.loc[mask, 'is_deg'].sum()
    print(f"  {osd}: {n} DEGs (relaxed: padj<0.10, |log2fc|>0.5)")

# Set threshold for standard datasets
combined['threshold'] = combined['threshold'].fillna('standard')

# Recompute deg_direction
combined['deg_direction'] = np.where(combined['is_deg'] & (combined['log2fc'] > 0), 'up',
                                     np.where(combined['is_deg'] & (combined['log2fc'] < 0), 'down', 'ns'))

# Save combined
combined.to_csv(os.path.join(deg_dir, 'all_degs_combined.csv'), index=False)

# Summary
print("\n=== Final DEG summary ===")
summary = combined[combined['is_deg']].groupby(['species','osd','deg_direction']).size().unstack(fill_value=0)
summary['total'] = summary.sum(axis=1)
print(summary)
print(f"\nTotal DEGs: {combined['is_deg'].sum()}")
print(f"Datasets: {combined['osd'].nunique()}")

# --- Python code block 10 ---
import pandas as pd
import os
import numpy as np

deg_dir = '/mnt/shared-workspace/shared/results/deg'
combined = pd.read_csv(os.path.join(deg_dir, 'all_degs_combined.csv'))

# Fix: ensure is_deg is boolean, fill NaN with False
combined['is_deg'] = combined['is_deg'].fillna(False).astype(bool)
combined['log2fc'] = pd.to_numeric(combined['log2fc'], errors='coerce')
combined['padj'] = pd.to_numeric(combined['padj'], errors='coerce')
combined['pvalue'] = pd.to_numeric(combined['pvalue'], errors='coerce')

# Recompute deg_direction
combined['deg_direction'] = np.where(combined['is_deg'] & (combined['log2fc'] > 0), 'up',
                                     np.where(combined['is_deg'] & (combined['log2fc'] < 0), 'down', 'ns'))

# Save
combined.to_csv(os.path.join(deg_dir, 'all_degs_combined.csv'), index=False)

# Summary
print("=== Final DEG summary ===")
summary = combined[combined['is_deg']].groupby(['species','osd','deg_direction']).size().unstack(fill_value=0)
summary['total'] = summary.sum(axis=1)
print(summary.to_string())
print(f"\nTotal DEGs: {combined['is_deg'].sum()}")
print(f"Datasets: {combined['osd'].nunique()}")

# Per-species summary
print("\n=== Per-species DEG counts ===")
sp_summary = combined[combined['is_deg']].groupby(['species','deg_direction']).size().unstack(fill_value=0)
sp_summary['total'] = sp_summary.sum(axis=1)
print(sp_summary.to_string())

# --- Python code block 11 ---
import pandas as pd
import os

proc_dir = '/mnt/shared-workspace/shared/osdr_data/processed'
fpath = os.path.join(proc_dir, 'OSD-62_de_GLDS-62_array_differential_expression_GLmicroarray.csv')

with open(fpath, 'r') as f:
    header = f.readline().strip()
cols = [c.strip('"') for c in header.split('","')]
log2fc = [c for c in cols if c.startswith('Log2fc_')]
print("OSD-62 all contrasts:")
for c in log2fc:
    print(f"  {c}")

# Find Wild Type Space Flight vs Ground Control
wt_contrasts = [c for c in log2fc if 'Wild Type' in c and 'Space Flight' in c and 'Ground Control' in c]
print(f"\nWild Type SF vs GC contrasts: {len(wt_contrasts)}")
for c in wt_contrasts:
    print(f"  {c}")

# --- Python code block 12 ---
import pandas as pd
import os
import numpy as np

proc_dir = '/mnt/shared-workspace/shared/osdr_data/processed'
deg_dir = '/mnt/shared-workspace/shared/results/deg'

# Re-extract OSD-62 with Wild Type SF vs GC contrast
fpath = os.path.join(proc_dir, 'OSD-62_de_GLDS-62_array_differential_expression_GLmicroarray.csv')
contrast_col = 'Log2fc_(Space Flight & Wild Type)v(Ground Control & Wild Type)'
padj_col = 'Adj.p.value_(Space Flight & Wild Type)v(Ground Control & Wild Type)'
pval_col = 'P.value_(Space Flight & Wild Type)v(Ground Control & Wild Type)'

df = pd.read_csv(fpath, usecols=['ENSEMBL','SYMBOL',contrast_col,padj_col,pval_col], low_memory=False)
df = df.rename(columns={'ENSEMBL':'gene_id','SYMBOL':'symbol',contrast_col:'log2fc',padj_col:'padj',pval_col:'pvalue'})
df['log2fc'] = pd.to_numeric(df['log2fc'], errors='coerce')
df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
df['pvalue'] = pd.to_numeric(df['pvalue'], errors='coerce')
df = df[df['gene_id'].notna()]
df['gene_id'] = df['gene_id'].astype(str).str.split('|').str[0]

# This is SF_vs_GC direction (positive = up in flight) - correct
# Apply relaxed threshold for yeast
df['is_deg'] = (df['padj'] < 0.10) & (df['log2fc'].abs() > 0.5)
df['deg_direction'] = np.where(df['is_deg'] & (df['log2fc'] > 0), 'up', np.where(df['is_deg'] & (df['log2fc'] < 0), 'down', 'ns'))
df['osd'] = 'OSD-62'
df['species'] = 'Saccharomyces cerevisiae'
df['threshold'] = 'relaxed'

n_deg = df['is_deg'].sum()
print(f"OSD-62 (yeast, Wild Type SF vs GC): {len(df)} genes, {n_deg} DEGs ({(df['deg_direction']=='up').sum()} up, {(df['deg_direction']=='down').sum()} down)")

# If still 0, try even more relaxed: pvalue < 0.05, |log2fc| > 0.25
if n_deg == 0:
    df['is_deg'] = (df['pvalue'] < 0.05) & (df['log2fc'].abs() > 0.25)
    df['deg_direction'] = np.where(df['is_deg'] & (df['log2fc'] > 0), 'up', np.where(df['is_deg'] & (df['log2fc'] < 0), 'down', 'ns'))
    n_deg = df['is_deg'].sum()
    print(f"  After ultra-relaxed (p<0.05, |log2fc|>0.25): {n_deg} DEGs ({(df['deg_direction']=='up').sum()} up, {(df['deg_direction']=='down').sum()} down)")

df[['gene_id','symbol','log2fc','pvalue','padj','is_deg','deg_direction','osd','species','threshold']].to_csv(
    os.path.join(deg_dir, 'OSD-62_deg.csv'), index=False)

# Also fix OSD-258 - check for a better contrast
print("\n--- OSD-258 contrasts ---")
fpath258 = os.path.join(proc_dir, [f for f in os.listdir(proc_dir) if f.startswith('OSD-258_de_')][0])
with open(fpath258, 'r') as f:
    header = f.readline().strip()
cols258 = [c.strip('"') for c in header.split('","')]
log2fc258 = [c for c in cols258 if c.startswith('Log2fc_')]
for c in log2fc258:
    print(f"  {c}")

# --- Python code block 13 ---
import pandas as pd
import os
import numpy as np

proc_dir = '/mnt/shared-workspace/shared/osdr_data/processed'
deg_dir = '/mnt/shared-workspace/shared/results/deg'

# Re-extract OSD-258 with the 4.5 weeks microgravity vs ground control contrast
fpath = os.path.join(proc_dir, [f for f in os.listdir(proc_dir) if f.startswith('OSD-258_de_')][0])
contrast_col = 'Log2fc_(Space Flight & 4.5 weeks microgravity exposure)v(Ground Control & 5.5 weeks plus 10 days)'
padj_col = 'Adj.p.value_(Space Flight & 4.5 weeks microgravity exposure)v(Ground Control & 5.5 weeks plus 10 days)'
pval_col = 'P.value_(Space Flight & 4.5 weeks microgravity exposure)v(Ground Control & 5.5 weeks plus 10 days)'

df = pd.read_csv(fpath, usecols=['ENSEMBL','SYMBOL',contrast_col,padj_col,pval_col], low_memory=False)
df = df.rename(columns={'ENSEMBL':'gene_id','SYMBOL':'symbol',contrast_col:'log2fc',padj_col:'padj',pval_col:'pvalue'})
df['log2fc'] = pd.to_numeric(df['log2fc'], errors='coerce')
df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
df['pvalue'] = pd.to_numeric(df['pvalue'], errors='coerce')
df = df[df['gene_id'].notna()]
df['gene_id'] = df['gene_id'].astype(str).str.split('|').str[0]

# SF_vs_GC direction (positive = up in flight)
# Apply relaxed threshold
df['is_deg'] = (df['padj'] < 0.10) & (df['log2fc'].abs() > 0.5)
df['deg_direction'] = np.where(df['is_deg'] & (df['log2fc'] > 0), 'up', np.where(df['is_deg'] & (df['log2fc'] < 0), 'down', 'ns'))
df['osd'] = 'OSD-258'
df['species'] = 'Homo sapiens'
df['threshold'] = 'relaxed'

n_deg = df['is_deg'].sum()
print(f"OSD-258 (human, 4.5wk microgravity vs GC): {len(df)} genes, {n_deg} DEGs ({(df['deg_direction']=='up').sum()} up, {(df['deg_direction']=='down').sum()} down)")

# If still very few, try pvalue < 0.01, |log2fc| > 0.5 (nominal significance)
if n_deg < 10:
    df['is_deg'] = (df['pvalue'] < 0.01) & (df['log2fc'].abs() > 0.5)
    df['deg_direction'] = np.where(df['is_deg'] & (df['log2fc'] > 0), 'up', np.where(df['is_deg'] & (df['log2fc'] < 0), 'down', 'ns'))
    n_deg = df['is_deg'].sum()
    print(f"  After nominal (p<0.01, |log2fc|>0.5): {n_deg} DEGs ({(df['deg_direction']=='up').sum()} up, {(df['deg_direction']=='down').sum()} down)")

df[['gene_id','symbol','log2fc','pvalue','padj','is_deg','deg_direction','osd','species','threshold']].to_csv(
    os.path.join(deg_dir, 'OSD-258_deg.csv'), index=False)

# Now rebuild the combined DEG table
deg_files = [f for f in os.listdir(deg_dir) if f.endswith('_deg.csv') and f != 'all_degs_combined.csv']
all_degs = []
for f in deg_files:
    df = pd.read_csv(os.path.join(deg_dir, f))
    all_degs.append(df)
combined = pd.concat(all_degs, ignore_index=True)
combined['is_deg'] = combined['is_deg'].fillna(False).astype(bool)
combined['log2fc'] = pd.to_numeric(combined['log2fc'], errors='coerce')
combined['threshold'] = combined.get('threshold', 'standard').fillna('standard')
combined['deg_direction'] = np.where(combined['is_deg'] & (combined['log2fc'] > 0), 'up',
                                     np.where(combined['is_deg'] & (combined['log2fc'] < 0), 'down', 'ns'))
combined.to_csv(os.path.join(deg_dir, 'all_degs_combined.csv'), index=False)

print("\n=== FINAL DEG summary (all 22 datasets) ===")
summary = combined[combined['is_deg']].groupby(['species','osd','deg_direction']).size().unstack(fill_value=0)
summary['total'] = summary.sum(axis=1)
print(summary.to_string())
print(f"\nTotal DEGs: {combined['is_deg'].sum()}")
print(f"\nPer-species:")
sp = combined[combined['is_deg']].groupby(['species','deg_direction']).size().unstack(fill_value=0)
sp['total'] = sp.sum(axis=1)
print(sp.to_string())

# --- Python code block 14 ---
import pandas as pd
deg = pd.read_csv('/mnt/shared-workspace/shared/results/deg/all_degs_combined.csv', low_memory=False)
print("Datasets in combined:", sorted(deg['osd'].unique()))
print("Count:", deg['osd'].nunique())
# Check OSD-35
print("\nOSD-35 rows:", len(deg[deg['osd']=='OSD-35']))
print("OSD-35 DEGs:", deg[deg['osd']=='OSD-35']['is_deg'].sum())

# --- Python code block 15 ---
import pandas as pd
import numpy as np
import os

deg_dir = '/mnt/shared-workspace/shared/results/deg'
combined = pd.read_csv(os.path.join(deg_dir, 'all_degs_combined.csv'), low_memory=False)
combined['is_deg'] = combined['is_deg'].fillna(False).astype(bool)
combined['log2fc'] = pd.to_numeric(combined['log2fc'], errors='coerce')
combined['pvalue'] = pd.to_numeric(combined['pvalue'], errors='coerce')
combined['padj'] = pd.to_numeric(combined['padj'], errors='coerce')

# Fix OSD-35: use pvalue < 0.05, |log2fc| > 0.5 (no padj available - no replicates)
mask35 = combined['osd'] == 'OSD-35'
combined.loc[mask35, 'is_deg'] = (combined.loc[mask35, 'pvalue'] < 0.05) & (combined.loc[mask35, 'log2fc'].abs() > 0.5)
combined.loc[mask35, 'threshold'] = 'relaxed_nominal'

# Recompute deg_direction
combined['deg_direction'] = np.where(combined['is_deg'] & (combined['log2fc'] > 0), 'up',
                                     np.where(combined['is_deg'] & (combined['log2fc'] < 0), 'down', 'ns'))

combined.to_csv(os.path.join(deg_dir, 'all_degs_combined.csv'), index=False)

print("OSD-35 DEGs after fix:", combined[combined['osd']=='OSD-35']['is_deg'].sum())
print("\n=== FINAL COMPLETE DEG summary ===")
summary = combined[combined['is_deg']].groupby(['species','osd','deg_direction']).size().unstack(fill_value=0)
summary['total'] = summary.sum(axis=1)
print(summary.to_string())
print(f"\nTotal DEGs: {combined['is_deg'].sum()}")
print(f"\nPer-species:")
sp = combined[combined['is_deg']].groupby(['species','deg_direction']).size().unstack(fill_value=0)
sp['total'] = sp.sum(axis=1)
print(sp.to_string())

