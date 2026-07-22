# ======================================================================
# Step 1: Data Acquisition from NASA OSDR
# 
# Description: Query NASA OSDR Biological Data API, select datasets, download processed data
# 
# Inputs: NASA OSDR API (https://visualization.osdr.nasa.gov/biodata/api/v2/query/)
# Outputs: data/dataset_selection.csv, data/contrast_selection.csv, data/all_sample_metadata.csv
# 
# Language: PYTHON
# Extracted from analysis notebook (worker-0.ipynb)
# See METHODS.md for full parameter details.
# ======================================================================

# --- Python code block 1 ---
import pandas as pd
df = pd.read_csv('/tmp/osdr_assays.csv')
print("Total rows:", len(df))
print("\nUnique assay technology types:")
print(df['investigation.study assays.study assay technology type'].value_counts())
print("\nUnique organisms:")
print(df['study.characteristics.organism'].value_counts())

# --- Python code block 2 ---
import pandas as pd
df = pd.read_csv('/tmp/osdr_assays.csv')
# Filter to RNA-Seq (bulk) only
rnaseq = df[df['investigation.study assays.study assay technology type'] == 'RNA Sequencing (RNA-Seq)'].copy()
print("RNA-Seq rows:", len(rnaseq))
print("\nRNA-Seq organisms:")
print(rnaseq['study.characteristics.organism'].value_counts())

# Our 6 target species
targets = ['Mus musculus', 'Homo sapiens', 'Drosophila melanogaster', 
           'Caenorhabditis elegans', 'Saccharomyces cerevisiae', 'Arabidopsis thaliana']
print("\n=== RNA-Seq datasets for target species ===")
for org in targets:
    sub = rnaseq[rnaseq['study.characteristics.organism'] == org]
    datasets = sub['id.accession'].unique()
    factors = sub['study.factor value.spaceflight'].unique()
    print(f"\n{org}: {len(datasets)} datasets -> {list(datasets)}")
    print(f"  spaceflight factors: {list(factors)}")

# --- Python code block 3 ---
import pandas as pd
df = pd.read_csv('/tmp/osdr_assays.csv')

# Check yeast across ALL assay types
yeast_all = df[df['study.characteristics.organism'] == 'Saccharomyces cerevisiae']
print("=== Saccharomyces cerevisiae - ALL assay types ===")
print(yeast_all['investigation.study assays.study assay technology type'].value_counts())
print("Datasets:", yeast_all['id.accession'].unique().tolist())

# Check C. elegans across ALL assay types
ce_all = df[df['study.characteristics.organism'] == 'Caenorhabditis elegans']
print("\n=== Caenorhabditis elegans - ALL assay types ===")
print(ce_all['investigation.study assays.study assay technology type'].value_counts())
print("Datasets:", ce_all['id.accession'].unique().tolist())
print(ce_all[['id.accession','investigation.study assays.study assay technology type','study.factor value.spaceflight']].to_string())

# --- Python code block 4 ---
import pandas as pd
import json

df = pd.read_csv('/tmp/osdr_assays.csv')
rnaseq = df[df['investigation.study assays.study assay technology type'] == 'RNA Sequencing (RNA-Seq)'].copy()
microarray = df[df['investigation.study assays.study assay technology type'] == 'DNA microarray'].copy()

# Build dataset selection
# For mammals (mouse, human): select ~5-6 datasets covering multiple tissues
# For Arabidopsis: select ~4 datasets
# For Drosophila: select 3-4 datasets
# For C. elegans: use microarray (OSD-35, 42, 112, 113)
# For yeast: use microarray OSD-62

selection = {
    'Mus musculus': {
        'rnaseq': ['OSD-47', 'OSD-48', 'OSD-104', 'OSD-105', 'OSD-252' if False else 'OSD-242', 'OSD-379', 'OSD-421'],
        'microarray': []
    },
    'Homo sapiens': {
        'rnaseq': ['OSD-258', 'OSD-323', 'OSD-481', 'OSD-530', 'OSD-571'],
        'microarray': []
    },
    'Arabidopsis thaliana': {
        'rnaseq': ['OSD-37', 'OSD-120', 'OSD-217', 'OSD-251'],
        'microarray': []
    },
    'Drosophila melanogaster': {
        'rnaseq': ['OSD-96', 'OSD-207', 'OSD-347'],
        'microarray': []
    },
    'Caenorhabditis elegans': {
        'rnaseq': [],
        'microarray': ['OSD-35', 'OSD-42', 'OSD-112', 'OSD-113']
    },
    'Saccharomyces cerevisiae': {
        'rnaseq': [],
        'microarray': ['OSD-62']
    }
}

# Verify all selected datasets exist in our query results
all_rnaseq_datasets = set(rnaseq['id.accession'].unique())
all_microarray_datasets = set(microarray['id.accession'].unique())

print("=== Dataset selection verification ===")
total = 0
for species, dsets in selection.items():
    valid_rnaseq = [d for d in dsets['rnaseq'] if d in all_rnaseq_datasets]
    valid_micro = [d for d in dsets['microarray'] if d in all_microarray_datasets]
    missing_rnaseq = [d for d in dsets['rnaseq'] if d not in all_rnaseq_datasets]
    missing_micro = [d for d in dsets['microarray'] if d not in all_microarray_datasets]
    n = len(valid_rnaseq) + len(valid_micro)
    total += n
    print(f"\n{species}: {n} datasets")
    print(f"  RNA-seq: {valid_rnaseq}" + (f"  [MISSING: {missing_rnaseq}]" if missing_rnaseq else ""))
    print(f"  Microarray: {valid_micro}" + (f"  [MISSING: {missing_micro}]" if missing_micro else ""))
print(f"\nTotal datasets selected: {total}")

# Save selection
with open('/mnt/shared-workspace/shared/osdr_data/dataset_selection.json', 'w') as f:
    json.dump(selection, f, indent=2)
print("Saved selection to dataset_selection.json")

# --- Python code block 5 ---
import pandas as pd
import json
import time
import urllib.request
import urllib.parse

with open('/mnt/shared-workspace/shared/osdr_data/dataset_selection.json') as f:
    selection = json.load(f)

# Flatten all selected datasets
all_datasets = []
for species, dsets in selection.items():
    for d in dsets['rnaseq']:
        all_datasets.append({'species': species, 'osd': d, 'assay': 'RNA-Seq'})
    for d in dsets['microarray']:
        all_datasets.append({'species': species, 'osd': d, 'assay': 'microarray'})

print(f"Querying metadata for {len(all_datasets)} datasets...")

results = []
for ds in all_datasets:
    osd = ds['osd']
    # Query metadata: organism, tissue, strain, spaceflight factor, assay name
    url = f"https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/?id.accession={osd}&study.characteristics.organism&study.characteristics.tissue&study.characteristics.strain&study.factor%20value.spaceflight&study.characteristics.age&study.characteristics.sex&format=csv"
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=30) as resp:
            text = resp.read().decode('utf-8')
        from io import StringIO
        meta = pd.read_csv(StringIO(text))
        meta['selected_species'] = ds['species']
        meta['selected_assay'] = ds['assay']
        results.append(meta)
        tissues = meta['study.characteristics.tissue'].dropna().unique().tolist()
        factors = meta['study.factor value.spaceflight'].dropna().unique().tolist()
        strains = meta['study.characteristics.strain'].dropna().unique().tolist()
        n_samples = len(meta['id.sample name'].unique())
        print(f"  {osd} ({ds['species']}): {n_samples} samples, tissues={tissues[:3]}, factors={factors[:3]}")
    except Exception as e:
        print(f"  {osd}: ERROR {e}")
    time.sleep(0.3)

all_meta = pd.concat(results, ignore_index=True)
all_meta.to_csv('/mnt/shared-workspace/shared/osdr_data/all_sample_metadata.csv', index=False)
print(f"\nTotal samples: {len(all_meta)}")
print(f"Saved to all_sample_metadata.csv")

# Summary table
summary = all_meta.groupby(['selected_species','selected_assay','id.accession']).agg(
    n_samples=('id.sample name','nunique'),
    tissues=('study.characteristics.tissue', lambda x: '|'.join(sorted(set(x.dropna())))),
    strains=('study.characteristics.strain', lambda x: '|'.join(sorted(set(x.dropna())))[:60]),
    factors=('study.factor value.spaceflight', lambda x: '|'.join(sorted(set(x.dropna()))))
).reset_index()
summary.to_csv('/mnt/shared-workspace/shared/osdr_data/dataset_summary.csv', index=False)
print("\n=== Dataset summary ===")
print(summary.to_string())

# --- Python code block 6 ---
import pandas as pd
import json
import urllib.request
from io import StringIO

# Check full metadata structure for one mouse dataset to find tissue field
url = "https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/?id.accession=OSD-379&study.characteristics&format=csv"
with urllib.request.urlopen(url, timeout=30) as resp:
    text = resp.read().decode('utf-8')
meta = pd.read_csv(StringIO(text))
print("OSD-379 columns:")
for c in meta.columns:
    vals = meta[c].dropna().unique()[:3]
    print(f"  {c}: {list(vals)}")

# --- Python code block 7 ---
import pandas as pd
import json
import urllib.request
import time
from io import StringIO

with open('/mnt/shared-workspace/shared/osdr_data/dataset_selection.json') as f:
    selection = json.load(f)

all_datasets = []
for species, dsets in selection.items():
    for d in dsets['rnaseq']:
        all_datasets.append({'species': species, 'osd': d, 'assay': 'RNA-Seq'})
    for d in dsets['microarray']:
        all_datasets.append({'species': species, 'osd': d, 'assay': 'microarray'})

# Query material type + factor for each
print("=== Tissue (material type) + spaceflight factor per dataset ===")
refined = []
for ds in all_datasets:
    osd = ds['osd']
    url = f"https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/?id.accession={osd}&study.characteristics.material%20type&study.factor%20value.spaceflight&study.characteristics.organism&format=csv"
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            text = resp.read().decode('utf-8')
        meta = pd.read_csv(StringIO(text))
        tissues = meta['study.characteristics.material type'].dropna().unique().tolist()
        factors = meta['study.factor value.spaceflight'].dropna().unique().tolist()
        has_flt = any('Space Flight' in str(f) for f in factors)
        has_gc = any('Ground Control' in str(f).lower().replace('control','Control') for f in factors) or any('Ground control' in str(f) for f in factors)
        usable = has_flt and has_gc
        print(f"  {osd} ({ds['species']}, {ds['assay']}): usable={usable}")
        print(f"    tissues: {tissues}")
        print(f"    factors: {factors}")
        refined.append({'species': ds['species'], 'osd': osd, 'assay': ds['assay'], 
                       'tissues': tissues, 'factors': factors, 'usable': usable,
                       'n_samples': len(meta['id.sample name'].unique())})
    except Exception as e:
        print(f"  {osd}: ERROR {e}")
    time.sleep(0.3)

refined_df = pd.DataFrame(refined)
refined_df.to_csv('/mnt/shared-workspace/shared/osdr_data/dataset_refined.csv', index=False)
print(f"\nUsable datasets: {refined_df['usable'].sum()}/{len(refined_df)}")
print("\nUnusable:")
print(refined_df[~refined_df['usable']][['species','osd','assay','factors']].to_string())

# --- Python code block 8 ---
import pandas as pd
import json

refined_df = pd.read_csv('/mnt/shared-workspace/shared/osdr_data/dataset_refined.csv')
import ast
refined_df['factors'] = refined_df['factors'].apply(ast.literal_eval)

# Correct usable check: has Space Flight AND Ground Control (case-insensitive)
def is_usable(factors):
    has_flt = any('space flight' in str(f).lower() for f in factors)
    has_gc = any('ground control' in str(f).lower() for f in factors)
    return has_flt and has_gc

refined_df['usable'] = refined_df['factors'].apply(is_usable)
print(f"Usable datasets: {refined_df['usable'].sum()}/{len(refined_df)}")
print("\nUnusable (no flight+ground):")
print(refined_df[~refined_df['usable']][['species','osd','assay','factors']].to_string())

# Drop unusable
usable = refined_df[refined_df['usable']].copy()
print(f"\n=== Final usable selection: {len(usable)} datasets ===")
print(usable.groupby('species').size())
print("\nTissues per species:")
for sp in usable['species'].unique():
    sub = usable[usable['species']==sp]
    tissues = []
    for t in sub['tissues']:
        tissues.extend(ast.literal_eval(t) if isinstance(t,str) else [])
    print(f"  {sp}: {sorted(set(tissues))}")

usable.to_csv('/mnt/shared-workspace/shared/osdr_data/dataset_final.csv', index=False)

# --- Python code block 9 ---
import pandas as pd
df = pd.read_csv('/tmp/osdr_assays.csv')
rnaseq = df[df['investigation.study assays.study assay technology type'] == 'RNA Sequencing (RNA-Seq)'].copy()
human_rnaseq = rnaseq[rnaseq['study.characteristics.organism'] == 'Homo sapiens']
print("All human RNA-seq datasets:", human_rnaseq['id.accession'].unique().tolist())
# We already have OSD-258, OSD-323 (both Cells, flight+ground)
# Let's check OSD-431, OSD-483, OSD-635, OSD-684, OSD-781, OSD-811, OSD-863, OSD-867, OSD-871, OSD-903, OSD-936, OSD-937
import urllib.request, time
from io import StringIO
candidates = ['OSD-431','OSD-483','OSD-635','OSD-684','OSD-781','OSD-811','OSD-863','OSD-867','OSD-871','OSD-903','OSD-936','OSD-937']
for osd in candidates:
    url = f"https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/?id.accession={osd}&study.characteristics.material%20type&study.factor%20value.spaceflight&format=csv"
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            text = resp.read().decode('utf-8')
        meta = pd.read_csv(StringIO(text))
        tissues = meta['study.characteristics.material type'].dropna().unique().tolist()
        factors = meta['study.factor value.spaceflight'].dropna().unique().tolist()
        n = len(meta['id.sample name'].unique())
        print(f"  {osd}: n={n}, tissues={tissues}, factors={factors}")
    except Exception as e:
        print(f"  {osd}: ERR {e}")
    time.sleep(0.3)

# --- Python code block 10 ---
import pandas as pd
import ast

refined_df = pd.read_csv('/mnt/shared-workspace/shared/osdr_data/dataset_refined.csv')
refined_df['factors'] = refined_df['factors'].apply(ast.literal_eval)
refined_df['tissues'] = refined_df['tissues'].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else [])

def is_usable(factors):
    has_flt = any('space flight' in str(f).lower() for f in factors)
    has_gc = any('ground control' in str(f).lower() for f in factors)
    return has_flt and has_gc
refined_df['usable'] = refined_df['factors'].apply(is_usable)
usable = refined_df[refined_df['usable']].copy()

# Add 2 more human datasets
new_human = pd.DataFrame([
    {'species':'Homo sapiens','osd':'OSD-684','assay':'RNA-Seq','tissues':['Myoblasts'],'factors':['Ground Control','Space Flight'],'usable':True,'n_samples':12},
    {'species':'Homo sapiens','osd':'OSD-863','assay':'RNA-Seq','tissues':['Cells, Cultured'],'factors':['Ground Control','Space Flight'],'usable':True,'n_samples':19},
])
final = pd.concat([usable, new_human], ignore_index=True)
final = final.sort_values(['species','osd']).reset_index(drop=True)
final.to_csv('/mnt/shared-workspace/shared/osdr_data/dataset_final.csv', index=False)

print(f"=== FINAL SELECTION: {len(final)} datasets ===")
print(final.groupby('species').size())
print("\nTissues per species:")
for sp in final['species'].unique():
    sub = final[final['species']==sp]
    tissues = []
    for t in sub['tissues']:
        tissues.extend(t if isinstance(t,list) else [])
    print(f"  {sp} ({len(sub)} datasets): {sorted(set(tissues))}")
print("\nFull table:")
print(final[['species','osd','assay','n_samples','tissues']].to_string())

# --- Python code block 11 ---
import pandas as pd
import json
import urllib.request
import urllib.parse
import time
import re
from io import StringIO

final = pd.read_csv('/mnt/shared-workspace/shared/osdr_data/dataset_final.csv')

# For each dataset, find the differential expression table filename
# Pattern: GLDS-{n}_{rnaseq|array}_differential_expression_{pipeline}.csv
# Query the file metadata to find DE table + sample table + normalized expression

de_files = {}
sampletable_files = {}
norm_files = {}

for _, row in final.iterrows():
    osd = row['osd']
    # Query file names + data types
    url = f"https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/?id.accession={osd}&file.file%20name&file.data%20type&format=csv"
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            text = resp.read().decode('utf-8')
        fmeta = pd.read_csv(StringIO(text))
        # Get unique file names
        fnames = fmeta['file.file name'].dropna().unique().tolist()
        ftypes = fmeta['file.data type'].dropna().unique().tolist()
        
        # Find DE table (differential_expression)
        de_fn = [f for f in fnames if 'differential_expression' in f.lower() or 'differential expression' in f.lower()]
        st_fn = [f for f in fnames if 'sampletable' in f.lower() or 'sample_table' in f.lower() or 'sample table' in f.lower()]
        norm_fn = [f for f in fnames if 'normalized_expression' in f.lower() or 'normalized_expression' in f.lower() or ('normalized' in f.lower() and 'expression' in f.lower())]
        # Also look for raw counts for RNA-seq
        counts_fn = [f for f in fnames if 'raw_counts' in f.lower() or 'unnormalized' in f.lower() or ('count' in f.lower() and 'raw' in f.lower())]
        
        de_files[osd] = de_fn[0] if de_fn else None
        sampletable_files[osd] = st_fn[0] if st_fn else None
        norm_files[osd] = norm_fn[0] if norm_fn else None
        
        print(f"{osd}: DE={de_fn[0] if de_fn else 'NONE'}, ST={'Y' if st_fn else 'N'}, NORM={'Y' if norm_fn else 'N'}, COUNTS={'Y' if counts_fn else 'N'}")
    except Exception as e:
        print(f"{osd}: ERROR {e}")
    time.sleep(0.3)

# Save file mapping
file_map = pd.DataFrame([
    {'osd': osd, 'de_file': de_files.get(osd), 'sampletable_file': sampletable_files.get(osd), 'norm_file': norm_files.get(osd)}
    for osd in final['osd']
])
file_map.to_csv('/mnt/shared-workspace/shared/osdr_data/file_mapping.csv', index=False)
print(f"\nDE tables found: {sum(1 for v in de_files.values() if v)}/{len(final)}")
print(f"Sample tables found: {sum(1 for v in sampletable_files.values() if v)}/{len(final)}")

# --- Python code block 12 ---
import pandas as pd
import urllib.request
from io import StringIO
import time

# Check OSD-35 and OSD-96 file lists in detail
for osd in ['OSD-35', 'OSD-96']:
    url = f"https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/?id.accession={osd}&file.file%20name&file.data%20type&format=csv"
    with urllib.request.urlopen(url, timeout=30) as resp:
        text = resp.read().decode('utf-8')
    fmeta = pd.read_csv(StringIO(text))
    fnames = fmeta['file.file name'].dropna().unique().tolist()
    print(f"\n=== {osd} files ({len(fnames)}) ===")
    for f in fnames:
        print(f"  {f}")
    time.sleep(0.3)

# --- Python code block 13 ---
import pandas as pd
import urllib.request
import urllib.parse
import os
import time

final = pd.read_csv('/mnt/shared-workspace/shared/osdr_data/dataset_final.csv')
file_map = pd.read_csv('/mnt/shared-workspace/shared/osdr_data/file_mapping.csv')

# Build download list: DE tables + sample tables + count/expression matrices for OSD-35, OSD-96
download_dir = '/mnt/shared-workspace/shared/osdr_data/processed'
os.makedirs(download_dir, exist_ok=True)

# Special files for OSD-35 (normalized expression) and OSD-96 (RSEM counts)
special_files = {
    'OSD-35': 'GLDS-35_array_normalized_expression_probeset_GLmicroarray.csv',
    'OSD-96': 'GLDS-96_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv',
}

downloads = []
for _, row in file_map.iterrows():
    osd = row['osd']
    if pd.notna(row['de_file']):
        downloads.append((osd, row['de_file'], 'de'))
    if pd.notna(row['sampletable_file']):
        downloads.append((osd, row['sampletable_file'], 'sampletable'))
    if osd in special_files:
        downloads.append((osd, special_files[osd], 'counts_or_expr'))

print(f"Total files to download: {len(downloads)}")

# Download function
def download_file(osd, filename, category):
    outpath = os.path.join(download_dir, f"{osd}_{category}_{filename}")
    if os.path.exists(outpath) and os.path.getsize(outpath) > 100:
        return outpath, 'cached'
    url = f"https://osdr.nasa.gov/geode-py/ws/studies/{osd}/download?file={urllib.parse.quote(filename)}"
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = resp.read()
        with open(outpath, 'wb') as f:
            f.write(data)
        return outpath, f'{len(data)} bytes'
    except Exception as e:
        return outpath, f'ERROR: {e}'

# Download all
results = []
for osd, fname, cat in downloads:
    path, status = download_file(osd, fname, cat)
    sz = os.path.getsize(path) if os.path.exists(path) else 0
    results.append({'osd': osd, 'category': cat, 'filename': fname, 'path': path, 'size_mb': round(sz/1e6,2), 'status': status})
    print(f"  {osd} {cat}: {status} ({round(sz/1e6,2)} MB)")
    time.sleep(0.2)

dl_df = pd.DataFrame(results)
dl_df.to_csv('/mnt/shared-workspace/shared/osdr_data/download_log.csv', index=False)
n_ok = sum(1 for r in results if 'ERROR' not in str(r.get('status', '')))
print(f"\nDownloaded: {n_ok} files")
print(f"Total size: {round(dl_df['size_mb'].sum(),2)} MB")

# --- Python code block 14 ---
import pandas as pd
import urllib.request
import urllib.parse
import os
import time

final = pd.read_csv('/mnt/shared-workspace/shared/osdr_data/dataset_final.csv')
file_map = pd.read_csv('/mnt/shared-workspace/shared/osdr_data/file_mapping.csv')

download_dir = '/mnt/shared-workspace/shared/osdr_data/processed'
os.makedirs(download_dir, exist_ok=True)

special_files = {
    'OSD-35': 'GLDS-35_array_normalized_expression_probeset_GLmicroarray.csv',
    'OSD-96': 'GLDS-96_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv',
}

downloads = []
for _, row in file_map.iterrows():
    osd = row['osd']
    if pd.notna(row['de_file']):
        downloads.append((osd, row['de_file'], 'de'))
    if pd.notna(row['sampletable_file']):
        downloads.append((osd, row['sampletable_file'], 'sampletable'))
    if osd in special_files:
        downloads.append((osd, special_files[osd], 'counts_or_expr'))

print(f"Total files to download: {len(downloads)}")

def download_file(osd, filename, category):
    outpath = os.path.join(download_dir, f"{osd}_{category}_{filename}")
    if os.path.exists(outpath) and os.path.getsize(outpath) > 100:
        return outpath, 'cached'
    url = f"https://osdr.nasa.gov/geode-py/ws/studies/{osd}/download?file={urllib.parse.quote(filename)}"
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = resp.read()
        with open(outpath, 'wb') as f:
            f.write(data)
        return outpath, f'{len(data)} bytes'
    except Exception as e:
        return outpath, f'ERROR: {e}'

results = []
for osd, fname, cat in downloads:
    path, status = download_file(osd, fname, cat)
    sz = os.path.getsize(path) if os.path.exists(path) else 0
    results.append({'osd': osd, 'category': cat, 'filename': fname, 'path': path, 'size_mb': round(sz/1e6,2), 'status': status})
    print(f"  {osd} {cat}: {status} ({round(sz/1e6,2)} MB)")
    time.sleep(0.2)

dl_df = pd.DataFrame(results)
dl_df.to_csv('/mnt/shared-workspace/shared/osdr_data/download_log.csv', index=False)
n_ok = sum(1 for r in results if 'ERROR' not in str(r[3]))
print(f"\nDownloaded OK: {n_ok}/{len(results)} files")
print(f"Total size: {round(dl_df['size_mb'].sum(),2)} MB")

