process GENERATE_COORDINATES {
    publishDir { "${params.outdir}/${publish_subdir}" },
        mode: 'copy'
    label 'process_gpu'
    conda "${moduleDir}/environment.yml"
    // container "containers/umap/umap.sif"

    input:
        tuple val(embeddings_dirs), val(excluded_genofeatures), val(nn), val(md), val(publish_subdir)
        
    output:
        path "${publish_subdir}", emit: coordinates_files
        path "versions.yml", emit: versions

    script:
    """
#!/usr/bin/env python3
import os
import torch
import numpy as np
import umap
import random
import pandas as pd
import json
from pathlib import Path

# Create output directories
os.makedirs("${publish_subdir}", exist_ok=True)

# SET SEED for reproducibility
os.environ['PYTHONHASHSEED'] = '42'
random.seed(42)
np.random.seed(42)
torch.manual_seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed(42)
    torch.cuda.manual_seed_all(42)

def load_embedding(file_path):
    try:
        embedding_data = torch.load(file_path, map_location='cpu')
        embedding = embedding_data["mean_representations"][36].numpy()
        embedding_id = embedding_data.get("label")
        return embedding, embedding_id
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None, None

def find_pt_files(base_dir):
    pt_files = []
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith('.pt'):
                pt_files.append(os.path.join(root, file))
    return pt_files

# Load embeddings
embeddings = []
embedding_ids = []

base_dirs = [${embeddings_dirs.collect { "'${it}'" }.join(',')}]

# Find all .pt files
pt_files = []
for base_dir in base_dirs:
    if os.path.isdir(base_dir):
        pt_files.extend(find_pt_files(base_dir))
print(f"Found {len(pt_files)} .pt files")

# Process each .pt file
excluded_raw = "${excluded_genofeatures}".strip()
excluded_subdirs = []
if excluded_raw and excluded_raw.lower() not in ['null', 'none']:
    excluded_subdirs = [
        s.strip().strip("'").strip('"')
        for s in excluded_raw.split(',')
        if s.strip()
    ]
for file_path in pt_files:
    genofeature = Path(file_path).parent.parent.name
    if genofeature.lower() not in [s.lower() for s in excluded_subdirs]:
        print(f"Processing: {file_path}")
        embedding, embedding_id = load_embedding(file_path)
        if embedding is not None and embedding_id is not None:
            embeddings.append(embedding)
            embedding_ids.append(embedding_id)
            print(f"Successfully loaded embedding from {file_path}")
    else:
        print(f"Skipping excluded genofeature: {genofeature} in {file_path}")

if len(embeddings) > 0:
    print(f"Processing {len(embeddings)} embeddings")
    embeddings = np.stack(embeddings)

    groups = {}
    for i, eid in enumerate(embedding_ids):
        contig_id = '_'.join(eid.split('_')[:-3])  # Parse embedding ID
        if contig_id not in groups:
            groups[contig_id] = []
        groups[contig_id].append(i)
    
    nn = int(${nn})
    md = float(${md})
    # Generate UMAP for parameter
    print(f"Generating UMAP with nn={nn}, md={md}")
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=nn,
        min_dist=md,
        metric='cosine',
        random_state=42,
        transform_seed=42,
        n_epochs=200
    )
    coordinates = reducer.fit_transform(embeddings)
    
    # Save outputs with parameter-specific names
    md_int = int(md * 10)
    coord_file = f"${publish_subdir}/coordinates_nn{nn}_md{md_int}.tsv"

    coord_df = pd.DataFrame({
        'embedding_id': embedding_ids,
        'x': coordinates[:, 0],
        'y': coordinates[:, 1]
    })
    coord_df.to_csv(coord_file, sep='\\t', index=False)
            
    connections = []
    for contig_id, indices in groups.items():
        if len(indices) > 1:
            for i in range(len(indices) - 1):
                for j in range(i + 1, len(indices)):
                    id1 = embedding_ids[indices[i]]
                    id2 = embedding_ids[indices[j]]
                    connections.append([id1, id2])
    # Save connections as TSV
    connections_df = pd.DataFrame(connections, columns=['id1', 'id2'])
    connections_file = f"${publish_subdir}/connections.tsv"
    connections_df.to_csv(connections_file, sep='\\t', index=False)

print("Writing versions.yml")
import platform
from importlib import import_module

versions = [
    f'"${task.process}":',
    f'    python: {platform.python_version()}'
]

for package_name in ['numpy', 'pandas', 'torch', 'umap']:
    try:
        module = import_module(package_name)
        version = getattr(module, '__version__', None)
        if version:
            versions.append(f'    {package_name}: {version}')
    except Exception:
        pass

with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write('\\n'.join(versions) + '\\n')
    """
}