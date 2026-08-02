process PLOT_UMAP {
    tag "nn${nn}_md${md_tag}"
    label "process_medium"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_pandas_matplotlib_pillow:bc19f6df598f92e8' :
        'community.wave.seqera.io/library/python_numpy_pandas_matplotlib_pillow:352bfae9d89bced3' }"
    
    input:
        tuple val(nn), val(md_tag), path(coord_file), path(conns_file)
        path filtered_tsv    // TSV file with metadata
        path metadata_file   // Metadata file with colors/markers
        
    output:
        path "umap_nn${nn}_md${md_tag}.png", emit: plot

    script:
    """
#!/usr/bin/env python3
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import pandas as pd

module_df = pd.read_csv("${filtered_tsv}", sep='\\t')
metadata_df = pd.read_csv("${metadata_file}", sep='\\t')

coord_file = "${coord_file}"
conns_file = "${conns_file}"
output_file = "umap_nn${nn}_md${md_tag}.png"

try:
    coord_df = pd.read_csv(coord_file, sep='\\t')
    embedding_ids = coord_df['embedding_id'].tolist()
    embedding_2d = coord_df[['x', 'y']].values

    print(f"Loaded coordinates from {coord_file}")
    print(f"Shape: {embedding_2d.shape}, IDs: {len(embedding_ids)}")

    module_df["normalized_orf_id"] = module_df["orf_id"].str.replace(".", "-", regex=False)

    colors = []
    markers = []
    for eid in embedding_ids:
        genofeature = module_df.loc[module_df["normalized_orf_id"] == eid, "genofeature"].iloc[0] if eid in module_df["normalized_orf_id"].values else None
        if genofeature is not None and genofeature in metadata_df["genofeature"].values:
            color = metadata_df.loc[metadata_df["genofeature"] == genofeature, "color"].iloc[0]
            marker = metadata_df.loc[metadata_df["genofeature"] == genofeature, "marker"].iloc[0]
        else:
            color = "#808080"
            marker = "."
        colors.append(color)
        markers.append(marker if pd.notna(marker) else ".")

    connections_df = pd.read_csv(conns_file, sep='\\t')
    if not connections_df.empty:
        print(f"Loaded {len(connections_df)} connections")
        id_to_idx = {eid: idx for idx, eid in enumerate(embedding_ids)}
        connections = []
        for _, row in connections_df.iterrows():
            if row['id1'] in id_to_idx and row['id2'] in id_to_idx:
                connections.append([id_to_idx[row['id1']], id_to_idx[row['id2']]])
        connections = np.array(connections)
    else:
        connections = []
        print("connections.tsv present but empty")

    fig, ax = plt.subplots(figsize=(10, 10))

    if len(connections) > 0:
        for conn in connections:
            x_coords = [embedding_2d[conn[0]][0], embedding_2d[conn[1]][0]]
            y_coords = [embedding_2d[conn[0]][1], embedding_2d[conn[1]][1]]
            ax.plot(x_coords, y_coords, color='#CCCCCC', alpha=0.5, linewidth=0.2, zorder=1)

    for coord, color, marker in zip(embedding_2d, colors, markers):
        ax.scatter(coord[0], coord[1], c=color, marker=marker, s=10, alpha=0.7)

    ax.set_xticks([])
    ax.set_yticks([])
    plt.grid(False)

    plt.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close()

except Exception as e:
    print(f"Error processing plot: {e}")
    print("Available columns in metadata_df:", metadata_df.columns.tolist())
    raise
    """
}