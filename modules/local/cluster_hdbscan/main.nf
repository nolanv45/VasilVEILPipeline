process CLUSTER_HDBSCAN {
    tag "nn${nn}_md${md_tag}_mc${mc}"
    label "process_medium"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pytorch_numpy_pandas_pruned:81ab5dec161369ec' :
        'community.wave.seqera.io/library/python_pytorch_numpy_pandas_pruned:6c2364c205cfbeb9' }"
  
    publishDir "${params.outdir}/04_parameter_selection/hdbscan/clusters_csv",
        mode: 'copy',
        pattern: "*.csv"
        
    input:
        tuple val(nn), val(md_tag), path(coord_file), val(mc)
        val embeddings_dirs
        path filtered_tsv
        path metadata_file
        
    output:        
        path "hdbscan_nn${nn}_md${md_tag}_minclust${mc}.png", emit: plot
        path "hdbscan_nn${nn}_md${md_tag}_minclust${mc}_clusters.csv", emit: clusters_csv

    script:
    """
#!/usr/bin/env python3
import os
import re
os.environ.setdefault('MPLCONFIGDIR', '/tmp/matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import hdbscan
import torch

nn = ${nn}
md = float("${md_tag}".replace("p", "."))
mc = ${mc}

module_df = pd.read_csv("${filtered_tsv}", sep='\\t')
metadata_df = pd.read_csv("${metadata_file}", sep='\\t')

def find_pt_files(base_dir):
    pt_files = []
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith('.pt'):
                pt_files.append(os.path.join(root, file))
    return pt_files

def normalize_input_dirs(raw_value):
    raw_value = str(raw_value).strip()
    if not raw_value:
        return []
    cleaned = raw_value.strip('[]')
    return [item for item in re.split(r'[\\s,]+', cleaned) if item]

def load_raw_embeddings(embeddings_dirs_str):
    embeddings = []
    embedding_ids = []
    base_dirs = normalize_input_dirs(embeddings_dirs_str) if isinstance(embeddings_dirs_str, str) else embeddings_dirs_str
    pt_files = []
    for base_dir in base_dirs:
        pt_files.extend(find_pt_files(base_dir))
    print(f"Found {len(pt_files)} .pt files")
    for file_path in pt_files:
        try:
            data = torch.load(file_path, map_location='cpu')
            emb = data["mean_representations"][36].numpy()
            eid = data.get("label")
            if emb is not None and eid is not None:
                embeddings.append(emb)
                embedding_ids.append(eid)
        except Exception as e:
            print(f"Warning: could not load {file_path}: {e}")
    return np.stack(embeddings), embedding_ids

raw_embeddings, raw_embedding_ids = load_raw_embeddings("${embeddings_dirs}")

coord_file = "${coord_file}"
output_path = "hdbscan_nn${nn}_md${md_tag}_minclust${mc}.png"
cluster_output_path = "hdbscan_nn${nn}_md${md_tag}_minclust${mc}_clusters.csv"

try:
    coord_df = pd.read_csv(coord_file, sep='\\t')
    embedding_ids_2d = coord_df['embedding_id'].tolist()
    embedding_2d = coord_df[['x', 'y']].values

    id_to_raw = {eid: emb for eid, emb in zip(raw_embedding_ids, raw_embeddings)}
    aligned_embeddings = np.stack([id_to_raw[eid] for eid in embedding_ids_2d if eid in id_to_raw])
    aligned_ids = [eid for eid in embedding_ids_2d if eid in id_to_raw]
    aligned_2d = np.stack([embedding_2d[i] for i, eid in enumerate(embedding_ids_2d) if eid in id_to_raw])

    print(f"Clustering on {aligned_embeddings.shape[1]}-dimensional embeddings, plotting on 2D coords")

    labels = hdbscan.HDBSCAN(min_samples=2, min_cluster_size=mc).fit_predict(aligned_embeddings)
    clustered = (labels >= 0)

    module_df["normalized_orf_id"] = module_df["orf_id"].str.replace(".", "-", regex=False)

    cluster_df = pd.DataFrame({
        'embedding_id': aligned_ids,
        'cluster_label': labels,
        'genofeature': [module_df.loc[module_df["normalized_orf_id"] == eid, "genofeature"].iloc[0]
            if eid in module_df["normalized_orf_id"].values else "unknown"
            for eid in aligned_ids]
    })
    cluster_df.to_csv(cluster_output_path, index=False)

    fig, ax = plt.subplots(figsize=(10, 10))
    scatter = ax.scatter(aligned_2d[clustered, 0], aligned_2d[clustered, 1],
                         c=labels[clustered], s=10, alpha=0.7, cmap="Spectral")
    ax.scatter(aligned_2d[~clustered, 0], aligned_2d[~clustered, 1],
               color="gray", s=10, alpha=0.5, label="not clustered")

    unique_labels = np.unique(labels[labels >= 0])
    for cluster in unique_labels:
        cluster_points = aligned_2d[labels == cluster]
        cx, cy = cluster_points[:, 0].mean(), cluster_points[:, 1].mean()
        plt.text(cx, cy, str(cluster), fontsize=10, ha='center', va='center')

    ax.set_xticks([])
    ax.set_yticks([])
    plt.grid(False)
    plt.title(f"UMAP (nn={nn}, md={md}) with HDBSCAN (minclust={mc})")
    plt.colorbar(scatter, label="Cluster Labels")
    plt.tight_layout()
    plt.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close()

except Exception as e:
    print(f"Error processing plot: {e}")
    raise
    """
}
