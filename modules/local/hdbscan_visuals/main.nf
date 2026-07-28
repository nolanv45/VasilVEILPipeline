process HDBSCAN_VISUALS {
    label "process_single"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib:49fe58d99685085a' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib:81ea56d306702341' }"

    publishDir "${params.outdir}/05_final_analysis/hdbscan/visuals",
        mode: 'copy'

    input:
        path cluster_info_tsv
        path metadata_file

    output:
        path "hdbscan_cluster_composition_heatmap.png", emit: composition_heatmap
        path "hdbscan_dataset_composition_stacked_bar.png", emit: dataset_stacked
        path "hdbscan_genofeature_composition_stacked_bar.png", emit: genofeature_stacked
        path "hdbscan_genofeature_composition_size_scaled_stacked_bar.png", emit: genofeature_size_scaled
        path "versions.yml", emit: versions

    script:
    """
#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


def save_placeholder(output_file, title, message):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axis('off')
    ax.set_title(title, fontsize=12)
    ax.text(0.5, 0.5, message, ha='center', va='center', fontsize=10, wrap=True)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


cluster_info = pd.read_csv("${cluster_info_tsv}", sep='\t')

metadata_df = pd.read_csv("${metadata_file}", sep='\t')
if 'genofeature' in metadata_df.columns and 'color' in metadata_df.columns:
    metadata_color_map = dict(zip(metadata_df['genofeature'], metadata_df['color']))
else:
    metadata_color_map = {}

if cluster_info.empty:
    save_placeholder(
        "hdbscan_cluster_composition_heatmap.png",
        "Cluster Composition Heatmap",
        "No clustered points available (all points are unclustered, cluster_label = -1)."
    )
    save_placeholder(
        "hdbscan_dataset_composition_stacked_bar.png",
        "Dataset Composition (Stacked)",
        "No clustered points available (all points are unclustered, cluster_label = -1)."
    )
    save_placeholder(
        "hdbscan_genofeature_composition_stacked_bar.png",
        "Genofeature Composition (Stacked)",
        "No clustered points available (all points are unclustered, cluster_label = -1)."
    )
    save_placeholder(
        "hdbscan_genofeature_composition_size_scaled_stacked_bar.png",
        "Genofeature Composition (Size-Scaled)",
        "No clustered points available (all points are unclustered, cluster_label = -1)."
    )
    raise SystemExit(0)


cluster_info = cluster_info.sort_values(by='cluster_id').reset_index(drop=True)
cluster_labels = cluster_info['cluster_id'].astype(str).tolist()


# -----------------------------------------------------
# Plot 1: Cluster composition heatmap (genofeature counts)
# -----------------------------------------------------
genofeature_cols = [
    col for col in cluster_info.columns
    if col.endswith('_count') and not col.startswith('dataset__')
]

if genofeature_cols:
    heatmap_data = cluster_info[genofeature_cols].to_numpy(dtype=float)
    row_sums = heatmap_data.sum(axis=1, keepdims=True)
    heatmap_prop = np.divide(
        heatmap_data,
        row_sums,
        out=np.zeros_like(heatmap_data),
        where=row_sums > 0
    )

    fig_w = max(8, min(24, 0.45 * len(genofeature_cols) + 4))
    fig_h = max(4, min(20, 0.35 * len(cluster_labels) + 3))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(heatmap_prop, aspect='auto', cmap='viridis')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Within-cluster proportion')

    ax.set_title('Cluster Composition Heatmap (Genofeature Proportions)')
    ax.set_xlabel('Genofeature')
    ax.set_ylabel('Cluster ID')
    ax.set_xticks(np.arange(len(genofeature_cols)))
    ax.set_xticklabels([col[:-6] for col in genofeature_cols], rotation=90)
    ax.set_yticks(np.arange(len(cluster_labels)))
    ax.set_yticklabels(cluster_labels)
    plt.tight_layout()
    plt.savefig('hdbscan_cluster_composition_heatmap.png', dpi=300)
    plt.close()
else:
    save_placeholder(
        'hdbscan_cluster_composition_heatmap.png',
        'Cluster Composition Heatmap',
        'No genofeature count columns found in cluster info TSV.'
    )


# -----------------------------------------------------
# Plot 2: Dataset composition stacked bar chart
# -----------------------------------------------------
dataset_cols = [
    col for col in cluster_info.columns
    if col.startswith('dataset__') and col.endswith('_count')
]

if dataset_cols:
    dataset_names = [col[len('dataset__'):-len('_count')] for col in dataset_cols]
    dataset_counts = cluster_info[dataset_cols].to_numpy(dtype=float)
    row_sums = dataset_counts.sum(axis=1, keepdims=True)
    dataset_prop = np.divide(
        dataset_counts,
        row_sums,
        out=np.zeros_like(dataset_counts),
        where=row_sums > 0
    )

    fig_w = max(9, min(24, 0.6 * len(cluster_labels) + 4))
    fig, ax = plt.subplots(figsize=(fig_w, 6))
    bottom = np.zeros(len(cluster_labels))
    colors = plt.cm.tab20(np.linspace(0, 1, max(1, len(dataset_names))))

    for i, dataset in enumerate(dataset_names):
        values = dataset_prop[:, i]
        ax.bar(cluster_labels, values, bottom=bottom, width=0.8, color=colors[i], label=dataset)
        bottom += values

    ax.set_title('Dataset Composition by Cluster (Stacked Proportions)')
    ax.set_xlabel('Cluster ID')
    ax.set_ylabel('Proportion within cluster')
    ax.set_ylim(0, 1)
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title='Dataset', bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
    plt.tight_layout()
    plt.savefig('hdbscan_dataset_composition_stacked_bar.png', dpi=300)
    plt.close()
else:
    save_placeholder(
        'hdbscan_dataset_composition_stacked_bar.png',
        'Dataset Composition (Stacked)',
        'No dataset count columns found in cluster info TSV.'
    )


# -----------------------------------------------------
# Plot 3: Genofeature composition stacked bar chart
# -----------------------------------------------------
if genofeature_cols:
    genofeature_counts = cluster_info[genofeature_cols].to_numpy(dtype=float)
    row_sums = genofeature_counts.sum(axis=1, keepdims=True)
    genofeature_prop = np.divide(
        genofeature_counts,
        row_sums,
        out=np.zeros_like(genofeature_counts),
        where=row_sums > 0
    )

    genofeature_names = [col[:-6] for col in genofeature_cols]
    fig_w = max(9, min(24, 0.6 * len(cluster_labels) + 4))
    fig, ax = plt.subplots(figsize=(fig_w, 6))
    bottom = np.zeros(len(cluster_labels))
    fallback_colors = plt.cm.tab20(np.linspace(0, 1, max(1, len(genofeature_names))))

    for i, genofeature in enumerate(genofeature_names):
        values = genofeature_prop[:, i]
        color = metadata_color_map.get(genofeature, fallback_colors[i])
        ax.bar(cluster_labels, values, bottom=bottom, width=0.8, color=color, label=genofeature)
        bottom += values

    ax.set_title('Genofeature Composition by Cluster (Stacked Proportions)')
    ax.set_xlabel('Cluster ID')
    ax.set_ylabel('Proportion within cluster')
    ax.set_ylim(0, 1)
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title='Genofeature', bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
    plt.tight_layout()
    plt.savefig('hdbscan_genofeature_composition_stacked_bar.png', dpi=300)
    plt.close()
else:
    save_placeholder(
        'hdbscan_genofeature_composition_stacked_bar.png',
        'Genofeature Composition (Stacked)',
        'No genofeature count columns found in cluster info TSV.'
    )


# -----------------------------------------------------
# Plot 4: Genofeature composition stacked by cluster-size proportion
# -----------------------------------------------------
if genofeature_cols:
    genofeature_counts = cluster_info[genofeature_cols].to_numpy(dtype=float)
    cluster_sizes = genofeature_counts.sum(axis=1)
    total_clustered = cluster_sizes.sum()

    if total_clustered > 0:
        genofeature_global_prop = genofeature_counts / total_clustered
    else:
        genofeature_global_prop = np.zeros_like(genofeature_counts)

    genofeature_names = [col[:-6] for col in genofeature_cols]
    fig_w = max(9, min(24, 0.6 * len(cluster_labels) + 4))
    fig, ax = plt.subplots(figsize=(fig_w, 6))
    bottom = np.zeros(len(cluster_labels))
    fallback_colors = plt.cm.tab20(np.linspace(0, 1, max(1, len(genofeature_names))))

    for i, genofeature in enumerate(genofeature_names):
        values = genofeature_global_prop[:, i]
        color = metadata_color_map.get(genofeature, fallback_colors[i])
        ax.bar(cluster_labels, values, bottom=bottom, width=0.8, color=color, label=genofeature)
        bottom += values

    ax.set_title('Genofeature Composition by Cluster (Height Scaled by Cluster Size)')
    ax.set_xlabel('Cluster ID')
    ax.set_ylabel('Proportion of all clustered members')
    ax.set_ylim(0, max(0.05, bottom.max() * 1.1 if len(bottom) else 0.05))
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title='Genofeature', bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
    plt.tight_layout()
    plt.savefig('hdbscan_genofeature_composition_size_scaled_stacked_bar.png', dpi=300)
    plt.close()
else:
    save_placeholder(
        'hdbscan_genofeature_composition_size_scaled_stacked_bar.png',
        'Genofeature Composition (Size-Scaled)',
        'No genofeature count columns found in cluster info TSV.'
    )

import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    numpy: {np.__version__}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    handle.write(f'    matplotlib: {matplotlib.__version__}\\n')
"""
}