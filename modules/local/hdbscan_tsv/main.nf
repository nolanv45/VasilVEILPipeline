process HDBSCAN_TSV {
    label "process_single"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/05_final_analysis/hdbscan",
        mode: 'copy'
    
    input:
        path modified_clusters  // TSV from MODIFY_CLUSTERS with embedding_id, cluster_label, genofeature, label_source
        path contig_file  // Path to contig_df.tsv or combined_datasets.tsv
        
    output:
        path "hdbscan_cluster_info.tsv", emit: cluster_info
        path "hdbscan_cluster_relationships.tsv", emit: cluster_relationships
        path "versions.yml", emit: versions

    script:
    """
#!/usr/bin/env python3
import pandas as pd
from itertools import combinations

# Load data
clusters_df = pd.read_csv("${modified_clusters}", sep="\t")
contig_df = pd.read_csv("${contig_file}", sep="\t")

print("Clusters DataFrame shape:", clusters_df.shape)
print("Contig DataFrame shape:", contig_df.shape)

# The clusters_df already has embedding_id, cluster_label, genofeature, label_source
# We need to add dataset and contig_id information

# If contig_df is combined_datasets.tsv, it has orf_id = embedding_id
if 'orf_id' in contig_df.columns:
    # This is combined_datasets.tsv
    clusters_df = clusters_df.merge(
        contig_df[['contig_id', 'orf_id', 'dataset']].drop_duplicates(),
        left_on='embedding_id',
        right_on='orf_id',
        how='left'
    )
else:
    # Fallback: try to extract dataset and contig_id from embedding_id pattern
    # embedding_id format: DATASET_NODE_..._position_..._...
    def extract_dataset_from_embedding(emb_id):
        parts = emb_id.split('_')
        return parts[0]
    
    clusters_df['dataset'] = clusters_df['embedding_id'].apply(extract_dataset_from_embedding)

print("After merging with contig info:")
print(clusters_df.head())

# Filter out unclustered (-1) for most analyses
clustered_df = clusters_df[clusters_df['cluster_label'] >= 0].copy()

# Get unique clusters (excluding -1)
unique_clusters = sorted(clustered_df['cluster_label'].unique())
print(f"\\nFound {len(unique_clusters)} clusters")

# =====================================================
# TABLE 1: CLUSTER INFO (one row per cluster)
# =====================================================
print("\\n=== Generating Cluster Info Table ===")

# First: Genofeature counts for pivoting
genofeature_count = clustered_df.groupby(['cluster_label', 'genofeature']).size().reset_index(name='count')
genofeature_count_pivot = genofeature_count.pivot_table(
    index='cluster_label',
    columns='genofeature',
    values='count',
    fill_value=0
)

# Dataset counts for stacked composition plotting
dataset_count = clustered_df.groupby(['cluster_label', 'dataset']).size().reset_index(name='count')
dataset_count_pivot = dataset_count.pivot_table(
    index='cluster_label',
    columns='dataset',
    values='count',
    fill_value=0
)

# Build cluster info
cluster_info_list = []

for cluster in unique_clusters:
    cluster_data = clustered_df[clustered_df['cluster_label'] == cluster]
    
    num_embeddings = len(cluster_data)
    num_contigs = len(cluster_data['contig_id'].dropna().unique())
    num_genofeatures = len(cluster_data['genofeature'].unique())
    datasets_in_cluster = sorted(cluster_data['dataset'].dropna().unique())
    num_datasets = len(datasets_in_cluster)
    
    # Most common genofeatures
    top_genofeatures = cluster_data['genofeature'].value_counts().to_dict()
    top_gf_str = ', '.join([f"{gf}({cnt})" for gf, cnt in sorted(top_genofeatures.items(), key=lambda x: x[1], reverse=True)])
    
    row = {
        'cluster_id': cluster,
        'num_embeddings': num_embeddings,
        'num_contigs': num_contigs,
        'num_genofeatures': num_genofeatures,
        'num_datasets': num_datasets,
        'datasets': '|'.join(datasets_in_cluster),
        'top_genofeatures': top_gf_str
    }
    
    # Add genofeature counts as columns
    for gf in genofeature_count_pivot.columns:
        row[f"{gf}_count"] = int(genofeature_count_pivot.loc[cluster, gf]) if cluster in genofeature_count_pivot.index else 0

    # Add dataset counts as columns
    for dataset in dataset_count_pivot.columns:
        row[f"dataset__{dataset}_count"] = int(dataset_count_pivot.loc[cluster, dataset]) if cluster in dataset_count_pivot.index else 0
    
    cluster_info_list.append(row)

cluster_info_columns = [
    'cluster_id',
    'num_embeddings',
    'num_contigs',
    'num_genofeatures',
    'num_datasets',
    'datasets',
    'top_genofeatures'
]
cluster_info_columns += [f"{gf}_count" for gf in genofeature_count_pivot.columns]
cluster_info_columns += [f"dataset__{dataset}_count" for dataset in dataset_count_pivot.columns]

cluster_info_df = pd.DataFrame(cluster_info_list, columns=cluster_info_columns)
cluster_info_df.to_csv("hdbscan_cluster_info.tsv", sep="\t", index=False)
print("Saved: hdbscan_cluster_info.tsv")
print(f"Shape: {cluster_info_df.shape}")

# =====================================================
# TABLE 2: CLUSTER RELATIONSHIPS (one row per cluster pair)
# =====================================================
print("\\n=== Generating Cluster Relationships Table ===")

relationships_data = []

for cluster1, cluster2 in combinations(unique_clusters, 2):
    # Get contigs in each cluster
    contigs_c1 = set(clustered_df[clustered_df['cluster_label'] == cluster1]['contig_id'].dropna())
    contigs_c2 = set(clustered_df[clustered_df['cluster_label'] == cluster2]['contig_id'].dropna())
    
    # Overlap
    overlap_contigs = len(contigs_c1 & contigs_c2)
    
    # Proportion relative to each cluster (asymmetric)
    prop_c1_in_c2 = overlap_contigs / len(contigs_c1) if len(contigs_c1) > 0 else 0
    prop_c2_in_c1 = overlap_contigs / len(contigs_c2) if len(contigs_c2) > 0 else 0
    
    # Symmetric overlap proportion (balanced view)
    union_contigs = len(contigs_c1 | contigs_c2)
    symmetric_overlap = overlap_contigs / union_contigs if union_contigs > 0 else 0
    
    # Jaccard index (intersection / union)
    jaccard = overlap_contigs / union_contigs if union_contigs > 0 else 0
    
    relationships_data.append({
        'cluster_1': cluster1,
        'cluster_2': cluster2,
        'shared_contigs': overlap_contigs,
        'contigs_in_cluster_1': len(contigs_c1),
        'contigs_in_cluster_2': len(contigs_c2),
        'total_unique_contigs': union_contigs,
        'proportion_c1_in_c2': round(prop_c1_in_c2, 4),
        'proportion_c2_in_c1': round(prop_c2_in_c1, 4),
        'symmetric_overlap': round(symmetric_overlap, 4),
        'jaccard_index': round(jaccard, 4)
    })

relationships_columns = [
    'cluster_1',
    'cluster_2',
    'shared_contigs',
    'contigs_in_cluster_1',
    'contigs_in_cluster_2',
    'total_unique_contigs',
    'proportion_c1_in_c2',
    'proportion_c2_in_c1',
    'symmetric_overlap',
    'jaccard_index'
]

relationships_df = pd.DataFrame(relationships_data, columns=relationships_columns)
relationships_df.to_csv("hdbscan_cluster_relationships.tsv", sep="\t", index=False)
print("Saved: hdbscan_cluster_relationships.tsv")
print(f"Shape: {relationships_df.shape}")

print("\\n=== HDBSCAN_TSV Complete ===")
print(f"Total clusters: {len(unique_clusters)}")
print(f"Total unclustered embeddings: {len(clusters_df[clusters_df['cluster_label'] == -1])}")

import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
"""
}

