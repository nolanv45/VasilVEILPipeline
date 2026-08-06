process MODIFY_CLUSTERS {
    label "process_single"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas:2a5ca2e7dd4ced9c' :
        'community.wave.seqera.io/library/python_pandas:eb68d9296e3f036a' }"
        
    publishDir "${params.outdir}/05_final_analysis/hdbscan/clusters",
        mode: 'copy'
    
    input:
        path cluster_csv_dir // Directory containing the cluster CSV files
        
    output:
        path "hdbscan_modified_cluster_labels.tsv", emit: modified_clusters
        path "versions.yml", emit: versions

    script:
    """
#!/usr/bin/env python3
import os
import pandas as pd

nn = "${params.final_nn}"
md_tag = "${params.final_md}".replace('.', 'p')   # e.g. "0.5" -> "0p5" — must match CLUSTER_HDBSCAN's tag

file1 = os.path.join("clusters_csv", f"hdbscan_nn{nn}_md{md_tag}_minclust30_clusters.csv")
file2 = os.path.join("clusters_csv", f"hdbscan_nn{nn}_md{md_tag}_minclust10_clusters.csv")

df1 = pd.read_csv(file1, sep=",", header=0)
df2 = pd.read_csv(file2, sep=",", header=0)

max_cluster_label = df1['cluster_label'].max()
next_cluster = max_cluster_label + 1

df1['label_source'] = 'minclust30'

df1_unclustered = df1[df1['cluster_label'] == -1]
embedding_ids_unclustered = df1_unclustered['embedding_id']
matching_rows = df2[df2['embedding_id'].isin(embedding_ids_unclustered)]
unmatched_rows = df2[~df2['embedding_id'].isin(embedding_ids_unclustered)]
matching_rows = matching_rows.sort_values(by='cluster_label')

for cluster_label in matching_rows['cluster_label'].unique():
    if cluster_label > -1:
        if cluster_label not in unmatched_rows['cluster_label'].values:
            cluster_mask = matching_rows['cluster_label'] == cluster_label
            embedding_ids = matching_rows[cluster_mask]['embedding_id']
            for emb_id in embedding_ids:
                df1.loc[df1['embedding_id'] == emb_id, 'cluster_label'] = next_cluster
                df1.loc[df1['embedding_id'] == emb_id, 'label_source'] = f"minclust10_{cluster_label}"
            next_cluster += 1

df1['label_source'] = df1['label_source'].fillna('minclust30')
df1_sorted = df1.sort_values(by='cluster_label')
df1_sorted.to_csv("hdbscan_modified_cluster_labels.tsv", index=False, sep="\\t")

import platform
process_name = "${task.process}"
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"{process_name}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    """
}