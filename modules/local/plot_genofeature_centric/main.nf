process PLOT_GENOFEATURE_CENTRIC {
    label "process_medium"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib_pillow:3f1442a23078f4ef' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib_pillow:e3f9f64a94d792de' }"

    publishDir "${params.outdir}/05_final_analysis/umap/genofeature_plots",
        mode: 'copy',
        pattern: "*.png"
    
    input:
        val feature
        path module_file  // Path to the module file with genofeature information
        path coordinates_file
        path connections_file  // Path to the connections file
        path metadata_file  // Path to the metadata file

        
    output:
        path "umap_${feature}_recolor_mqc.png", emit: overview, optional: true
        path "umap_${feature}_by_dataset_*_mqcS.png", emit: by_dataset, optional: true


    script:
    """
#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

feature = "${feature}"

coord_df = pd.read_csv("${coordinates_file}", sep='\\t')
embedding_id_to_index = {eid: idx for idx, eid in enumerate(coord_df['embedding_id'])}
embedding_2d = coord_df[['x', 'y']].values

module_df = pd.read_csv("${module_file}", sep='\\t')
metadata_df = pd.read_csv("${metadata_file}", sep='\\t')
color_map = dict(zip(metadata_df["genofeature"], metadata_df["color"]))
display_map = dict(zip(metadata_df["genofeature"], metadata_df["display_name"]))
marker_map = dict(zip(metadata_df["genofeature"], metadata_df["marker"]))
genofeature_cols = [col for col in module_df.columns if col not in ['dataset', 'contig_id', 'module']]

connections = []
if os.path.exists("${connections_file}"):
    connections_df = pd.read_csv("${connections_file}", sep='\\t')
    connections = list(zip(connections_df['id1'], connections_df['id2']))

# 1. Find all contigs where this feature is in the module (split by '_')
contigs = module_df[module_df['module'].str.split('_').apply(lambda x: feature in x)]['contig_id'].unique()

if len(contigs) == 0:
    print(f"Skipping {feature}: not present in data.")
else:
    embedding_info = []
    for contig in contigs:
        rows = module_df[module_df['contig_id'] == contig]
        for _, row in rows.iterrows():
            for gf in genofeature_cols:
                eid = row[gf]
                if pd.notnull(eid) and str(eid).strip():
                    embedding_info.append((eid, gf, contig))

    highlight_set = {eid for eid, gf, _ in embedding_info if gf == feature}
    context_set = {eid for eid, gf, _ in embedding_info if gf != feature}
    all_eids = highlight_set | context_set

    fig, ax = plt.subplots(figsize=(10, 10))

    for id1, id2 in connections:
        if id1 in all_eids and id2 in all_eids:
            idx1 = embedding_id_to_index.get(id1)
            idx2 = embedding_id_to_index.get(id2)
            if idx1 is not None and idx2 is not None:
                ax.plot([embedding_2d[idx1, 0], embedding_2d[idx2, 0]],
                        [embedding_2d[idx1, 1], embedding_2d[idx2, 1]],
                        color="gray", alpha=0.2, linewidth=0.5, zorder=0)

    for eid, idx in embedding_id_to_index.items():
        if eid not in context_set and eid not in highlight_set:
            x, y = embedding_2d[idx]
            ax.scatter(x, y, c="#808080", marker="o", s=10, alpha=0.2, zorder=0)

    for eid in context_set:
        idx = embedding_id_to_index.get(eid)
        if idx is not None:
            gf = next((gf for e, gf, _ in embedding_info if e == eid), "unknown")
            row = metadata_df[metadata_df["genofeature"] == gf]
            if not row.empty:
                row = row.iloc[0]
                color = row["color"]
                marker = row["marker"]
            else:
                color = "#808080"
                marker = "o"
            ax.scatter(
                embedding_2d[idx, 0], embedding_2d[idx, 1],
                c=color, marker=marker, s=10, alpha=0.8, zorder=1, linewidths=0.5
            )

    for eid in highlight_set:
        idx = embedding_id_to_index.get(eid)
        if idx is not None:
            row = metadata_df[metadata_df["genofeature"] == feature]
            if not row.empty:
                row = row.iloc[0]
                color = row["color"]
                marker = row["marker"]
            else:
                color = "#FF0000"
                marker = "o"
            ax.scatter(
                embedding_2d[idx, 0], embedding_2d[idx, 1],
                c=color, marker=marker, s=10, alpha=0.8, zorder=2, linewidths=0.5
            )

    ax.set_xticks([])
    ax.set_yticks([])
    plt.figtext(0.5, 0.01, feature, ha='center', va='bottom', fontsize=20)

    feature_set = set(gf for _, gf, _ in embedding_info)
    handles = [
        mlines.Line2D([], [], color=color_map.get(f, "#808080"), marker=marker_map.get(f, "."),
                      linestyle="None", markersize=8, label=display_map.get(f, f))
        for f in metadata_df["genofeature"] if f in feature_set
    ]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5),
              frameon=False, fontsize=10, title="Genofeatures", title_fontsize=12)

    plt.savefig(f"umap_{feature}_recolor_mqc.png", dpi=600, bbox_inches="tight")
    plt.close()

    # --------- PER-DATASET PLOTS ---------
    embedding_to_dataset = {}
    for _, row in module_df.iterrows():
        for gf in genofeature_cols:
            eid = row[gf]
            if pd.notnull(eid) and str(eid).strip():
                embedding_to_dataset[eid] = row['dataset']

    datasets_for_feature = sorted({embedding_to_dataset[eid] for eid in all_eids if eid in embedding_to_dataset})

    for dataset in datasets_for_feature:
        fig, ax = plt.subplots(figsize=(10, 10))

        for id1, id2 in connections:
            if id1 in all_eids and id2 in all_eids:
                idx1 = embedding_id_to_index.get(id1)
                idx2 = embedding_id_to_index.get(id2)
                if idx1 is not None and idx2 is not None:
                    ax.plot([embedding_2d[idx1, 0], embedding_2d[idx2, 0]],
                            [embedding_2d[idx1, 1], embedding_2d[idx2, 1]],
                            color="gray", alpha=0.2, linewidth=0.5, zorder=0)

        for eid, idx in embedding_id_to_index.items():
            ds = embedding_to_dataset.get(eid, None)
            x, y = embedding_2d[idx]
            if ds == dataset:
                if eid in highlight_set:
                    row = metadata_df[metadata_df["genofeature"] == feature]
                    if not row.empty:
                        row = row.iloc[0]
                        color = row["color"]
                        marker = row["marker"]
                    else:
                        color = "#FF0000"
                        marker = "o"
                    ax.scatter(x, y, c=color, marker=marker, s=10, alpha=0.8, zorder=2, linewidths=0.5)
                elif eid in context_set:
                    gf = next((gf for e, gf, _ in embedding_info if e == eid), "unknown")
                    row = metadata_df[metadata_df["genofeature"] == gf]
                    if not row.empty:
                        row = row.iloc[0]
                        color = row["color"]
                        marker = row["marker"]
                    else:
                        color = "#FF0000"
                        marker = "o"
                    ax.scatter(x, y, c=color, marker=marker, s=10, alpha=0.8, zorder=1, linewidths=0.5)
                else:
                    ax.scatter(x, y, c="#808080", marker="o", s=10, alpha=0.2, zorder=0)

        ax.set_xticks([])
        ax.set_yticks([])
        label = f"{feature} | Dataset: {dataset}"
        plt.figtext(0.5, 0.01, label, ha='center', va='bottom', fontsize=20)
        plt.savefig(f"umap_{feature}_by_dataset_{dataset}_mqcS.png", dpi=600, bbox_inches="tight")
        plt.close()
    """
}