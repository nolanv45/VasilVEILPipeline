process GENOFEATURE_CENTRIC {
    label "process_medium"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib_pillow:3f1442a23078f4ef' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib_pillow:e3f9f64a94d792de' }"

    publishDir "${params.outdir}/05_final_analysis/umap",
        mode: 'copy'
    
    input:
        path module_file  // Path to the module file with genofeature information
        path coordinates_file
        path connections_file  // Path to the connections file
        path metadata_file  // Path to the metadata file

        
    output:
        path "genofeature_plots/", emit: plots
        path "versions.yml", emit: versions


    script:
    """
#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from collections import defaultdict
from PIL import Image
import math

# Load data
coords_file = "${coordinates_file}"  # Path to coordinates file
module_file = "${module_file}"  # Path to module file with genofeature information
metadata_file = "${metadata_file}"  # Path to metadata file
connections_file = "${connections_file}"

coord_df = pd.read_csv(coords_file, sep='\\t')
embedding_id_to_index = {eid: idx for idx, eid in enumerate(coord_df['embedding_id'])}
embedding_2d = coord_df[['x', 'y']].values

module_df = pd.read_csv(module_file, sep='\\t')
metadata_df = pd.read_csv(metadata_file, sep="\\t")
color_map = dict(zip(metadata_df["genofeature"], metadata_df["color"]))
display_map = dict(zip(metadata_df["genofeature"], metadata_df["display_name"]))
marker_map = dict(zip(metadata_df["genofeature"], metadata_df["marker"]))
genofeature_cols = [col for col in module_df.columns if col not in ['dataset', 'contig_id', 'module']]

connections = []
if os.path.exists(connections_file):
    connections_df = pd.read_csv(connections_file, sep='\\t')
    connections = list(zip(connections_df['id1'], connections_df['id2']))

plot_output_dir = "genofeature_plots"
os.makedirs(plot_output_dir, exist_ok=True)

excluded = "${params.excluded_genofeatures.join(',')}".split(',') if "${params.excluded_genofeatures.join(',')}" else []
visible_features_list = [f for f in genofeature_cols if f not in excluded]

for feature in visible_features_list:
    # 1. Find all contigs where any visible_feature is in the module (split by '_')
    contigs = module_df[module_df['module'].str.split('_').apply(lambda x: feature in x)]['contig_id'].unique()
    if len(contigs) == 0:
        print(f"Skipping {feature}: not present in data.")
        continue

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

    fig, ax = plt.subplots(figsize=(10, 7))

    # Plot connections
    for id1, id2 in connections:
        if id1 in all_eids and id2 in all_eids:
            idx1 = embedding_id_to_index.get(id1)
            idx2 = embedding_id_to_index.get(id2)
            if idx1 is not None and idx2 is not None:
                ax.plot([embedding_2d[idx1, 0], embedding_2d[idx2, 0]],
                        [embedding_2d[idx1, 1], embedding_2d[idx2, 1]],
                        color="gray", alpha=0.2, linewidth=0.5, zorder=0)

    # Plot all points as faded background
    for eid, idx in embedding_id_to_index.items():
        if eid not in context_set and eid not in highlight_set:
            x, y = embedding_2d[idx]
            ax.scatter(x, y, c="#808080", marker="o", s=10, alpha=0.2, zorder=0)


    # Plot context points with metadata-driven style
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

    # Plot highlighted points with metadata-driven style
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

    # Add visible_feature label below the plot
    label = feature
    plt.figtext(0.5, 0.01, label, ha='center', va='bottom', fontsize=20)


    feature_set = set(gf for _, gf, _ in embedding_info)
    handles = [
        mlines.Line2D([], [], color=color_map.get(f, "#808080"), marker=marker_map.get(f, "."), linestyle="None", markersize=8, label=display_map.get(f, f))
        for f in metadata_df["genofeature"] if f in feature_set
    ]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=10, title="Genofeatures", title_fontsize=12)

    plt.savefig(os.path.join(plot_output_dir, f"umap_{feature}_recolor_mqc.png"), dpi=600, bbox_inches="tight")
    plt.close()

    # --------- PER-DATASET PLOTS ---------
    # Build mapping from embedding_id to dataset
    embedding_to_dataset = {}
    for _, row in module_df.iterrows():
        for gf in genofeature_cols:
            eid = row[gf]
            if pd.notnull(eid) and str(eid).strip():
                embedding_to_dataset[eid] = row['dataset']

    # Get all unique datasets for this feature
    datasets_for_feature = sorted({embedding_to_dataset[eid] for eid in all_eids if eid in embedding_to_dataset})

    for dataset in datasets_for_feature:
        fig, ax = plt.subplots(figsize=(10, 10))

        # Plot connections (same as above)
        for id1, id2 in connections:
            if id1 in all_eids and id2 in all_eids:
                idx1 = embedding_id_to_index.get(id1)
                idx2 = embedding_id_to_index.get(id2)
                if idx1 is not None and idx2 is not None:
                    ax.plot([embedding_2d[idx1, 0], embedding_2d[idx2, 0]],
                            [embedding_2d[idx1, 1], embedding_2d[idx2, 1]],
                            color="gray", alpha=0.2, linewidth=0.5, zorder=0)

        # Plot all points, using dataset and context/highlight logic
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
                    # In selected dataset, but not context or highlight
                    ax.scatter(x, y, c="#808080", marker="o", s=10, alpha=0.2, zorder=0)
        ax.set_xticks([])
        ax.set_yticks([])
        label = f"{feature} | Dataset: {dataset}"
        plt.figtext(0.5, 0.01, label, ha='center', va='bottom', fontsize=20)
        plt.savefig(os.path.join(plot_output_dir, f"umap_{feature}_by_dataset_{dataset}_mqcS.png"), dpi=600, bbox_inches="tight")
        plt.close()

# --------- AUTOMATIC TILING ---------
def tile_images(image_dir, output_file, padding=10):
    png_files = [os.path.join(image_dir, f) for f in os.listdir(image_dir) 
                 if f.lower().endswith("_mqc.png")]
    png_files.sort()
    if not png_files:
        print("No plot PNG files found for tiling.")
        return
    print(f"Tiling {len(png_files)} plot images...")
    images = [Image.open(f) for f in png_files]
    heights = [img.height for img in images]
    target_height = max(heights)
    resized_images = []
    resized_sizes = []
    for img in images:
        w, h = img.size
        new_w = int(w * (target_height / h))
        resized_img = img.resize((new_w, target_height), Image.LANCZOS)
        resized_images.append(resized_img)
        resized_sizes.append((new_w, target_height))
    n = len(resized_images)
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)
    col_widths = [0] * cols
    row_heights = [target_height] * rows
    for idx, (w, h) in enumerate(resized_sizes):
        row = idx // cols
        col = idx % cols
        if w > col_widths[col]:
            col_widths[col] = w
    x_offsets = [0]
    for w in col_widths[:-1]:
        x_offsets.append(x_offsets[-1] + w + padding)
    y_offsets = [0]
    for h in row_heights[:-1]:
        y_offsets.append(y_offsets[-1] + h + padding)
    canvas_width = sum(col_widths) + padding * (cols - 1)
    canvas_height = sum(row_heights) + padding * (rows - 1)
    canvas = Image.new('RGBA', (canvas_width, canvas_height), (255, 255, 255, 255))
    for idx, img in enumerate(resized_images):
        row = idx // cols
        col = idx % cols
        x = x_offsets[col]
        y = y_offsets[row]
        canvas.paste(img, (x, y))
    canvas.save(output_file)
    print(f"Tiled image saved as {output_file}")

# Create tiled image of all plots
tiled_output_file = os.path.join(plot_output_dir, "genofeature_highlighting_tiled_all_mqc.png")
tile_images(plot_output_dir, tiled_output_file)
print(f"All plots completed!")

import platform
import PIL

with open("versions.yml", "w", encoding="utf-8") as versions_handle:
    versions_handle.write('"${task.process}":\\n')
    versions_handle.write(f'    python: {platform.python_version()}\\n')
    versions_handle.write(f'    numpy: {np.__version__}\\n')
    versions_handle.write(f'    pandas: {pd.__version__}\\n')
    versions_handle.write(f'    matplotlib: {plt.matplotlib.__version__}\\n')
    versions_handle.write(f'    pillow: {PIL.__version__}\\n')
    """
}