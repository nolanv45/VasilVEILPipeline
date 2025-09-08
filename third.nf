nextflow.enable.dsl=2

process GENERATE_COORDINATES_2 {
    publishDir "${params.outdir}/second_coordinates",
        mode: 'copy'
        
    label 'gpu'  
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    memory { 50.GB }
    
    input:
        path(embeddings)   
        
    output:
        path "second_coordinates/coordinates_nn*.tsv", emit: coordinates_tsv
        path "second_coordinates/connections.tsv", emit: connections_tsv
        
    script:
    def selected_list = params.selected_genofeatures_2.collect { "\'${it}\'" }.join(', ')
    """
    #!/usr/bin/env python3
    import os
    import torch
    import numpy as np
    import umap
    import random
    import pandas as pd
    from pathlib import Path

    # Create output directories
    os.makedirs("second_coordinates", exist_ok=True)
   
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
            if 'batch_0' in root:
                for file in files:
                    if file.endswith('.pt'):
                        pt_files.append(os.path.join(root, file))
        return pt_files

    # Load embeddings
    embeddings = []
    embedding_ids = []
    
    # Find all .pt files
    pt_files = find_pt_files("${embeddings}")
    print(f"Found {len(pt_files)} .pt files")
    
    # Process each .pt file
    selected_subdirs = [${selected_list}]
    for file_path in pt_files:
        genofeature = Path(file_path).parent.parent.name
        print(genofeature)
        if genofeature.lower() in [s.lower() for s in selected_subdirs]:
            print(f"Processing: {file_path}")
            embedding, embedding_id = load_embedding(file_path)
            if embedding is not None and embedding_id is not None:
                embeddings.append(embedding)
                embedding_ids.append(embedding_id)
                print(f"Successfully loaded embedding from {file_path}")

    if len(embeddings) > 0:
        print(f"Processing {len(embeddings)} embeddings")
        embeddings = np.stack(embeddings)

        groups = {}
        for i, eid in enumerate(embedding_ids):
            contig_id = '_'.join(eid.split('_')[:-3])  # Parse embedding ID
            if contig_id not in groups:
                groups[contig_id] = []
            groups[contig_id].append(i)
        
        # Generate UMAP for each parameter combination
        nn = ${params.final_nn}
        md = ${params.final_md}


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
        
        md_int = int(md * 10)
        coord_file = f"second_coordinates/coordinates_nn{nn}_md{md_int}.tsv"

        coord_df = pd.DataFrame({
            'embedding_id': embedding_ids,
            'x': coordinates[:, 0],
            'y': coordinates[:, 1]
        })
        coord_df.to_csv(coord_file, sep='\t', index=False)
                
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
        connections_file = f"second_coordinates/connections.tsv"
        connections_df.to_csv(connections_file, sep='\t', index=False)
    """
}

process MODULE_FILE {
    
    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}", 
        mode: 'copy',
        pattern: "*.tsv"

    input:
    path(tsv_file)
    path(metadata_file)

    output:
    path("contig_df.tsv"), emit: standardized

    script:
    def selected_list = params.selected_genofeatures_2.collect { "\'${it}\'" }.join(', ')
    """
#!/usr/bin/env python3
# will take tsv files that output from each dataset per orf basis.
# will filter genofeatures based on metadata file only, so user input required
import pandas as pd
import glob
import os

df = pd.read_csv("${tsv_file}", sep='\t')
metadata = pd.read_csv("${metadata_file}", sep='\t')
genofeatures = []
for feature in metadata['genofeature']:
    if feature in [${selected_list}]:
        genofeatures.append(feature)

df = df[df['genofeature'].isin(genofeatures)]
contig_df = pd.DataFrame({'contig_id': df['contig_id'].unique()})

for genofeature in genofeatures:
    feature_data = df[df['genofeature'] == genofeature]
    if len(feature_data) > 0:
        feature_orfs = feature_data.set_index('contig_id')['orf_id']
        contig_df[genofeature] = contig_df['contig_id'].map(
            feature_orfs.groupby('contig_id').first()
        )
    else:
        contig_df[genofeature] = None
genofeature_columns = [col for col in contig_df.columns if col != 'contig_id']
contig_df['module'] = contig_df[genofeature_columns].apply(
    lambda row: '_'.join([col for col, val in zip(genofeature_columns, row) if pd.notna(val)]), 
    axis=1
)
contig_df['dataset'] = contig_df['contig_id'].map(
    df.groupby('contig_id')['dataset'].first()
)
cols = ['dataset', 'contig_id'] + genofeature_columns + ['module']
contig_df = contig_df[cols]

contig_df.to_csv("contig_df.tsv", sep='\t', index=False)

# stats file
module_counts = contig_df['module'].value_counts().reset_index()
module_counts.columns = ['module', 'count']
module_counts.to_csv("module_stats.tsv", sep='\t', index=False)
    """
}


process MODIFY_CLUSTERS {
    publishDir "${params.outdir}/clusters",
        mode: 'copy'
        
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    input:
        path cluster_csv_dir // Directory containing the cluster CSV files
        
    output:
        path "hdbscan_modified_cluster_labels.tsv", emit: modified_clusters

    script:
    """
#!/usr/bin/env python3
import os
import pandas as pd

# Original Version from Barb Ferrell 05/2025
# Alicia Holk aholk@udel.edu

#nn 100 and clusters 30 and 10
print()

# Hardcoded paths and parameters
file1 = os.path.join("${cluster_csv_dir}", "hdbscan_nn100_md7_minclust30_clusters.csv")
file2 = os.path.join("${cluster_csv_dir}", "hdbscan_nn100_md7_minclust10_clusters.csv")

# Create a dataframe from each file
df1 = pd.read_csv(file1, sep=",", header=0)
df2 = pd.read_csv(file2, sep=",", header=0)


# Get the maximum cluster_label in minclust30 (to avoid overlap), assign next_cluster as max + 1
max_cluster_label = df1['cluster_label'].max()
next_cluster = max_cluster_label + 1

# Add a new column 'label_source' in df1 to track the origin of the cluster label
df1['label_source'] = 'minclust30'  # Default all to minclust30

# Get all rows in df1 where cluster_label == -1 (unassigned)
df1_unclustered = df1[df1['cluster_label'] == -1]

# Get the corresponding rows in df2 with the same embedding_ids
embedding_ids_unclustered = df1_unclustered['embedding_id']
matching_rows = df2[df2['embedding_id'].isin(embedding_ids_unclustered)]

# Get the unmatched rows from df2 (rows that don't have embedding_ids in df1_unclustered)
unmatched_rows = df2[~df2['embedding_id'].isin(embedding_ids_unclustered)]

# Sort matching_rows by cluster_label
matching_rows = matching_rows.sort_values(by='cluster_label')

# For each unique cluster_label in matching_rows
for cluster_label in matching_rows['cluster_label'].unique():
    if cluster_label > -1:  # Only process non-negative cluster labels
        # Check if cluster_label exists in unmatched_rows
        if cluster_label not in unmatched_rows['cluster_label'].values:
            # If it doesn't exist, assign the next available cluster_label in df1
            cluster_mask = matching_rows['cluster_label'] == cluster_label
            embedding_ids = matching_rows[cluster_mask]['embedding_id']

            # Update the cluster_label for all matching embedding_ids in df1
            for emb_id in embedding_ids:
                df1.loc[df1['embedding_id'] == emb_id, 'cluster_label'] = next_cluster
                df1.loc[df1['embedding_id'] == emb_id, 'label_source'] = f"minclust10_{cluster_label}"

            # Increment the cluster counter for next assignment
            next_cluster += 1

# Ensure all remaining 'label_source' values for unprocessed rows are "minclust30"
df1['label_source'].fillna('minclust30', inplace=True)

# Sort df1 by 'cluster_label' before saving
df1_sorted = df1.sort_values(by='cluster_label')

# Save the modified dataframe to a new file
df1_sorted.to_csv("hdbscan_modified_cluster_labels.tsv", index=False, sep="\t")
    """
}

process ZEROFIVEC {
    publishDir "${params.outdir}/clusters_imgs",
        mode: 'copy'
        
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    input:
        path coordinates_dir  // Directory containing the pre-generated coordinates
        path modified_clusters  // Modified clusters from previous process
        
    output:
        path "umap_nn*"

    script:
    """
#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import umap
import pandas as pd
import matplotlib.lines as mlines
import random  # <-- Add this import

#----------v_25_06_25_umap_coordinates---------------------------------
# refers to plotted coordinates from 04a umap script to cluster and map providing consistency and efficiency
# updated to match 04a umap section for consistency. Map is plotting differently than 04a and 05a for some reason
#--------- v_25_06_17------------------------------------
# added script from Zach to plot more consistently
# sections updated include SET SEED and Perform UMAP 
# Alicia Holk aholk@udel.edu
# ----------------------------------------------------



# Define parameter ranges
umap_nn = 100
umap_md = 0.7


def load_cluster_metadata(metadata_file):

    try:
        metadata_df = pd.read_csv(metadata_file, sep="\t")
        if "embedding_id" not in metadata_df.columns or "cluster_label" not in metadata_df.columns:
            raise ValueError("Metadata file must contain 'embedding_id' and 'cluster_label' columns.")
        # Remove exact duplicate rows
        metadata_df = metadata_df.drop_duplicates()

        # Identify duplicated embedding_ids
        duplicated_ids = metadata_df["embedding_id"].duplicated(keep=False)

        # Remove only the rows where embedding_id is duplicated AND cluster_label is -1
        metadata_df = metadata_df[~(duplicated_ids & (metadata_df["cluster_label"] == -1))]

        # Reset index after cleanup
        metadata_df = metadata_df.reset_index(drop=True)
        return metadata_df
    except Exception as e:
        print(f"Error loading metadata file: {e}")
        sys.exit(1)


def plot_umap_with_metadata(umap_nn, umap_md, cluster_df, output_file):
    
    # Load UMAP coordinates from 04a (instead of recomputing)
    coords_file = f"${coordinates_dir}/coordinates_nn{umap_nn}_md{int(umap_md * 10)}.tsv"

    # Load coordinates and IDs from 04a
    coord_df = pd.read_csv(coords_file, sep='\t')
    umap_ids = coord_df['embedding_id'].tolist()
    umap_coords = coord_df[['x', 'y']].values
    
    # Create mapping from embedding_id to coordinates
    id_to_coords = {eid: coord for eid, coord in zip(umap_ids, umap_coords)}
    
    if len(common_ids) == 0:
        print("Error: No common embedding IDs found between 04a coordinates and current embeddings")
        return
    
    # Get coordinates for common IDs in the same order as they appear in common_ids
    mapper = np.array([id_to_coords[eid] for eid in common_ids])
    
    print(f"Using {len(common_ids)} embeddings with coordinates from 04a")
    
    # Convert common embedding IDs to a dataframe for merging
    embedding_df = pd.DataFrame({"embedding_id": common_ids})
    # Merge with metadata to assign cluster labels
    merged_df = embedding_df.merge(cluster_df, on="embedding_id", how="left")
    labels = merged_df["cluster_label"].fillna(-1).values  # Assign -1 to unclustered points

    clustered = labels >= 0  # Identify clustered points

    # Plot UMAP without colorbar
    fig, ax = plt.subplots(figsize=(10, 7))
    scatter = ax.scatter(
        mapper[clustered, 0],
        mapper[clustered, 1],
        c=labels[clustered],
        s=10, alpha=0.5, cmap="Spectral"
    )

    ax.scatter(
        mapper[~clustered, 0],
        mapper[~clustered, 1],
        color="gray", s=10, alpha=0.5, label="not clustered"
    )

    # Annotate cluster centers
    unique_labels = np.unique(labels[labels >= 0])
    for cluster in unique_labels:
        cluster_points = mapper[labels == cluster]
        center_x, center_y = cluster_points[:, 0].mean(), cluster_points[:, 1].mean()
        ax.text(center_x, center_y, str(cluster), fontsize=10, ha='center', va='center')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # plt.title(f"UMAP with merged HDBSCAN Clustering (nn={umap_nn}, umap_md={umap_md})")
    plt.grid(False)
    plt.tight_layout()

    # Save main plot (without colorbar)
    plt.savefig(f"{output_file}.png", dpi=600)
    plt.savefig(f"{output_file}.svg", transparent=True)
    plt.close()
    print(f"Saved: {output_file}")

    # Save colorbar as a separate figure
    colorbar_fig, colorbar_ax = plt.subplots(figsize=(1.2, 4))  # Adjust size as needed

    norm = plt.Normalize(vmin=labels[clustered].min(), vmax=labels[clustered].max())
    cbar = plt.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap="Spectral"),
        cax=colorbar_ax,
        orientation='vertical'
    )

    cbar.set_label("Cluster Labels")
    cbar.ax.hlines(90, 0, 1, color='black', linestyle='--', linewidth=2)

    colorbar_fig.savefig(f"{output_file}_colorbar.png", dpi=300, bbox_inches="tight")
    colorbar_fig.savefig(f"{output_file}_colorbar.svg", dpi=300, bbox_inches="tight")
    plt.close(colorbar_fig)


if __name__ == "__main__":
    cluster_metadata = "${modified_clusters}"  # Metadata file with colors/markers

    # Load cluster labels
    cluster_df = load_cluster_metadata(cluster_metadata)
    cluster_df.to_csv("hdbscan_modified_cluster_labels_no_duplicates.tsv")

    embedding_ids = cluster_df["embedding_id"].tolist()

    if len(embedding_ids) > 0:
        output_file = f"umap_nn{umap_nn}_umapmd{int(umap_md * 10)}_with_combined_hdbscan_clusters"
        plot_umap_with_metadata(umap_nn, umap_md, cluster_df, output_file)
    else:
        print("No valid embeddings found.")
    """
}



process GENOFEATURE_CENTRIC {
    publishDir "${params.outdir}/",
        mode: 'copy'
        
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    input:
        path module_file  // Path to the module file with genofeature information
        path coordinates_file
        path connections_file  // Path to the connections file

        
    output:
        path "genofeature_plots/"


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
metadata_file = "${params.genofeature_metadata}"  # Path to metadata file
connections_file = "${connections_file}"

coord_df = pd.read_csv(coords_file, sep='\t')
embedding_id_to_index = {eid: idx for idx, eid in enumerate(coord_df['embedding_id'])}
embedding_2d = coord_df[['x', 'y']].values

module_df = pd.read_csv(module_file, sep='\t')
metadata_df = pd.read_csv(metadata_file, sep="\t")
color_map = dict(zip(metadata_df["genofeature"], metadata_df["color"]))
display_map = dict(zip(metadata_df["genofeature"], metadata_df["display_name"]))
marker_map = dict(zip(metadata_df["genofeature"], metadata_df["marker"]))
genofeature_cols = [col for col in module_df.columns if col not in ['dataset', 'contig_id', 'module']]

connections = []
if os.path.exists(connections_file):
    connections_df = pd.read_csv(connections_file, sep='\t')
    connections = list(zip(connections_df['id1'], connections_df['id2']))

visible_features_list = ['rPolB', 'pPolB', 'piPolB']

plot_output_dir = "genofeature_plots"
os.makedirs(plot_output_dir, exist_ok=True)

visible_features_list = "${params.selected_genofeatures.join(',')}".split(',')

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

    plt.savefig(os.path.join(plot_output_dir, f"umap_{feature}_recolor_plot.png"), dpi=600, bbox_inches="tight")
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
        fig, ax = plt.subplots(figsize=(10, 7))

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
        plt.savefig(os.path.join(plot_output_dir, f"umap_{feature}_by_dataset_{dataset}.png"), dpi=600, bbox_inches="tight")
        plt.close()

# --------- AUTOMATIC TILING ---------
def tile_images(image_dir, output_file, padding=10):
    png_files = [os.path.join(image_dir, f) for f in os.listdir(image_dir) 
                 if f.lower().endswith("_plot.png")]
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
tiled_output_file = os.path.join(plot_output_dir, "genofeature_highlighting_tiled_all.png")
tile_images(plot_output_dir, tiled_output_file)
print(f"All plots completed!")
    """
}


workflow {
    ch_cluster_dir = Channel.fromPath(params.clusters_dir)
    ch_filtered_tsv = Channel.fromPath(params.second_run)
    ch_metadata = Channel.fromPath(params.genofeature_metadata)

    ch_module_file = MODULE_FILE(
        ch_filtered_tsv,
        ch_metadata
    )

    GENERATE_COORDINATES_2(params.embeddings)
    ch_coordinates = GENERATE_COORDINATES_2.out.coordinates_tsv
    ch_connections = GENERATE_COORDINATES_2.out.connections_tsv
    // MODIFY_CLUSTERS(ch_cluster_dir)
    // ZEROFIVEC(GENERATE_COORDINATES_2.coordinates_tsv, MODIFY_CLUSTERS.out)

    GENOFEATURE_CENTRIC(ch_module_file, ch_coordinates, ch_connections)

}