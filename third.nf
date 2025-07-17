nextflow.enable.dsl=2

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
import os
import pandas as pd

# Original Version from Barb Ferrell 05/2025
# Alicia Holk aholk@udel.edu

#nn 100 and clusters 30 and 10

# Hardcoded paths and parameters
input_dir = "${cluster_csv_dir}"
file1 = os.path.join(input_dir, "hdbscan_nn100_umapmd7_md0_minclust30_clusters.csv")
file2 = os.path.join(input_dir, "hdbscan_nn100_umapmd7_md0_minclust10_clusters.csv")


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

process 05c{
    publishDir "${params.outdir}/clusters",
        mode: 'copy'
        
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    input:
        path modified_clusters  // Modified clusters from previous process
        
    output:
        path "plots/*.{png,svg}"

    script:
    """
import os
import sys
import torch
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

# SET SEED for packages to default all randomization (Include at start of script)
os.environ['PYTHONHASHSEED'] = '42'
random.seed(42)
np.random.seed(42)
torch.manual_seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed(42)
    torch.cuda.manual_seed_all(42)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

# Define parameter ranges
umap_nn = 100
umap_md = 0.7

# Hardcoded paths and parameters
input_dir = "${embeddings_dir}"  # Updated to filtered embeddings directory
coords_dir = "${coordinates_dir}"  # Directory with saved coordinates from 04a
metadata_file = "${metadata_file}"  # Metadata file with colors/markers
output_dir = "${output_dir}"  # Output directory for UMAP plots
os.makedirs(output_dir, exist_ok=True)


def load_embedding(file_path):
    try:
        embedding_data = torch.load(file_path)
        embedding = embedding_data["mean_representations"][36].numpy()
        embedding_id = embedding_data.get("label")
        if embedding_id is None:
            raise ValueError(f"Embedding ID not found in file: {os.path.basename(file_path)}")
        return embedding, embedding_id
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return None, None


def load_all_embeddings(base_dir, subdirs=None):

    embeddings = []
    embedding_ids = []

    # If subdirs is None, scan all subdirectories
    if subdirs is None:
        subdirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]

    for subdir in subdirs:
        subdir_path = os.path.join(base_dir, subdir)
        if not os.path.isdir(subdir_path):
            print(f"Skipping {subdir_path}, not a valid directory")
            continue

        # Recursively search for .pt files
        for root, _, files in os.walk(subdir_path):  # Walk through all levels
            for file_name in files:
                if file_name.endswith(".pt"):
                    file_path = os.path.join(root, file_name)
                    embedding, embedding_id = load_embedding(file_path)
                    if embedding is not None and embedding_id is not None:
                        embeddings.append(embedding)
                        embedding_ids.append(embedding_id)

    return embeddings, embedding_ids


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


def plot_umap_with_metadata(embedding_ids, umap_nn, umap_md, cluster_df, output_file, coords_dir):
    
    # Load UMAP coordinates from 04a (instead of recomputing)
    coords_file = os.path.join(coords_dir, f"embedding_2d_nn{umap_nn}_md{int(umap_md * 10)}.npy")
    ids_file = os.path.join(coords_dir, f"embedding_ids_nn{umap_nn}_md{int(umap_md * 10)}.txt")
    
    if not os.path.exists(coords_file) or not os.path.exists(ids_file):
        print(f"Required coordinate files not found. Please run 04a script first.")
        print(f"Looking for: {coords_file}")
        print(f"Looking for: {ids_file}")
        return
    
    # Load coordinates and IDs from 04a
    umap_coords = np.load(coords_file)
    with open(ids_file, 'r') as f:
        umap_ids = [line.strip() for line in f]
    
    print(f"Loaded {len(umap_coords)} UMAP coordinates from 04a")
    print("First 10 embedding_ids from 04a:", umap_ids[:10])
    
    # Create mapping from embedding_id to coordinates
    id_to_coords = {eid: coord for eid, coord in zip(umap_ids, umap_coords)}
    
    # Filter to only include IDs that are in both the coordinate data and our embedding_ids
    common_ids = set(embedding_ids) & set(umap_ids)
    common_ids = list(common_ids)
    
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
    # Load cluster labels
    cluster_df = load_cluster_metadata(metadata_file)
    cluster_df.to_csv(os.path.join(output_dir, "hdbscan_modified_cluster_labels_no_duplicates.tsv"), sep="\t", index=False, header=True)
    
    # Load embedding IDs only (we don't need full embeddings since we're using coordinates from 04a)
    selected_subdirs = ["POL", "RNR", "PolB", "PolC", "RNR_classIbeta", "helicase", "RNR_classIII"]
    embeddings, embedding_ids = load_all_embeddings(input_dir, selected_subdirs)

    if len(embedding_ids) > 0:
        output_file = os.path.join(output_dir, f"umap_nn{umap_nn}_umapmd{int(umap_md * 10)}_with_combined_hdbscan_clusters")
        plot_umap_with_metadata(embedding_ids, umap_nn, umap_md, cluster_df, output_file, coords_dir)
    else:
        print("No valid embeddings found.")
    """
}

process 05d {
    publishDir "${params.outdir}/umap_cluster",
        mode: 'copy'
        
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    input:
        path coordinates_dir  // Directory containing the pre-generated coordinates
        path filtered_tsv    // TSV file with metadata
        path metadata_file   // Metadata file with colors/markers
        
    output:
        path "plots/*.{png,svg}"

    script:
    """
import os
import torch    
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from collections import defaultdict
import random  

#---------v_25_06_25------------------------------------
# UPDATED: This script now loads UMAP coordinates from 04a script instead of recalculating them
# This ensures consistency and efficiency across all downstream scripts
# UMAP coordinates are loaded from files saved by 04a_umap_filter_25_06_25.py

#--------- v_25_06_19-------------
# saves umap coordinates and order of embedding ids
# puts edges behind nodes in network plot (z-order=1)
# added script from Zach to plot more consistently
# sections updated include SET SEED and Perform UMAP 

# output directory name is now "05d_umap_network" instead of "05d_umap_network_v4run"
# Alicia Holk aholk@udel.edu
# ----------------------------------------------------

# SET SEED for packages to default all randomization (Include at start of script)
os.environ['PYTHONHASHSEED'] = '42'
random.seed(42)
np.random.seed(42)
torch.manual_seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed(42)
    torch.cuda.manual_seed_all(42)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

# Define parameter ranges
n_neighbors_list = [100]
min_dist_list = [0.7]

# Hardcoded paths and parameters
input_dir = "/mnt/VEIL/users/aholk/ENA_2025/03c_embeddings/"  
module_file = "/mnt/VEIL/users/aholk/ENA_2025/02e_orf_ctg_summary/orf_df.tsv"
metadata_file = "/mnt/VEIL/users/aholk/ENA_2025/bin/protein_metadata_all.txt"
column_name = "genofeature"
coords_dir = "/mnt/VEIL/users/aholk/ENA_2025/04a_umap"  # Directory containing UMAP coordinates from 04a

output_dir = "/mnt/VEIL/users/aholk/ENA_2025/05d_umap_network"
os.makedirs(output_dir, exist_ok=True)

def parse_embedding_id(embedding_id): #taking orf ID and returning the contig ID
    parts = embedding_id.split('_')
    if len(parts) >= 3:
        return '_'.join(parts[:-3])
    return embedding_id

def load_embedding(file_path):
    """Loads embeddings from a .pt file."""
    try:
        embedding_data = torch.load(file_path)
        embedding = embedding_data["mean_representations"][36].numpy()
        embedding_id = embedding_data.get("label")
        if embedding_id is None:
            raise ValueError(f"Embedding ID not found in file: {os.path.basename(file_path)}")
        return embedding, embedding_id
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return None, None

def load_module(module_file):
    """Loads module data from a TSV file."""
    try:
        module_df = pd.read_csv(module_file, sep='\t')
        module_df["normalized_orf_id"] = module_df["orf_id"].str.replace(".", "-", regex=False)
        return module_df
    except Exception as e:
        print(f"Error loading module file: {e}")
        return None


def load_metadata(metadata_file):
    """Loads metadata (including colors and shapes) from a file into dictionaries."""
    try:
        metadata_df = pd.read_csv(metadata_file, sep="\t")  # Ensure tab-delimited format
        # print(metadata_df.head())
        color_map = dict(zip(metadata_df["genofeature"], metadata_df["color"]))  # Feature → Color
        display_map = dict(zip(metadata_df["genofeature"], metadata_df["display_name"]))  # Feature → Display Name
        marker_map = dict(zip(metadata_df["genofeature"], metadata_df["marker"]))  # Feature → Marker Shape
        return color_map, display_map, marker_map
    except Exception as e:
        print(f"Error loading metadata file: {e}")
        return {}, {}, {}

def load_all_embedding_ids(base_dir, subdirs=None):
    """
    Recursively loads embedding IDs from .pt files in the specified subdirectories.
    Only loads IDs, not the actual embeddings.

    Parameters:
        base_dir (str): The root directory containing subdirectories with .pt files.
        subdirs (list, optional): A list of subdirectories to process. If None, all subdirectories are processed.

    Returns:
        list: A list of embedding IDs.
    """
    embedding_ids = []

    # If subdirs is None, scan all subdirectories
    if subdirs is None:
        subdirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]

    for subdir in subdirs:
        subdir_path = os.path.join(base_dir, subdir)
        if not os.path.isdir(subdir_path):
            print(f"Skipping {subdir_path}, not a valid directory")
            continue

        # Recursively search for .pt files
        for root, _, files in os.walk(subdir_path):  # Walk through all levels
            for file_name in files:
                if file_name.endswith(".pt"):
                    file_path = os.path.join(root, file_name)
                    embedding, embedding_id = load_embedding(file_path)
                    if embedding_id is not None:
                        embedding_ids.append(embedding_id)

    return embedding_ids


def plot_umap(embedding_ids, module_df, nn, md, output_file, display_names):
    """Plots UMAP embedding with colors mapped from metadata using coordinates from 04a and saves legend as separate image."""
    
    # Load UMAP coordinates from 04a
    md_int = int(md * 10)  # Convert min_dist to integer for filename
    coord_file = os.path.join(coords_dir, f"embedding_2d_nn{nn}_md{md_int}.npy")
    coord_ids_file = os.path.join(coords_dir, f"embedding_ids_nn{nn}_md{md_int}.txt")
    
    try:
        # Load coordinates and IDs from 04a
        embedding_2d = np.load(coord_file)
        with open(coord_ids_file, 'r') as f:
            coord_embedding_ids = [line.strip() for line in f]
        
        print(f"Loaded UMAP coordinates from {coord_file}")
        print(f"Coordinate shape: {embedding_2d.shape}")
        print(f"Number of coordinate IDs: {len(coord_embedding_ids)}")
        
    except FileNotFoundError as e:
        print(f"Error: Could not find coordinate files from 04a: {e}")
        print(f"Expected files: {coord_file} and {coord_ids_file}")
        return
    
    # Align the embedding IDs with coordinate IDs
    coord_id_to_index = {eid: i for i, eid in enumerate(coord_embedding_ids)}
    
    # Filter to only embedding_ids that have coordinates
    valid_indices = []
    valid_embedding_ids = []
    
    for eid in embedding_ids:
        if eid in coord_id_to_index:
            valid_indices.append(coord_id_to_index[eid])
            valid_embedding_ids.append(eid)
    
    if not valid_indices:
        print("Warning: No embedding IDs match coordinate IDs")
        return
    
    # Get coordinates for valid embeddings
    valid_embedding_2d = embedding_2d[valid_indices]
    
    print(f"Aligned {len(valid_indices)} embeddings with coordinates")
    print("First 10 valid embedding_ids:", valid_embedding_ids[:10])
   
    # Normalize orf_name for consistent lookup
    module_df["normalized_orf_id"] = module_df["orf_id"].str.replace(".", "-", regex=False)

    # Create lookup dictionaries for colors and markers
    id_to_color = dict(zip(module_df["normalized_orf_id"], module_df["manual_color"]))
    id_to_marker = dict(zip(module_df["normalized_orf_id"], module_df["manual_marker"]))

    # Assign colors and markers to embedding IDs
    colors = [id_to_color.get(embedding_id, "#808080") for embedding_id in valid_embedding_ids]  # Default: gray
    # markers = [id_to_marker.get(embedding_id, ".") for embedding_id in embedding_ids]      # Default: '.'
    markers = [id_to_marker.get(embedding_id, ".") if pd.notna(id_to_marker.get(embedding_id)) else '.' for embedding_id in valid_embedding_ids] #if NaN or not in list, put a "." as marker style

    # Group points for connecting lines - finding things that occur on the same contig
    groups = defaultdict(list)
    for i, embedding_id in enumerate(valid_embedding_ids):
        parsed_id = parse_embedding_id(embedding_id)
        groups[parsed_id].append(i)

    # Create UMAP figure (without legend)
    fig, ax = plt.subplots(figsize=(10, 7))

    # Draw connections
    for parsed_id, indices in groups.items():
        if len(indices) > 1:
            for i in range(len(indices) - 1):
                for j in range(i + 1, len(indices)):
                    idx1, idx2 = indices[i], indices[j]
                    x1, y1 = valid_embedding_2d[idx1, 0], valid_embedding_2d[idx1, 1]
                    x2, y2 = valid_embedding_2d[idx2, 0], valid_embedding_2d[idx2, 1]
                    ax.plot([x1, x2], [y1, y2], color='#CCCCCC', alpha=0.5, linewidth=0.2,zorder=1)#only change in v3 from original script

    # Scatter plot points
    for embedding, color, marker in zip(valid_embedding_2d, colors, markers):
        ax.scatter(embedding[0], embedding[1], c=color, marker=marker, s=10, alpha=0.7)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.grid(False)
    plt.tight_layout()

    # Save UMAP figure (without legend)
    plt.savefig(f"{output_file}.png", dpi=600, bbox_inches="tight")
    plt.savefig(f"{output_file}.svg", transparent=True)
    print(f"Saved: {output_file}")

    # ----- Separate Legend Creation -----

    # Drop duplicates for one legend entry per genofeature
    unique_features = module_df.drop_duplicates(subset=["genofeature"], keep="first")

    # Filter to only those with known display names
    existing_keys = unique_features["genofeature"].isin(display_names.keys())
    unique_features = unique_features[existing_keys]

    # Add display names
    unique_features["display_name"] = unique_features["genofeature"].map(display_names)

    # Order and sort
    unique_features["genofeature"] = pd.Categorical(
        unique_features["genofeature"],
        categories=list(display_names.keys()),
        ordered=True
    )
    unique_features = unique_features.sort_values("genofeature")

    # Create legend handles
    legend_handles = [
        mlines.Line2D([], [], color=row["manual_color"], marker=row["manual_marker"], linestyle="None", markersize=8,
                      label=row["display_name"])
        for _, row in unique_features.iterrows()
    ]

    # Create a new figure just for the legend
    legend_fig = plt.figure(figsize=(4, len(legend_handles) * 0.3))
    legend_ax = legend_fig.add_subplot(111)
    legend_ax.axis("off")
    legend_ax.legend(handles=legend_handles, title="Genofeature", loc="center left", frameon=False)

    legend_fig.savefig(f"{output_file}_legend.png", dpi=300, bbox_inches="tight")
    legend_fig.savefig(f"{output_file}_legend.svg", dpi=300, bbox_inches="tight")
    print(f"Saved legend: {output_file}_legend")

    # Close figures
    plt.close(fig)
    plt.close(legend_fig)


if __name__ == "__main__":
    # Load data
    # Load embedding IDs only (we don't need full embeddings since we're using coordinates from 04a)
    selected_subdirs = ["POL", "RNR", "PolB", "PolC", "RNR_classIbeta", "helicase", "RNR_classIII"]

    embedding_ids = load_all_embedding_ids(input_dir, selected_subdirs)

    module_df = load_module(module_file)
    color_map, display_names, marker_map = load_metadata(metadata_file)
    # print(color_map)
    # print(display_names)
    # print(marker_map)

    if module_df is not None:
        # Map colors to metadata
        module_df['manual_color'] = module_df[column_name].map(color_map)
        module_df['display_names'] = module_df[column_name].map(display_names)
        module_df['manual_marker'] = module_df[column_name].map(marker_map)
        print("Assigned Colors:", color_map)
        print("Assigned Labels:", display_names)
        print("Assigned Shapes:", marker_map)

    print(module_df.head())

    if len(embedding_ids) > 0:
        for nn in n_neighbors_list:
            for md in min_dist_list:
                md_int = int(md * 10)  # Convert min_dist to integer for filename
                output_file = os.path.join(output_dir, f"umap_nn{nn}_md{md_int}")
                plot_umap(embedding_ids, module_df, nn, md, output_file, display_names)
    else:
        print("No valid embeddings found.")

    """
}

process 05e {
    publishDir "${params.outdir}/umap_projection",
        mode: 'copy'
        
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    input:

        
    output:


    script:
    """
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from collections import defaultdict
from PIL import Image
import math

#----------v_25_06_26----------
# Added automatic tiling feature at the end to create a combined image of all plots
#----------v_25_06_25----------
# UPDATED: Now loads UMAP coordinates directly from 04a directory for consistency
# Uses the same coordinate files as all other downstream scripts
# Pointed to updated ctg_df_filtered.tsv input file name
#----------v_25_06_20-----------
# IMPORTANT UPDATE: 
# previous scripts when coloring GP4 would also pull GP41 - This script has been modified to find "exact" matches rather than containing matches.
#----------------------------------------------------
# This script recolors UMAP embeddings based on specified 
# genofeatures using the consistent UMAP coordinates from 04a.
# This ensures all scripts use the exact same coordinate system.

# UPDATED v25_06_25: Now loads coordinates directly from 04a_umap directory
# instead of from 05d to ensure consistency across all downstream scripts.
# Alicia Holk aholk@udel.edu
#----------------------------------------------------

# Update paths as needed
base_dir = "/mnt/VEIL/users/aholk/ENA_2025"
coords_dir = os.path.join(base_dir, "04a_umap")  # Use coordinates directly from 04a for consistency
# UMAP parameters to match other scripts
umap_nn = 100
umap_md = 0.7
md_int = int(umap_md * 10)  # Convert to integer for filename
embedding_2d_file = os.path.join(coords_dir, f"embedding_2d_nn{umap_nn}_md{md_int}.npy")
embedding_ids_file = os.path.join(coords_dir, f"embedding_ids_nn{umap_nn}_md{md_int}.txt")
output_dir = os.path.join(base_dir, "05e_genofeature_highlighting")
os.makedirs(output_dir, exist_ok=True)

module_file = os.path.join(base_dir, "02e_orf_ctg_summary", "ctg_df.tsv")
metadata_file = os.path.join(base_dir, "bin", "protein_metadata_all.txt")

selected_subdirs = ["POL", "RNR", "PolB", "PolC", "RNR_classIbeta", "helicase", "RNR_classIII"]


# --------- PARAMETERS ---------
visible_features_list = [
    ['RDKF'],
    ['RDKY'],
    ['RDKH'],
    ['RDKL'],
    ['TDKY'],
    ['RDNF'],
    ['EDKY'],
    ['KDKL'],
    ['SDKL'],
    ['DDKL'],
    ['CDKY'],
    ['NCECI'],
    ['NCECV'],
    ['NCECP'],
    ['NCECA'],
    ['QCECL'],
    ['NCECL'],
    ['NCECM'],
    ['RNRclassIbeta'],
    ['RNRclassIII'],
    ['rPolB'],
    ['pPolB'],
    ['piPolB'],
    ['RecA'],
    ['UvsX'],
    ['RecB'],
    ['RecD'],
    ['PcrA'],
    ['UvrD'],
    ['Dda'],
    ['UvsW'],
    ['UvrB'],
    ['RecG'],
    ['SNF2'],
    ['Twinkle'],
    ['DnaB'],
    ['Gp4'],
    ['Gp41'],
    ['clamp'],
    ['primase']
 #['RDKF', 'RDKY'],  # Example: plot both together
    # Add more lists as needed
    # make sure it keeps all the points
]

# --------- LOAD DATA ---------
# Note: We no longer load embeddings since we use precomputed coordinates from 04a

# Load UMAP coordinates from 04a (always use precomputed coordinates for consistency)
if os.path.exists(embedding_2d_file) and os.path.exists(embedding_ids_file):
    print(f"Loading UMAP coordinates from 04a: {embedding_2d_file}")
    embedding_2d = np.load(embedding_2d_file)
    with open(embedding_ids_file) as f:
        embedding_ids = [line.strip() for line in f]
    print(f"Loaded {len(embedding_ids)} embedding IDs and coordinates from 04a")
else:
    print(f"Error: Required coordinate files not found in 04a directory:")
    print(f"  Expected: {embedding_2d_file}")
    print(f"  Expected: {embedding_ids_file}")
    print("Please run 04a_umap_filter_25_06_25.py first to generate coordinates.")
    exit(1)
embedding_id_to_index = {eid: idx for idx, eid in enumerate(embedding_ids)}

module_df = pd.read_csv(module_file, sep='\t')
genofeature_cols = [col for col in module_df.columns if 'genofeature' in col]
module_df["module"] = module_df[genofeature_cols].astype(str).apply(
    lambda row: "_".join(val for val in row if val != 'nan'), axis=1
)

metadata_df = pd.read_csv(metadata_file, sep="\t")
color_map = dict(zip(metadata_df["genofeature"], metadata_df["color"]))
display_map = dict(zip(metadata_df["genofeature"], metadata_df["display_name"]))
marker_map = dict(zip(metadata_df["genofeature"], metadata_df["marker"]))

def map_orf_to_metadata(module_df, color_map, display_map, marker_map):
    orf_metadata = {}
    for col in module_df.columns:
        if 'orf_id' in col:
            genofeature_col = col.replace('orf_id', 'genofeature')
            if genofeature_col in module_df.columns:
                for idx, orf_id in enumerate(module_df[col]):
                    genofeature = module_df[genofeature_col].iloc[idx]
                    if genofeature in color_map:
                        orf_metadata[orf_id] = {
                            'color': color_map[genofeature],
                            'name': display_map[genofeature],
                            'marker': marker_map[genofeature]
                        }
    return orf_metadata

orf_metadata = map_orf_to_metadata(module_df, color_map, display_map, marker_map)

def define_connections(module_df, visible_features):
    connections = []
    visible_orf_ids = set()
    legend_features = set()
    for _, row in module_df.iterrows():
        module_features = set(row['module'].split('_'))
        if any(feature in module_features for feature in visible_features):
            orf_ids = [row[col] for col in module_df.columns if 'orf_id' in col]
            visible_orf_ids.update(orf_ids)
            legend_features.update(module_features)
            for i in range(len(orf_ids)):
                for j in range(i + 1, len(orf_ids)):
                    connections.append((orf_ids[i], orf_ids[j]))
    return connections, visible_orf_ids, legend_features
# def define_connections(module_df, visible_features):
#     connections = []
#     visible_orf_ids = set()
#     legend_features = set()
#     for _, row in module_df.iterrows():
#         if any(feature in row['module'] for feature in visible_features):
#             orf_ids = [row[col] for col in module_df.columns if 'orf_id' in col]
#             visible_orf_ids.update(orf_ids)
#             module_features = row['module'].split('_')
#             legend_features.update(module_features)
#             for i in range(len(orf_ids)):
#                 for j in range(i + 1, len(orf_ids)):
#                     connections.append((orf_ids[i], orf_ids[j]))
#     return connections, visible_orf_ids, legend_features

# --------- PLOTTING ---------
for visible_features in visible_features_list:
    # Determine which ORFs/contigs to highlight
    connections, visible_orf_ids, legend_features = define_connections(module_df, visible_features)

    # Skip if none of the features are present
    if not visible_orf_ids:
        print(f"Skipping {visible_features}: not present in data.")
        continue

    fig, ax = plt.subplots(figsize=(10, 7))
    # 1. Plot grayed-out points (background)
    for eid in embedding_ids:

        # Skip if the embedding ID is not in the mapping
        if pd.isna(eid) or eid not in embedding_id_to_index:
            continue

        if eid not in visible_orf_ids:
            idx = embedding_id_to_index[eid]
            metadata = orf_metadata.get(eid, {'color': '#808080', 'marker': '.'})
            ax.scatter(
                embedding_2d[idx, 0], embedding_2d[idx, 1],
                c="#D3D3D3", marker=metadata['marker'], s=10, alpha=0.05, zorder=1
            )

    # 2. Plot highlighted edges (middle)
    for orf1, orf2 in connections:
        if orf1 in visible_orf_ids and orf2 in visible_orf_ids:
            idx1 = embedding_id_to_index.get(orf1)
            idx2 = embedding_id_to_index.get(orf2)
            if idx1 is not None and idx2 is not None:
                ax.plot(
                    [embedding_2d[idx1, 0], embedding_2d[idx2, 0]],
                    [embedding_2d[idx1, 1], embedding_2d[idx2, 1]],
                    color="gray", alpha=0.2, linewidth=0.5, zorder=2
                )

    # 3. Plot highlighted points (foreground, after edges)
    for eid in visible_orf_ids:

        # Skip if the embedding ID is not in the mapping
        if pd.isna(eid) or eid not in embedding_id_to_index:
           continue
        
        idx = embedding_id_to_index[eid]
        metadata = orf_metadata.get(eid, {'color': '#808080', 'marker': '.'})
        ax.scatter(
            embedding_2d[idx, 0], embedding_2d[idx, 1],
            c=metadata['color'], marker=metadata['marker'], s=10, alpha=0.7, zorder=3
        )

    ax.set_xticks([])
    ax.set_yticks([])

    # Add visible_features label below the plot
    if visible_features:
        label = ", ".join(visible_features)
        plt.figtext(0.5, -0.05, label, ha='center', va='top', fontsize=14)

    # Save plot
    visible_str = "_".join(visible_features)
    output_filename = f"umap_{visible_str}_recolor"
    output_file_base = os.path.join(output_dir, output_filename)
    plt.savefig(f"{output_file_base}_plot.png", dpi=600, bbox_inches="tight")
    plt.savefig(f"{output_file_base}_plot.svg", transparent=True, bbox_inches="tight")
    plt.close()

    # Legend
    master_genofeature_order = metadata_df["genofeature"].tolist()
    feature_handles = []
    feature_set = set(legend_features)
    for feature in master_genofeature_order:
        if feature in feature_set:
            feature_row = metadata_df[metadata_df["genofeature"] == feature]
            if not feature_row.empty:
                color = feature_row.iloc[0]["color"] if "color" in feature_row.columns else "#808080"
                marker = feature_row.iloc[0]["marker"] if "marker" in feature_row.columns else "."
                display_name = feature_row.iloc[0]["display_name"] if "display_name" in feature_row.columns else feature
                handle = mlines.Line2D([], [], color=color,
                                       marker=marker, linestyle="None",
                                       markersize=8, label=display_name)
                feature_handles.append(handle)
    fig, ax = plt.subplots(figsize=(3, max(1, len(feature_handles) * 0.3)))
    legend = ax.legend(handles=feature_handles,
                       loc="center left",
                       frameon=False,
                       handlelength=0.5,
                       handletextpad=1,
                       fontsize=6,
                       title="Feature Types",
                       title_fontsize=7,
                       ncol=1)
    ax.axis('off')
    plt.savefig(f"{output_file_base}_legend.png", dpi=600, bbox_inches="tight")
    plt.close()

# --------- AUTOMATIC TILING ---------
# Create a tiled image of all generated plots

def tile_images(image_dir, output_file, padding=10):
    """Tile all PNG plot files in the directory into a single image."""
    # Find all PNG plot files (excluding legend files)
    png_files = [os.path.join(image_dir, f) for f in os.listdir(image_dir) 
                 if f.lower().endswith("_plot.png")]
    png_files.sort()  # Sort alphabetically
    
    if not png_files:
        print("No plot PNG files found for tiling.")
        return
    
    print(f"Tiling {len(png_files)} plot images...")
    
    # Resize all images to the same height, keeping aspect ratio
    images = [Image.open(f) for f in png_files]
    heights = [img.height for img in images]
    target_height = max(heights)  # Use the height of the tallest image
    
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
    
    # Compute max width for each column
    col_widths = [0] * cols
    row_heights = [target_height] * rows  # all rows have the same height now
    for idx, (w, h) in enumerate(resized_sizes):
        row = idx // cols
        col = idx % cols
        if w > col_widths[col]:
            col_widths[col] = w
    
    # Compute x offsets for columns and y offsets for rows
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
    
    # Save the tiled image
    canvas.save(output_file)
    print(f"Tiled image saved as {output_file}")

# Create tiled image of all plots
tiled_output_file = os.path.join(output_dir, "genofeature_highlighting_tiled_all.png")
tile_images(output_dir, tiled_output_file)

print(f"\nAll plots completed! Tiled image available at: {tiled_output_file}")
    """
}


process UMAP_PROJECTION {
    publishDir "${params.outdir}/umap",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith("tiled_image.png")) filename
            else null
        }
        
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    input:
        path coordinates_dir  // Directory containing the pre-generated coordinates
        path filtered_tsv    // TSV file with metadata
        path metadata_file   // Metadata file with colors/markers
        
    output:
        path "plots/tiled_image.png", emit: tiled_image
        path "plots/*.png", emit: plots
        
    script:
    def selected_list = params.selected_genofeatures.collect { "\'${it}\'" }.join(', ')
    """
#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set non-interactive backend
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from collections import defaultdict
import platform
from PIL import Image, ImageDraw, ImageFont

# Create output directory
os.makedirs("plots", exist_ok=True)

def parse_embedding_id(embedding_id):
    parts = embedding_id.split('_')
    if len(parts) >= 3:
        return '_'.join(parts[:-3])
    return embedding_id

def plot_umap(module_df, nn, md, output_file, display_names):
    # Load pre-generated UMAP coordinates and IDs
    md_int = int(md * 10)
    coord_file = os.path.join("${coordinates_dir}", f"coords/coordinates_nn{nn}_md{md_int}.npy")
    id_file = os.path.join("${coordinates_dir}", f"ids/embedding_ids_nn{nn}_md{md_int}.txt")
    conns_file = os.path.join("${coordinates_dir}", f"connections/connections_nn{nn}_md{md_int}.npy")

    try:
        # Load coordinates and their corresponding IDs
        embedding_2d = np.load(coord_file)
        with open(id_file, 'r') as f:
            embedding_ids = [line.strip() for line in f]
        
        print(f"Loaded coordinates from {coord_file}")
        print(f"Shape: {embedding_2d.shape}, IDs: {len(embedding_ids)}")
        
        # Prepare metadata mappings
        module_df["normalized_orf_id"] = module_df["orf_id"].str.replace(".", "-", regex=False)
        id_to_genofeature = dict(zip(module_df["normalized_orf_id"], module_df["genofeature"]))
        genofeature_to_color = dict(zip(metadata_df["genofeature"], metadata_df["color"]))
        genofeature_to_marker = dict(zip(metadata_df["genofeature"], metadata_df["marker"]))
        genofeature_to_display = dict(zip(metadata_df["genofeature"], metadata_df["display_name"]))

        # Assign colors and markers
        colors = []
        markers = []
        for eid in embedding_ids:
            genofeature = id_to_genofeature.get(eid, None)
            colors.append(genofeature_to_color.get(genofeature, "#808080"))
            marker = genofeature_to_marker.get(genofeature, ".")
            markers.append(marker if pd.notna(marker) else ".")

        # Group by contig ID for connections
        groups = defaultdict(list)
        for i, eid in enumerate(embedding_ids):
            groups[parse_embedding_id(eid)].append(i)

        if os.path.exists(conns_file):
            connections = np.load(conns_file)
            print(f"Loaded {len(connections)} connections")
        else:
            connections = []
            print("No connections file found")

        # Create plot
        fig, ax = plt.subplots(figsize=(10, 7))

        # Draw connections first (background)
        if len(connections) > 0:
            # Plot all connections at once
            coords = embedding_2d[connections]
            ax.plot(coords[:, :, 0].T, coords[:, :, 1].T,
                   color='#CCCCCC', alpha=0.5, 
                   linewidth=0.2, zorder=1)

        # Plot points
        for coord, color, marker in zip(embedding_2d, colors, markers):
            ax.scatter(coord[0], coord[1], c=color, marker=marker, s=10, alpha=0.7)

        # Configure plot
        ax.set_xticks([])
        ax.set_yticks([])
        plt.grid(False)

        # Create legend handles
        unique_genofeatures = sorted(set(id_to_genofeature.values()))
        legend_handles = []
        
        for genofeature in unique_genofeatures:
            if genofeature in genofeature_to_color and genofeature in genofeature_to_marker:
                color = genofeature_to_color[genofeature]
                marker = genofeature_to_marker[genofeature]
                display_name = genofeature_to_display.get(genofeature, genofeature)
                
                handle = mlines.Line2D([], [], 
                                     color=color,
                                     marker=marker,
                                     linestyle='None',
                                     markersize=8,
                                     label=display_name)
                legend_handles.append(handle)

        # Save plot without legend
        plt.savefig(output_file, dpi=600, bbox_inches="tight")
        plt.close()
        
        return legend_handles
        
        
    except Exception as e:
        print(f"Error processing plot: {e}")
        print("Available columns in metadata_df:", metadata_df.columns.tolist())
        return False



# Main execution
print("Loading metadata...")
module_df = pd.read_csv("${filtered_tsv}", sep='\\t')
metadata_df = pd.read_csv("${metadata_file}", sep='\\t')
display_names = dict(zip(metadata_df["genofeature"], metadata_df["display_name"]))

# Generate plots for specified parameters
legend_handles = None
for nn in [75, 100, 125]:
    for md in [0, 0.3, 0.5, 0.7]:
        output_file = f"plots/umap_nn{nn}_md{int(md*10)}.png"
        handles = plot_umap(module_df, nn, md, output_file, display_names)
        if legend_handles is None:
            legend_handles = handles

from PIL import Image
import math

def tile_images(image_dir, output_file, legend_handles, padding=10, bg_color=(255,255,255,255)):
    # Find and sort PNG files
    png_files = [os.path.join(image_dir, f) for f in os.listdir(image_dir) if f.lower().endswith(".png")]
    
    # Extract parameter values using proper filename parsing
    nn_values = set()
    md_values = set()
    for filename in png_files:
        base = os.path.basename(filename)  # Get just the filename
        if base.startswith("umap_nn"):
            parts = base.replace(".png", "").split("_")  # Split on underscore after removing .png
            try:
                nn_part = parts[1]  # gets "nn125"
                nn = int(nn_part.replace("nn", ""))  # properly removes "nn" prefix
                md_part = parts[2]  # gets "md5"
                md = float(md_part.replace("md", "")) / 10
                nn_values.add(nn)
                md_values.add(md)
            except (IndexError, ValueError) as e:
                print(f"Skipping malformed filename: {filename}")
                continue

    # Sort the values
    nn_values = sorted(list(nn_values))
    md_values = sorted(list(md_values))

    def get_params(filename):
        base = os.path.basename(filename)
        parts = base.replace(".png", "").split("_")
        try:
            # Extract everything after "nn" until "_"
            nn_part = parts[1] 
            nn = int(nn_part.replace("nn", ""))  # properly removes "nn" prefix
            
            # Extract everything after "md" until ".png"
            md_part = parts[2]  
            md = float(md_part.replace("md", "")) / 10
            
            return (nn_values.index(nn), md_values.index(md))
        except (IndexError, ValueError) as e:
            print(f"Error parsing filename {filename}: {e}")
            return (0, 0)
    
    png_files.sort(key=get_params)

    # Open all images and get sizes
    images = [Image.open(f) for f in png_files]
    sizes = [img.size for img in images]

    cols = len(md_values)  # Number of columns based on md values
    rows = len(nn_values)  # Number of rows based on nn values

    # Calculate column widths and row heights
    col_widths = [0] * cols
    row_heights = [0] * rows
    for idx, (w, h) in enumerate(sizes):
        row = idx // cols
        col = idx % cols
        col_widths[col] = max(col_widths[col], w)
        row_heights[row] = max(row_heights[row], h)


    label_height = 500  # Height for labels
    label_width = 500  # Width for labels

    axis_title_height = 200 

    grid_width = sum(col_widths) + padding * (cols - 1)
    grid_height = sum(row_heights) + padding * (rows - 1)

    # Create canvas with space for labels
    total_width = grid_width + label_width + padding * 2 + int(grid_height * 0.3)
    total_height = grid_height + label_height * 2 + axis_title_height * 2  # Space for top and bottom labels
    canvas = Image.new('RGBA', (total_width, total_height), bg_color)

    # Create a drawing object
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 120)


    # Add X-axis title
    x_title = "Minimum Distance (md)"
    x_title_img = Image.new('RGBA', (grid_width, axis_title_height), bg_color)
    x_draw = ImageDraw.Draw(x_title_img)
    x_draw.text((grid_width//2, axis_title_height//2), x_title, fill='black', font=font, anchor='mm')
    
    # Paste X-axis title at the top
    x_title_x = label_width
    x_title_y = axis_title_height//2  # Moved higher up
    canvas.paste(x_title_img, (x_title_x, x_title_y))

    # Create and rotate Y-axis title
    y_title = "Nearest Neighbors (nn)"
    y_title_img = Image.new('RGBA', (grid_height, axis_title_height), bg_color)
    y_draw = ImageDraw.Draw(y_title_img)
    y_draw.text((grid_height//2, axis_title_height//2), y_title, fill='black', font=font, anchor='mm')
    y_title_img = y_title_img.rotate(90, expand=True)
    
    # Position Y-axis title to the left
    y_title_x = 50  # Keep it on the far left
    y_title_y = label_height + (grid_height // 2) - (y_title_img.height // 2)  # Center with grid
    canvas.paste(y_title_img, (y_title_x, y_title_y), y_title_img)



    # Add column labels (md values)
    md_labels = [f"md={md:.1f}" for md in md_values]
    for col, label in enumerate(md_labels):
        x = label_width + sum(col_widths[:col]) + padding * col + col_widths[col]//2 
        y = axis_title_height + label_height//2
        draw.text((x, y), label, fill='black', font=font, anchor='mm')

    nn_label_height = 500
    nn_label_width = 200

    # Add row labels (nn values)
    nn_labels = [f"nn={nn}" for nn in nn_values]
    for row, label in enumerate(nn_labels):
        # Create a new image for the rotated text
        label_img = Image.new('RGBA', (nn_label_height, nn_label_width), bg_color)
        label_draw = ImageDraw.Draw(label_img)
        
        # Draw the text vertically
        label_draw.text((nn_label_height//2, nn_label_width//2), label,
                        fill='black', font=font, anchor='mm')
        
        # Rotate the image
        label_img = label_img.rotate(90, expand=True)
        
        # Calculate position
        x = y_title_x + axis_title_height + 100  # Increased spacing after y-axis title
        y = label_height + sum(row_heights[:row]) + padding * row + row_heights[row]//2 - label_img.height//2
        
        # Paste the rotated text
        canvas.paste(label_img, (x, y), label_img)

    # Paste plots into grid with offset for labels
    x_offset = label_width
    y_offset = label_height
    for idx, img in enumerate(images):
        row = idx // cols
        col = idx % cols
        x = x_offset + sum(col_widths[:col]) + padding * col
        y = y_offset + sum(row_heights[:row]) + padding * row
        canvas.paste(img, (x, y))

    # Create and paste legend
    figure_width = 16
    figure_height = 36
    
    legend_fig = plt.figure(figsize=(figure_width, figure_height))
    legend_ax = legend_fig.add_subplot(111)
    legend_ax.axis('off')

    if legend_handles:
        legend = legend_ax.legend(handles=legend_handles,
                                title="Genofeature",
                                loc="center left",
                                bbox_to_anchor=(0, 0.5),
                                frameon=False,
                                fontsize=30,              # Adjusted font size
                                title_fontsize=36,        # Adjusted title size
                                borderaxespad=2,
                                labelspacing=2.0,         # Adjusted spacing
                                handlelength=4,
                                handleheight=4,
                                markerscale=4)

    # Save initial legend
    legend_path = "temp_legend.png"
    legend_fig.savefig(legend_path, 
                      bbox_inches='tight', 
                      dpi=150,              # Reduced DPI to prevent size issues
                      facecolor='white',
                      edgecolor='none',
                      pad_inches=0.5)
    plt.close(legend_fig)

    # Load and scale legend
    legend_img = Image.open(legend_path)
    legend_height = grid_height
    legend_width = int(grid_height * 0.25)  # 25% of grid height
    legend_img = legend_img.resize((legend_width, legend_height), Image.Resampling.LANCZOS)

    # Paste legend
    legend_x = x_offset + grid_width + padding
    legend_y = y_offset
    canvas.paste(legend_img, (legend_x, legend_y))

    # Cleanup and save
    os.remove(legend_path)
    canvas.save(output_file)
    print(f"Tiled image with legend saved as {output_file}")
    #--------END Tiled Image Block----------------
tile_images("plots", "plots/tiled_image.png", legend_handles)
    """
    
}

workflow {
    MODIFY_CLUSTERS()
}