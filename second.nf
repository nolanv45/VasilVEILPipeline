nextflow.enable.dsl=2

process SPLIT_BY_GENOFEATURE {
    publishDir "${params.outdir}",
        mode: 'copy'

    input:
        path filtered_tsv 
        path input_fasta

    output:
        path "files_for_embeddings", emit: fastas_for_embeddings

    script:
    def selected_list = params.selected_genofeatures.collect { "\'${it}\'" }.join(', ')
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import pandas as pd
    import os

    # Create the base output directory
    os.makedirs("files_for_embeddings", exist_ok=True)
    df = pd.read_csv("${filtered_tsv}", sep='\\t')
    seq_dict = SeqIO.to_dict(SeqIO.parse("${input_fasta}", "fasta"))
    

    selected = [gf.lower() for gf in [${selected_list}]]

    # Group by protein and genofeature
    for protein in df['protein'].unique():
        protein_df = df[df['protein'] == protein]
        
        for genofeature in protein_df['genofeature'].unique():
            if genofeature.lower() in selected:
                # Create nested directory structure
                dir_path = os.path.join("files_for_embeddings", protein, genofeature)
                os.makedirs(dir_path, exist_ok=True)
                
                # Filter records for this genofeature
                genofeature_df = protein_df[protein_df['genofeature'] == genofeature]
                
                # Save filtered TSV
                genofeature_df.to_csv(f"{dir_path}/filtered.tsv", sep='\\t', index=False)
                
                # Extract matching sequences
                matching_seqs = []
                for orf_id in genofeature_df['orf_id']:
                    if orf_id in seq_dict:
                        matching_seqs.append(seq_dict[orf_id])
                
                # Write matching sequences to FASTA
                if matching_seqs:
                    SeqIO.write(matching_seqs, f"{dir_path}/{protein}_{genofeature}.fasta", "fasta")
                
                print(f"Processed {protein}/{genofeature}: {len(matching_seqs)} sequences")
    """
}

process EMBEDDINGS {
    publishDir "${params.outdir}",
        mode: 'copy'
        
    label 'gpu'  
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    memory { 50.GB }
    
    input:
        path split_dir
        
    output:
        path "*", emit: embeddings_dirs
        
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    mkdir -p embeddings

    extract_script="/mnt/VEIL/tools/embeddings/models/extract.py"
    model_path="/mnt/VEIL/tools/embeddings/models/esm2_t36_3B_UR50D.pt"
    ls -R "${split_dir}"

    # Process each protein directory
    for protein_dir in "${split_dir}"/*; do
        if [ -d "\$protein_dir" ]; then
            protein=\$(basename "\$protein_dir")
            echo "Processing protein directory: \$protein"
            mkdir -p "embeddings/\$protein"
            
            # Process genofeature directories
            for genofeature_dir in "\$protein_dir"/*; do
                if [ -d "\$genofeature_dir" ]; then
                    genofeature=\$(basename "\$genofeature_dir")
                    echo "Processing genofeature directory: \$genofeature"
                    mkdir -p "embeddings/\$protein/\$genofeature"
                    
                    # Process FASTA files
                    for fasta in "\$genofeature_dir"/*.fasta; do
                        if [ -f "\$fasta" ]; then
                            echo "Processing FASTA: \$fasta"
                            output_dir="embeddings/\$protein/\$genofeature"
                            
                            /home/nolanv/.conda/envs/esm-umap/bin/python \\
                                "\$extract_script" \\
                                "\$model_path" \\
                                "\$fasta" \\
                                "\$output_dir" \\
                                --repr_layers 36 \\
                                --include mean || exit 1
                        fi
                    done
                fi
            done
        fi
    done
    """
}


process HDBSCAN {
    publishDir "${params.outdir}/hdbscan",
        mode: 'copy'
        
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    input:
        path coordinates_dir  // Directory containing the pre-generated coordinates
        path filtered_tsv    // TSV file with metadata
        path metadata_file   // Metadata file with colors/markers
        
    output:
        path "plots/*.{png,svg}", emit: plots     // Plot outputs

    script:
    """
#!/usr/bin/env python3
import os
import torch
import numpy as np
import matplotlib.pyplot as plt
import umap
import pandas as pd
import hdbscan
import random 

#----------v_25_06_25---------------------------------
# refers to plotted coordinates from 04a umap script to cluster and map providing consistency and efficiency

# Define parameter ranges
umap_nn = 100
umap_md = 0.7
hdbscan_nn_list = [100]
hdbscan_min_dist_list = [0]
hdbscan_min_cluster_size_list = [10, 20, 30, 40] #usually 10, 20, 30, 40

# Create required directories
os.makedirs("plots", exist_ok=True)
os.makedirs("plots/clusters", exist_ok=True)

def plot_umap_hdbscan(module_df, metadata_df, nn, md, mc, output_path, iteration_output):
    try:
        md_int = int(md * 10)
        coord_file = os.path.join("${coordinates_dir}", f"coords/coordinates_nn{nn}_md{md_int}.npy")
        id_file = os.path.join("${coordinates_dir}", f"ids/embedding_ids_nn{nn}_md{md_int}.txt")

        embedding_2d = np.load(coord_file)
        with open(id_file, 'r') as f:
            embedding_ids = [line.strip() for line in f]

        print(f"Loaded coordinates from {coord_file}")
        print(f"Shape: {embedding_2d.shape}, IDs: {len(embedding_ids)}")

        module_df["normalized_orf_id"] = module_df["orf_id"].str.replace(".", "-", regex=False)
        id_to_genofeature = dict(zip(module_df["normalized_orf_id"], module_df["genofeature"]))
        genofeature_to_color = dict(zip(metadata_df["genofeature"], metadata_df["color"]))
        genofeature_to_marker = dict(zip(metadata_df["genofeature"], metadata_df["marker"]))
        genofeature_to_display = dict(zip(metadata_df["genofeature"], metadata_df["display_name"]))
      
        labels = hdbscan.HDBSCAN(min_samples=2, min_cluster_size=mc).fit_predict(embedding_2d)
        clustered = (labels >= 0)

        cluster_df = pd.DataFrame({
            'embedding_id': embedding_ids,
            'cluster_label': labels,
            'genofeature': [id_to_genofeature.get(eid, "unknown") for eid in embedding_ids]
        })

        cluster_output_path = os.path.join(iteration_output,
                                         f"hdbscan_nn{nn}_md{md_int}_minclust{mc}_clusters.csv")
        cluster_df.to_csv(cluster_output_path, index=False)
        print(f"Cluster membership saved: {cluster_output_path}")

        # Create plot
        fig, ax = plt.subplots(figsize=(10, 7))

        # Plot clustered points
        scatter = ax.scatter(embedding_2d[clustered, 0], embedding_2d[clustered, 1],
                           c=labels[clustered], s=10, alpha=0.7, cmap="Spectral")

        # Plot unclustered points
        ax.scatter(embedding_2d[~clustered, 0], embedding_2d[~clustered, 1],
                  color="gray", s=10, alpha=0.5, label="not clustered")

        # Label cluster centers
        unique_labels = np.unique(labels[labels >= 0])
        for cluster in unique_labels:
            cluster_points = embedding_2d[labels == cluster]
            center_x, center_y = cluster_points[:, 0].mean(), cluster_points[:, 1].mean()
            plt.text(center_x, center_y, str(cluster),
                    fontsize=10, weight='regular',
                    ha='center', va='center')

        # Configure plot
        ax.set_xticks([])
        ax.set_yticks([])
        plt.grid(False)

        plt.title(f"UMAP (nn={nn}, md={md}) with HDBSCAN (minclust={mc})")
        plt.colorbar(scatter, label="Cluster Labels")
        
        # Save plot
        plt.tight_layout()
        plt.savefig(output_path, dpi=600, bbox_inches="tight")
        print(f"Saved plot: {output_path}")
        
        plt.close()
        return True

    except Exception as e:
        print(f"Error processing plot: {e}")
        print("Debug info - metadata columns:", metadata_df.columns.tolist())
        return False

# Main execution
print("Loading metadata...")
module_df = pd.read_csv("${filtered_tsv}", sep='\t')
metadata_df = pd.read_csv("${metadata_file}", sep='\t')

# Generate plots
for nn in [100]:
    for md in [0.7]:
        for mc in [10, 20, 30, 40]:
            output_path = os.path.join("plots", f"hdbscan_nn{nn}_md{int(md*10)}_minclust{mc}.png")
            if plot_umap_hdbscan(module_df, metadata_df, nn, md, mc, output_path, "plots"):
                print(f"Successfully generated plot for nn={nn}, md={md}, mc={mc}")
            else:
                print(f"Failed to generate plot for nn={nn}, md={md}, mc={mc}")

    """

}



process UMAP_PROJECTION {
    publishDir "${params.outdir}/umap",
        mode: 'copy'
        
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    input:
        path coordinates_dir  // Directory containing the pre-generated coordinates
        path filtered_tsv    // TSV file with metadata
        path metadata_file   // Metadata file with colors/markers
        
    output:
        path "plots/*.{png,svg}", emit: plots     // Plot outputs
        
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

        # Create plot
        fig, ax = plt.subplots(figsize=(10, 7))

        # Draw connections first (background)
        for indices in groups.values():
            if len(indices) > 1:
                for i in range(len(indices) - 1):
                    for j in range(i + 1, len(indices)):
                        idx1, idx2 = indices[i], indices[j]
                        x1, y1 = embedding_2d[idx1, 0], embedding_2d[idx1, 1]
                        x2, y2 = embedding_2d[idx2, 0], embedding_2d[idx2, 1]
                        ax.plot([x1, x2], [y1, y2], color='#CCCCCC', alpha=0.5, 
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

        # Add legend to main plot
        if legend_handles:
            legend = ax.legend(handles=legend_handles,
                             title="Genofeature",
                             bbox_to_anchor=(1.05, 1),
                             loc='upper left',
                             frameon=False)
        #save main plot with legend
        plt.tight_layout()
        plt.savefig(output_file, dpi=600, bbox_inches="tight")
        print(f"Saved plot with legend: {output_file}")
        
        # Create separate legend figure
        legend_fig = plt.figure(figsize=(4, len(legend_handles) * 0.3))
        legend_ax = legend_fig.add_subplot(111)
        legend_ax.axis('off')
        
        # Add legend to separate figure
        legend_ax.legend(handles=legend_handles,
                        title="Genofeature",
                        loc="center left",
                        frameon=False)
        
        # Close figures
        plt.close(fig)
        plt.close(legend_fig)
        return True
        
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
for nn in [75, 100, 125]:
    for md in [0, 0.3, 0.5, 0.7]:
        output_file = f"plots/umap_nn{nn}_md{int(md*10)}.png"
        if plot_umap(module_df, nn, md, output_file, display_names):
            print(f"Successfully generated plot for nn={nn}, md={md}")
        else:
            print(f"Failed to generate plot for nn={nn}, md={md}")
from PIL import Image
import math

def tile_images(image_dir, output_file, padding=10, bg_color=(255,255,255,255)):
    # Find all PNG files in the output directory
    png_files = [os.path.join(image_dir, f) for f in os.listdir(image_dir) if f.lower().endswith(".png")]
    png_files.sort()
    if not png_files:
        print("No PNG files found to tile.")
        return

    # Open all images and get their sizes
    images = [Image.open(f) for f in png_files]
    sizes = [img.size for img in images]

    n = len(images)
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)

    # Compute max width for each column and max height for each row
    col_widths = [0] * cols
    row_heights = [0] * rows
    for idx, (w, h) in enumerate(sizes):
        row = idx // cols
        col = idx % cols
        if w > col_widths[col]:
            col_widths[col] = w
        if h > row_heights[row]:
            row_heights[row] = h

    # Compute x offsets for columns and y offsets for rows
    x_offsets = [0]
    for w in col_widths[:-1]:
        x_offsets.append(x_offsets[-1] + w + padding)
    y_offsets = [0]
    for h in row_heights[:-1]:
        y_offsets.append(y_offsets[-1] + h + padding)

    canvas_width = sum(col_widths) + padding * (cols - 1)
    canvas_height = sum(row_heights) + padding * (rows - 1)
    canvas = Image.new('RGBA', (canvas_width, canvas_height), bg_color)

    for idx, img in enumerate(images):
        row = idx // cols
        col = idx % cols
        x = x_offsets[col]
        y = y_offsets[row]
        canvas.paste(img, (x, y))

    canvas.save(output_file)
    print(f"Tiled image saved as {output_file}")
    #--------END Tiled Image Block----------------
tile_images("plots", "plots/tiled_image.png")
    """
}

process GENERATE_COORDINATES {
    publishDir "${params.outdir}/coordinates",
        mode: 'copy'
        
    label 'gpu'  
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    memory { 50.GB }
    
    input:
        path(embeddings)   
        path(filtered_tsv)
        
    output:
        path "*", emit: coordinates_files
        
    script:
    def selected_list = params.selected_genofeatures.collect { "\'${it}\'" }.join(', ')
    """
    #!/usr/bin/env python3
    import os
    import torch
    import numpy as np
    import umap
    import random
    from pathlib import Path

    # Create output directories
    os.makedirs("coords", exist_ok=True)
    os.makedirs("ids", exist_ok=True)

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
        
        # Generate UMAP for each parameter combination
        for nn in [75, 100, 125]:
            for md in [0, 0.3, 0.5, 0.7]:
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
                coord_file = f"coords/coordinates_nn{nn}_md{md_int}.npy"
                id_file = f"ids/embedding_ids_nn{nn}_md{md_int}.txt"
                
                np.save(coord_file, coordinates)
                with open(id_file, 'w') as f:
                    for eid in embedding_ids:
                        f.write(f"{eid}\\n")
                print(f"Saved coordinates for nn={nn}, md={md}")
    else:
        print("No embeddings found to process")
        # Create empty files to satisfy output requirements
        np.save("coords/empty_coordinates.npy", np.array([]))
        with open("ids/empty_ids.txt", "w") as f:
            pass

    """
}

process CLUSTER {
    publishDir "${params.outdir}/clusters",
        mode: 'copy'
        
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    input:
        path coordinates_dir  // Directory containing the pre-generated coordinates
        path filtered_tsv    // TSV file with metadata
        path metadata_file   // Metadata file with colors/markers
        
    output:
        path ""

    script:
    """
import os
import pandas as pd

# Original Version from Barb Ferrell 05/2025
# Alicia Holk aholk@udel.edu

#nn 100 and clusters 30 and 10

# Hardcoded paths and parameters
input_dir = "/mnt/VEIL/users/aholk/temphd/05a_phidra_pola_polb_rnr_helicase_nn100_mindist07_hdbscan/"
file1 = os.path.join(input_dir, "hdbscan_nn100_umapmd7_md0_minclust30_clusters.csv")
file2 = os.path.join(input_dir, "hdbscan_nn100_umapmd7_md0_minclust10_clusters.csv")
output_dir = "/mnt/VEIL/users/aholk/temphd/05b_phidra_pola_polb_rnr_helicase_cluster_assignment/"
os.makedirs(output_dir, exist_ok=True)

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
df1_sorted.to_csv(os.path.join(output_dir, "hdbscan_modified_cluster_labels.tsv"), index=False, sep="\t")
    """
}

process UMAP_CLUSTER{
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
input_dir = "/mnt/VEIL/users/aholk/temphd/03d_filtered_dnab_embeddings/"  # Updated to filtered embeddings directory
coords_dir = "/mnt/VEIL/users/aholk/temphd/04a_umap/"  # Directory with saved coordinates from 04a
metadata_file = "/mnt/VEIL/users/aholk/temphd/05b_phidra_pola_polb_rnr_helicase_cluster_assignment/hdbscan_modified_cluster_labels.tsv"
output_dir = "/mnt/VEIL/users/aholk/temphd/05c_phidra_pola_polb_rnr_helicase_umap_w_cluster_assignment/"
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


workflow {
    // Create input channels
    ch_filtered_tsv = Channel.fromPath(params.second_run)
    ch_metadata = Channel.fromPath(params.genofeature_metadata)
    
    // Generate coordinates
    ch_split_dir = SPLIT_BY_GENOFEATURE(
        ch_filtered_tsv,
        params.input_fasta
    )

    ch_embeddings = EMBEDDINGS(
        ch_split_dir
    )

    ch_coordinates = GENERATE_COORDINATES(
        ch_embeddings,
        ch_filtered_tsv
    )

    // Use duplicated channels in parallel processes
    UMAP_PROJECTION(
        ch_coordinates,
        ch_filtered_tsv,
        ch_metadata
    )

    HDBSCAN(
        ch_coordinates,
        ch_filtered_tsv,
        ch_metadata
    )
}
