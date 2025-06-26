nextflow.enable.dsl=2

process SPLIT_BY_GENOFEATURE {
    publishDir "${params.outdir}/files_for_embeddings",
        mode: 'copy'

    input:
        tuple path(filtered_tsv), path(input_fasta)

    output:
        path "*", emit: fastas_for_embeddings

    script:
    def selected_list = params.selected_genofeatures.collect { "\'${it}\'" }.join(', ')
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import pandas as pd
    import os

    df = pd.read_csv("${filtered_tsv}", sep='\\t')
    seq_dict = SeqIO.to_dict(SeqIO.parse("${input_fasta}", "fasta"))
    

    selected = [gf.lower() for gf in [${selected_list}]]

    # Group by Protein and Genofeature
    for protein in df['Protein'].unique():
        protein_df = df[df['Protein'] == protein]
        
        for genofeature in protein_df['genofeature'].unique():
            if genofeature.lower() in selected:
                # Create nested directory structure
                dir_path = os.path.join(protein, genofeature)
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
    publishDir "${params.outdir}/embeddings",
        mode: 'copy'
        
    label 'gpu'  
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    memory { 50.GB }
    
    input:
        path split_dir
        
    output:
        path "*/*", emit: embeddings_dirs
        
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    extract_script="/mnt/VEIL/tools/embeddings/models/extract.py"
    model_path="/mnt/VEIL/tools/embeddings/models/esm2_t36_3B_UR50D.pt"

    for protein_dir in */; do
        if [ -d "\$protein_dir" ]; then
            echo "Processing protein directory: \$protein_dir"
            protein=\$(basename "\$protein_dir")
            
            for genofeature_dir in "\${protein_dir}"*/; do
                if [ -d "\$genofeature_dir" ]; then
                    echo "Processing genofeature directory: \$genofeature_dir"
                    genofeature=\$(basename "\$genofeature_dir")
                    
                    for fasta in "\${genofeature_dir}"*.fasta; do
                        if [ -f "\$fasta" ]; then
                            # Create final output directory structure
                            output_dir="\${protein}/\${genofeature}"
                            mkdir -p "\$output_dir"
                            
                            echo "Processing FASTA: \$fasta -> \$output_dir"
                            
                            # ESM embeddings 
                            /home/nolanv/.conda/envs/esm-umap/bin/python \\
                                "\$extract_script" \\
                                "\$model_path" \\
                                "\$fasta" \\
                                "\$output_dir" \\
                                --repr_layers 36 \\
                                --include mean || exit 1

                            # Move batch directories and clean up
                            if [ -d "\$output_dir/batch_0" ]; then
                                # Keep only batch directories
                                find "\$output_dir" -maxdepth 1 -type f -delete
                            fi
                        fi
                    done
                fi
            done
        fi
    done
    """
}

process UMAP {
    publishDir "${params.outdir}/umap",
        mode: 'copy'
        
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    
    input:
        path(embeddings)   
        path(filtered_tsv)
        
    output:
        path "*", emit: umap_embeddings
        
    script:
    def selected_list = params.selected_genofeatures.collect { "\'${it}\'" }.join(', ')
    """
#!/usr/bin/env python3
import os
import torch
import numpy as np
import matplotlib.pyplot as plt
import umap # modified from umap.umap_ to umap
import pandas as pd
import matplotlib.lines as mlines
import random  # <-- Add this import
from pathlib import Path


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
n_neighbors_list = [75, 100, 125]
min_dist_list = [0.3, 0.5, 0.7]

# Hardcoded paths and parameters
input_dir = "${embeddings}"  # Directory containing subdirectories with .pt files
module_file = "${filtered_tsv}"  # Path to module file (TSV)
metadata_file = "${params.genofeature_metadata}"  # Path to metadata file
column_name = "genofeature"

output_dir = "${params.outdir}/umap"  # Output directory for UMAP plots
os.makedirs(output_dir, exist_ok=True)


def load_embedding(file_path):
    # Loads embeddings from a .pt file.
    try:
        embedding_data = torch.load(file_path)
        
        #grab mean representation from layer 36 (how the data structure is laid out so it represents the whole protein)
        embedding = embedding_data["mean_representations"][36].numpy()
        embedding_id = embedding_data.get("label")
        if embedding_id is None:
            raise ValueError(f"Embedding ID not found in file: {os.path.basename(file_path)}")
        return embedding, embedding_id
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return None, None

def load_module(module_file):
    # Loads module data from a TSV file.
    try:
        module_df = pd.read_csv(module_file, sep='\t')
        module_df["normalized_orf_id"] = module_df["orf_id"].str.replace(".", "-", regex=False)
        return module_df
    except Exception as e:
        print(f"Error loading module file: {e}")
        return None


def load_metadata(metadata_file):
    # Loads metadata (including colors and shapes) from a file into dictionaries.
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


def find_pt_files(base_dir):
    pt_files = []
    # Walk through all directories
    for root, _, files in os.walk(base_dir):
        # Only look in batch_0 directories
        if 'batch_0' in root:
            for file in files:
                if file.endswith('.pt'):
                    pt_files.append(os.path.join(root, file))
    return pt_files





def load_all_embeddings(base_dir, subdirs=None):
    embeddings = []
    embedding_ids = []
    
    print(f"Looking for embeddings in: {base_dir}")
    print(f"Selected genofeatures: {subdirs}")
    
    # Find all .pt files
    pt_files = find_pt_files(base_dir)
    print(f"Found {len(pt_files)} .pt files")
    
    # Process each .pt file
    for file_path in pt_files:
        # Check if file's parent directory matches any of our selected genofeatures
        parent_dir = Path(file_path).parent.parent.name
        if subdirs and parent_dir not in subdirs:
            continue
            
        print(f"Processing: {file_path}")
        embedding, embedding_id = load_embedding(file_path)
        if embedding is not None and embedding_id is not None:
            embeddings.append(embedding)
            embedding_ids.append(embedding_id)
            print(f"Successfully loaded embedding from {file_path}")
    
    return embeddings, embedding_ids


def plot_umap(embeddings, embedding_ids, module_df, nn, md, output_file):
    #  Plots UMAP embedding with colors mapped from metadata and adds a legend with both colors and marker shapes.
    # print(nn, md)
    # Normalize orf_name for consistent lookup (removing, ".", etc)
    module_df["normalized_orf_id"] = module_df["orf_id"].str.replace(".", "-", regex=False)

    # Create lookup dictionaries for colors and markers
    id_to_color = dict(zip(module_df["normalized_orf_id"], module_df["manual_color"]))
    id_to_marker = dict(zip(module_df["normalized_orf_id"], module_df["manual_marker"]))

    # Assign colors and markers to embedding IDs (if we can't find color - make it gray, if we can't find marker, make a .)
    colors = [id_to_color.get(embedding_id, "#808080") for embedding_id in embedding_ids]  # Default: gray
    # markers = [id_to_marker.get(embedding_id, ".") for embedding_id in embedding_ids]  # Default: 'o'
    markers = [id_to_marker.get(embedding_id, ".") if pd.notna(id_to_marker.get(embedding_id)) else '.' for embedding_id in embedding_ids] #if NaN or not in list, put a "." as marker style

    # Fit the UMAP model (store data in variable called mapper, call function umap, run with nn and min dist, random state =42 for consistency from run to run so if you rerun it runs the same random order again, n_components = 2 axes, based on euclidean distance and fit transform embedddings)
    # mapper = umap.UMAP(n_neighbors=nn, min_dist=md, random_state=42, n_components=2, metric="euclidean").fit_transform(embeddings)
    
    # Perform UMAP
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=nn,
        min_dist=md,
        metric='cosine',      # or "euclidean" if you prefer
        random_state=42,      # very important
        transform_seed=42,    # very important
        n_epochs=200,         # better drawing optimization but can sometimes skew data if set too high on smaller samples.
    )
    mapper = reducer.fit_transform(embeddings)

#Visualize
    # Create figure and scatter plot
    fig, ax = plt.subplots(figsize=(10, 7))

    for embedding, color, marker in zip(mapper, colors, markers):
        ax.scatter(embedding[0], embedding[1], c=color, marker=marker, s=10, alpha=0.5)#alpha = transparency of points

    # Drop duplicates based on "genofeature"
    unique_features = module_df.drop_duplicates(subset=["genofeature"], keep="first")

#Legend
    # Ensure only genofeature keys that exist in display_names are kept
    # unassigned features are not included in legend 
    existing_keys = unique_features["genofeature"].isin(display_names.keys())
    unique_features = unique_features[existing_keys]

    # Add display names from display_names dictionary
    unique_features["display_name"] = unique_features["genofeature"].map(display_names)

    # Ensure the order in the legend follows display_names keys order
    unique_features["genofeature"] = pd.Categorical(
        unique_features["genofeature"],
        categories=list(display_names.keys()),
        ordered=True
    )

    # Sort the dataframe by genofeature for correct ordering in the legend
    unique_features = unique_features.sort_values("genofeature")

    # Create legend handles using the display_name for the labels
    legend_handles = [
        mlines.Line2D([], [], color=row["manual_color"], marker=row["manual_marker"], linestyle="None", markersize=8,
                    label=row["display_name"])
        for _, row in unique_features.iterrows()
    ]

    # Add legend outside the plot
    ax.legend(handles=legend_handles, title="Genofeature", bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    plt.title(f"UMAP Projection (nn={nn}, md={md})")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")  # Save with space for legend
    plt.close()

# Use this block to tile images at the end of your script into a single image
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

if __name__ == "__main__":
    # Load data
    # to process all subdirectories:
    # embeddings, embedding_ids = load_all_embeddings(input_dir)
    # to process specific subdirectories, pass them as a list:
    selected_subdirs = [${selected_list}]
    embeddings, embedding_ids = load_all_embeddings(input_dir, selected_subdirs)
    print(f"Loaded {len(embeddings)} embeddings.")
    print("First 5 embedding IDs:", embedding_ids[:5])
    
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

#if we have emeddings in our list
    if len(embeddings) > 0:
        for nn in n_neighbors_list:
            for md in min_dist_list:
                md_int = int(md * 10)  # Convert min_dist to integer for filename
                output_path = os.path.join(output_dir, f"umap_nn{nn}_md{md_int}.png") #create output path to store it
                plot_umap(embeddings, embedding_ids, module_df, nn, md, output_path) #run the plot_umap function
                print(f"Saved: {output_path}")
    else:
        print("No valid embeddings found.")

    # Tile all UMAP images into one final image
    tiled_output = os.path.join(output_dir, "umap_tiled_all.png")
    tile_images(output_dir, tiled_output)
    """
}


workflow {
    // Create input channels
    ch_filtered_tsv = Channel.fromPath(params.second_run)
    ch_input_fasta = Channel.fromPath(params.input_fasta)
    
    // // Combine the input files into a tuple
    // ch_inputs = ch_filtered_tsv.combine(ch_input_fasta)

    // ch_split_genofeatures = SPLIT_BY_GENOFEATURE(ch_inputs)

    // Pass to UMAP process
    UMAP("/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/output_test/embeddings", ch_filtered_tsv)

}