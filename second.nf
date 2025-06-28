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

    # Group by protein and genofeature
    for protein in df['protein'].unique():
        protein_df = df[df['protein'] == protein]
        
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


// process HDBSCAN {
//     publishDir "${params.outdir}/hdbscan",
//         mode: 'copy'
        
//     conda "/home/nolanv/.conda/envs/esm-umap"
    
//     input:
//         path(umap_embeddings)
        
//     output:
//         path "*", emit: hdbscan_results
        
//     script:


// }



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


// Below is tile stuff i think
// # Use this block to tile images at the end of your script into a single image
// from PIL import Image
// import math

// def tile_images(image_dir, output_file, padding=10, bg_color=(255,255,255,255)):
//     # Find all PNG files in the output directory
//     png_files = [os.path.join(image_dir, f) for f in os.listdir(image_dir) if f.lower().endswith(".png")]
//     png_files.sort()
//     if not png_files:
//         print("No PNG files found to tile.")
//         return

//     # Open all images and get their sizes
//     images = [Image.open(f) for f in png_files]
//     sizes = [img.size for img in images]

//     n = len(images)
//     cols = math.ceil(math.sqrt(n))
//     rows = math.ceil(n / cols)

//     # Compute max width for each column and max height for each row
//     col_widths = [0] * cols
//     row_heights = [0] * rows
//     for idx, (w, h) in enumerate(sizes):
//         row = idx // cols
//         col = idx % cols
//         if w > col_widths[col]:
//             col_widths[col] = w
//         if h > row_heights[row]:
//             row_heights[row] = h

//     # Compute x offsets for columns and y offsets for rows
//     x_offsets = [0]
//     for w in col_widths[:-1]:
//         x_offsets.append(x_offsets[-1] + w + padding)
//     y_offsets = [0]
//     for h in row_heights[:-1]:
//         y_offsets.append(y_offsets[-1] + h + padding)

//     canvas_width = sum(col_widths) + padding * (cols - 1)
//     canvas_height = sum(row_heights) + padding * (rows - 1)
//     canvas = Image.new('RGBA', (canvas_width, canvas_height), bg_color)

//     for idx, img in enumerate(images):
//         row = idx // cols
//         col = idx % cols
//         x = x_offsets[col]
//         y = y_offsets[row]
//         canvas.paste(img, (x, y))

//     canvas.save(output_file)
//     print(f"Tiled image saved as {output_file}")
//     #--------END Tiled Image Block----------------

// if __name__ == "__main__":
//     # Load data
//     # to process all subdirectories:
//     # embeddings, embedding_ids = load_all_embeddings(input_dir)
//     # to process specific subdirectories, pass them as a list:
//     selected_subdirs = [${selected_list}]
//     embeddings, embedding_ids = load_all_embeddings(input_dir, selected_subdirs)
//     print(f"Loaded {len(embeddings)} embeddings.")
//     print("First 5 embedding IDs:", embedding_ids[:5])
    
//     module_df = load_module(module_file)
//     color_map, display_names, marker_map = load_metadata(metadata_file)
//     # print(color_map)
//     # print(display_names)
//     # print(marker_map)

//     # Check if embedding_ids are in the module_df
//     not_found_rows = []
//     normalized_orf_ids = set(module_df["normalized_orf_id"])
//     for eid in embedding_ids:
//         if eid not in normalized_orf_ids:
//             not_found_rows.append({"embedding_id": eid})

//     if not_found_rows:
//         not_found_df = pd.DataFrame(not_found_rows)
//         not_found_path = os.path.join(output_dir, "embeddings_not_found_in_df.tsv")
//         not_found_df.to_csv(not_found_path, sep="\t", index=False)
//         print(f"Embeddings not found in df saved to: {not_found_path}")
//     else:
//         print("All embedding_ids were found in the dataframe.")

//     # Check if module_df is not None before proceeding
//     if module_df is not None:
//         # Map colors to metadata
//         module_df['manual_color'] = module_df[column_name].map(color_map)
//         module_df['display_names'] = module_df[column_name].map(display_names)
//         module_df['manual_marker'] = module_df[column_name].map(marker_map)
//         print("Assigned Colors:", color_map)
//         print("Assigned Labels:", display_names)
//         print("Assigned Shapes:", marker_map)

//     print(module_df.head())

// #if we have emeddings in our list
//     if len(embeddings) > 0:
//         for nn in n_neighbors_list:
//             for md in min_dist_list:
//                 md_int = int(md * 10)  # Convert min_dist to integer for filename
//                 output_path = os.path.join(output_dir, f"umap_nn{nn}_md{md_int}.png") #create output path to store it
//                 plot_umap(embeddings, embedding_ids, module_df, nn, md, output_path, display_names) #run the plot_umap function
//                 print(f"Saved: {output_path}")
//     else:
//         print("No valid embeddings found.")

//     # Tile all UMAP images into one final image
//     tiled_output = os.path.join(output_dir, "umap_tiled_all.png")
//     """
// }




workflow {
    // Create input channels
    ch_filtered_tsv = Channel.fromPath(params.second_run)
    ch_metadata = Channel.fromPath(params.genofeature_metadata)
    
    // Generate coordinates
    // EMBEDDINGS("/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/output_test/files_for_embeddings/")
    

    // GENERATE_COORDINATES(
    //     "/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/output_test/embeddings",
    //     ch_filtered_tsv
    // )

    UMAP_PROJECTION(
        "/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/output_test/coordinates",
        ch_filtered_tsv,
        ch_metadata
    )
}