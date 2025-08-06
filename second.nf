nextflow.enable.dsl=2

process EMBEDDINGS {
    publishDir "${params.outdir}",
        mode: 'copy'
        
    label 'gpu'  
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    memory { 50.GB }
    
    input:
        val embedding_datasets  // list of datasets to process
        
    output:
        path "*", emit: embeddings_dirs
        
    script:
    def selected_list = params.selected_genofeatures.collect { "\'${it}\'" }.join(', ')
    def datasets_str = embedding_datasets.collect { "\"${it}\"" }.join(' ')
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    mkdir -p embeddings

    extract_script="/mnt/VEIL/tools/embeddings/models/extract.py"
    model_path="/mnt/VEIL/tools/embeddings/models/esm2_t36_3B_UR50D.pt"

    # Process each protein directory
    for dataset in ${datasets_str}; do
        dataset_path="${params.outdir}/\$dataset/"
        embedding_path="\$dataset_path/files_for_embeddings"

        for protein_dir in "\$embedding_path"/*; do
            if [ -d "\$protein_dir" ]; then
                protein=\$(basename "\$protein_dir")
                echo "Processing protein directory: \$protein"
                mkdir -p "embeddings/\$dataset/\$protein"

                # Process genofeature directories, but only if they're in selected_list
                for genofeature_dir in "\$protein_dir"/*; do
                    if [ -d "\$genofeature_dir" ]; then
                        genofeature=\$(basename "\$genofeature_dir")
                        # Check if genofeature is in selected_list (case-insensitive)
                        if echo "${selected_list}" | tr '[:upper:]' '[:lower:]' | grep -iq "\'\$(echo \$genofeature | tr '[:upper:]' '[:lower:]')\'"; then
                            echo "Processing genofeature directory: \$genofeature"
                            output_dir="embeddings/\$dataset/\$protein/\$genofeature"
                            mkdir -p "\$output_dir"
                
                            # Process FASTA files
                            for fasta in "\$genofeature_dir"/*.fasta; do
                                if [ -f "\$fasta" ]; then
                                    fasta_base=\$(basename "\$fasta" .fasta)
                                    
                                    # Check if embeddings exist
                                    if [ -d "${params.outdir}/embeddings/\$dataset/\$protein/\$genofeature/batch_0" ] && [ -n "\$(find "${params.outdir}/embeddings/\$dataset/\$protein/\$genofeature/batch_0" -name '*.pt' 2>/dev/null)" ]; then
                                        echo "Using existing embeddings for \$fasta_base"
                                        cp -r "${params.outdir}/embeddings/\$dataset/\$protein/\$genofeature" "embeddings/\$dataset/\$protein/"
                                    else
                                        echo "Generating embeddings for \$fasta_base"
                                        /home/nolanv/.conda/envs/esm-umap/bin/python \\
                                            "\$extract_script" \\
                                            "\$model_path" \\
                                            "\$fasta" \\
                                            "\$output_dir" \\
                                            --repr_layers 36 \\
                                            --include mean || exit 1
                                    fi
                                fi
                            done
                        fi
                    fi
                done
            fi
        done
    done
    """
}


process HDBSCAN {
    publishDir "${params.outdir}/hdbscan",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.startsWith("plots/tiled_image_md")) filename
            // else if (filename.endsWith(".png")) "plots/${filename}"
            else if (filename.endsWith(".csv")) "clusters_csv/${filename}"
            else null
        }
        
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    input:
        path coordinates_dir  // Directory containing the pre-generated coordinates
        path filtered_tsv    // TSV file with metadata
        path metadata_file   // Metadata file with colors/markers
        path plots
        
    output:
        path "plots/*.png", emit: plots
        path "plots/tiled_image_md*", emit: tiled_image
        path "plots/*.csv", emit: clusters_csv

    script:
    """
#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import hdbscan
import random 
from PIL import Image, ImageDraw, ImageFont
Image.MAX_IMAGE_PIXELS = None 

#----------v_25_06_25---------------------------------
# refers to plotted coordinates from 04a umap script to cluster and map providing consistency and efficiency


# Tile all plots together
def load_and_resize_image(filepath, max_width=2400):
    try:
        with Image.open(filepath) as img:
            # Calculate new dimensions maintaining aspect ratio
            ratio = max_width / img.size[0]
            new_height = int(img.size[1] * ratio)
            
            # Convert to RGB if needed and resize
            if img.mode in ('RGBA', 'P'):
                img = img.convert('RGB')
            
            return img.resize((max_width, new_height), Image.Resampling.LANCZOS)
    except Exception as e:
        print(f"Error loading image {filepath}: {e}")
        return None

def plot_umap_hdbscan(module_df, metadata_df, nn, md, mc, output_path, iteration_output):
    try:
        md_int = int(md * 10)
        coord_file = os.path.join("${coordinates_dir}", f"coords/coordinates_nn{nn}_md{md_int}.tsv")
        coord_df = pd.read_csv(coord_file, sep='\t')
        embedding_ids = coord_df['embedding_id'].tolist()
        embedding_2d = coord_df[['x', 'y']].values

        print(f"Loaded coordinates from {coord_file}")
        print(f"Shape: {embedding_2d.shape}, IDs: {len(embedding_ids)}")

        module_df["normalized_orf_id"] = module_df["orf_id"].str.replace(".", "-", regex=False)

        # id_to_genofeature = dict(zip(module_df["normalized_orf_id"], module_df["genofeature"]))
        # genofeature_to_color = dict(zip(metadata_df["genofeature"], metadata_df["color"]))
        # genofeature_to_marker = dict(zip(metadata_df["genofeature"], metadata_df["marker"]))
        # genofeature_to_display = dict(zip(metadata_df["genofeature"], metadata_df["display_name"]))
      
        labels = hdbscan.HDBSCAN(min_samples=2, min_cluster_size=mc).fit_predict(embedding_2d)
        clustered = (labels >= 0)

        cluster_df = pd.DataFrame({
            'embedding_id': embedding_ids,
            'cluster_label': labels,
            #'genofeature': [id_to_genofeature.get(eid, "unknown") for eid in embedding_ids]
            'genofeature': [module_df.loc[module_df["normalized_orf_id"] == eid, "genofeature"].iloc[0] 
                if eid in module_df["normalized_orf_id"].values else "unknown" 
                for eid in embedding_ids]
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


# Create required directories
os.makedirs("plots", exist_ok=True)

# Copy input UMAP plots to working directory
os.system("cp ${plots} plots/")

# Generate HDBSCAN plots for each UMAP plot
for umap_file in os.listdir("plots"):
    if umap_file.startswith("umap_nn"):
        # Extract parameters from filename
        parts = umap_file.replace(".png", "").split("_")
        nn = int(parts[1].replace("nn", ""))
        md = float(parts[2].replace("md", "")) / 10
        md_int = int(md * 10)
        
        # Load corresponding coordinates
        # coord_file = os.path.join("${coordinates_dir}", f"coords/coordinates_nn{nn}_md{md_int}.tsv")
        # coord_df = pd.read_csv(coord_file, sep='\t')
        
        # Generate HDBSCAN plots with different min_cluster_size values
        for mc in ${params.mc}:
            output_path = f"plots/hdbscan_nn{nn}_md{md_int}_minclust{mc}.png"
            plot_umap_hdbscan(module_df, metadata_df, nn, md, mc, output_path, "plots")


from PIL import Image
import math

def tile_images(image_dir, output_prefix, padding=10, bg_color=(255,255,255,255)):
    # Find PNG files
    png_files = []
    for f in os.listdir(image_dir):
        if f.endswith(".png"):
            if f.startswith("umap_nn") or f.startswith("hdbscan_nn"):
                png_files.append(os.path.join(image_dir, f))
    
    print(f"Found {len(png_files)} PNG files")
    
    # Extract parameter values
    nn_values = set()
    md_values = set()
    for filename in png_files:
        base = os.path.basename(filename)
        parts = base.replace(".png", "").split("_")
        try:
            nn = int(parts[1].replace("nn", ""))
            md = float(parts[2].replace("md", ""))
            nn_values.add(nn)
            md_values.add(md)
        except (IndexError, ValueError):
            continue

    nn_values = sorted(list(nn_values)) 
    md_values = sorted(list(md_values))
    mc_values = sorted(${params.mc})

    # Process each md value separately
    for md in md_values:
        md_int = int(md)
        images_by_row = {}  # Dictionary to store images for each row (nn value)
        
        # First, collect UMAP images for each nn value
        for nn in nn_values:
            umap_file = f"umap_nn{nn}_md{md_int}.png"
            umap_path = os.path.join(image_dir, umap_file)
            if os.path.exists(umap_path):
                print(f"Loading UMAP: {umap_file}")
                img = Image.open(umap_path)
                #if img.size[0] > 1200:
                 #   ratio = 1200 / img.size[0]
                  #  new_size = (1200, int(img.size[1] * ratio))
                   # img = img.resize(new_size, Image.Resampling.LANCZOS)
                images_by_row[nn] = [img] 
        
        # Then add HDBSCAN images for each nn value
        for nn in nn_values:
            if nn in images_by_row:  # Only process if we have a UMAP image for this nn
                for mc in mc_values:
                    hdbscan_file = f"hdbscan_nn{nn}_md{md_int}_minclust{mc}.png"
                    hdbscan_path = os.path.join(image_dir, hdbscan_file)
                    if os.path.exists(hdbscan_path):
                        print(f"Loading HDBSCAN: {hdbscan_file}")
                        img = Image.open(hdbscan_path)
                        if img.size[0] > 1200:
                            ratio = 1200 / img.size[0]
                            new_size = (1200, int(img.size[1] * ratio))
                            img = img.resize(new_size, Image.Resampling.LANCZOS)
                        images_by_row[nn].append(img)

        # Calculate grid dimensions
        rows = len(nn_values)
        cols = 1 + len(mc_values)  # UMAP + HDBSCAN plots
        
        max_width = max(img.size[0] for row_images in images_by_row.values() for img in row_images)
        max_height = max(img.size[1] for row_images in images_by_row.values() for img in row_images)

        grid_width = (max_width * cols) + (padding * (cols - 1))
        grid_height = (max_height * rows) + (padding * (rows - 1))
        
        label_height = 300
        label_width = 300  
        axis_title_height = 100

        # Create canvas
        total_width = grid_width + label_width + padding * 2
        total_height = grid_height + label_height + axis_title_height * 2
        canvas = Image.new('RGBA', (total_width, total_height), bg_color)

        draw = ImageDraw.Draw(canvas)
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 40)

        # Add main title
        title = f"UMAP and HDBSCAN Clustering (md={md/10:.1f})"
        bbox = draw.textbbox((0, 0), title, font=font)
        text_width = bbox[2] - bbox[0]
        draw.text((grid_width//2, axis_title_height//2), title, 
                fill='black', font=font, anchor='mm')

        # Add X-axis title
        x_title = "Minimum Cluster Size"
        x_title_img = Image.new('RGBA', (grid_width, axis_title_height), bg_color)
        x_draw = ImageDraw.Draw(x_title_img)
        x_draw.text((grid_width//2, axis_title_height//2), x_title, 
                fill='black', font=font, anchor='mm')
        canvas.paste(x_title_img, (label_width, axis_title_height))

       
        y_title = "Nearest Neighbors Value (nn)"
        y_title_img = Image.new('RGBA', (grid_height, axis_title_height), bg_color)
        y_draw = ImageDraw.Draw(y_title_img)
        y_draw.text((grid_height//2, axis_title_height//2), y_title, 
                fill='black', font=font, anchor='mm')
        y_title_img = y_title_img.rotate(90, expand=True)
        canvas.paste(y_title_img, (padding * 2, label_height + (grid_height // 2) - (y_title_img.height // 2)))

        # Add column labels
        draw.text((label_width + max_width//2, axis_title_height * 2), "UMAP",  # Positioned closer to title
                fill='black', font=font, anchor='mm')
        for col, mc in enumerate(mc_values, 1):
            x = label_width + col * (max_width + padding) + max_width//2
            draw.text((x, axis_title_height * 2), f"MC={mc}",  # Positioned closer to title
                    fill='black', font=font, anchor='mm')

        # Add row labels
        for row, nn in enumerate(nn_values):
            # Create new image for rotated text
            label_text = f"nn={nn}"
            label_img = Image.new('RGBA', (max_height, label_width//2), bg_color)
            label_draw = ImageDraw.Draw(label_img)
            
            # Draw text horizontally (it will be rotated later)
            label_draw.text((max_height//2, label_width//4), label_text,
                            fill='black', font=font, anchor='mm')
            
            # Rotate the image 90 degrees
            label_img = label_img.rotate(90, expand=True)
            
            # Calculate position with reduced spacing
            x = label_width//3  # Reduced from label_width//2
            y = label_height + row * (max_height + padding) + max_height//2 - label_img.height//2
            
            # Paste the rotated label
            canvas.paste(label_img, (x, y), label_img)

        # Paste images
        for row_idx, nn in enumerate(nn_values):
            if nn in images_by_row:
                for col_idx, img in enumerate(images_by_row[nn]):
                    x = label_width + col_idx * (max_width + padding)
                    y = label_height + row_idx * (max_height + padding)
                    canvas.paste(img, (x, y))
        
        # Save the tiled image for this md value
        output_file = f"{output_prefix}_md{md_int}.png"
        canvas.save(output_file, optimize=True, quality=95, dpi=(600, 600))
        print(f"Saved tiled image for md={md:.1f}: {output_file}")

# Call the function
tile_images("plots", "plots/tiled_image")
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
from PIL import Image, ImageDraw, ImageFont

# Create output directory
os.makedirs("plots", exist_ok=True)

def plot_umap(module_df, metadata_df, nn, md, output_file):
    # Load pre-generated UMAP coordinates and IDs
    md_int = int(md * 10)
    coord_file = os.path.join("${coordinates_dir}", f"coords/coordinates_nn{nn}_md{md_int}.tsv")
    conns_file = os.path.join("${coordinates_dir}", "connections.tsv")

    try:
        # Load coordinates and their corresponding IDs
        coord_df = pd.read_csv(coord_file, sep='\t')
        embedding_ids = coord_df['embedding_id'].tolist()
        embedding_2d = coord_df[['x', 'y']].values
        
        print(f"Loaded coordinates from {coord_file}")
        print(f"Shape: {embedding_2d.shape}, IDs: {len(embedding_ids)}")
        
        # Prepare metadata mappings
        module_df["normalized_orf_id"] = module_df["orf_id"].str.replace(".", "-", regex=False)
        # id_to_genofeature = dict(zip(module_df["normalized_orf_id"], module_df["genofeature"]))
        # genofeature_to_color = dict(zip(metadata_df["genofeature"], metadata_df["color"]))
        # genofeature_to_marker = dict(zip(metadata_df["genofeature"], metadata_df["marker"]))
        # genofeature_to_display = dict(zip(metadata_df["genofeature"], metadata_df["display_name"]))

        # Assign colors and markers
        colors = []
        markers = []
        for eid in embedding_ids:
            genofeature = module_df.loc[module_df["normalized_orf_id"] == eid, "genofeature"].iloc[0] if eid in module_df["normalized_orf_id"].values else None
            if genofeature is not None and genofeature in metadata_df["genofeature"].values:
                color = metadata_df.loc[metadata_df["genofeature"] == genofeature, "color"].iloc[0]
                marker = metadata_df.loc[metadata_df["genofeature"] == genofeature, "marker"].iloc[0]
            else:
                color = "#808080"
                marker = "."
            colors.append(color)
            markers.append(marker if pd.notna(marker) else ".")


        if os.path.exists(conns_file):
            connections_df = pd.read_csv(conns_file, sep='\t')
            print(f"Loaded {len(connections_df)} connections")
            
            # Create index mapping for connections
            id_to_idx = {eid: idx for idx, eid in enumerate(embedding_ids)}
            
            # Convert connections to index pairs
            connections = []
            for _, row in connections_df.iterrows():
                if row['id1'] in id_to_idx and row['id2'] in id_to_idx:
                    connections.append([id_to_idx[row['id1']], id_to_idx[row['id2']]])
            connections = np.array(connections)
        else:
            connections = []
            print("No connections file found")

        # Create plot
        fig, ax = plt.subplots(figsize=(10, 7))

        # Draw connections first (background)
        if len(connections) > 0:
            for conn in connections:
                x_coords = [embedding_2d[conn[0]][0], embedding_2d[conn[1]][0]]
                y_coords = [embedding_2d[conn[0]][1], embedding_2d[conn[1]][1]]
                ax.plot(x_coords, y_coords, color='#CCCCCC', alpha=0.5, linewidth=0.2, zorder=1)


        # Plot points
        for coord, color, marker in zip(embedding_2d, colors, markers):
            ax.scatter(coord[0], coord[1], c=color, marker=marker, s=10, alpha=0.7)

        # Configure plot
        ax.set_xticks([])
        ax.set_yticks([])
        plt.grid(False)

        # Create legend handles
        unique_genofeatures = sorted(metadata_df["genofeature"].unique())
        legend_handles = []
        
        for genofeature in unique_genofeatures:
            row = metadata_df.loc[metadata_df["genofeature"] == genofeature].iloc[0]
            handle = mlines.Line2D([], [], 
                                 color=row["color"],
                                 marker=row["marker"],
                                 linestyle='None',
                                 markersize=8,
                                 label=row["display_name"])
            legend_handles.append(handle)

        # Save plot
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

# Generate plots for specified parameters
legend_handles = None
for nn in ${params.nn}:
    for md in ${params.md}:
        output_file = f"plots/umap_nn{nn}_md{int(md*10)}.png"
        handles = plot_umap(module_df, metadata_df, nn, md, output_file)
        if legend_handles is None:
            legend_handles = handles

from PIL import Image

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
                nn_part = parts[1]  
                nn = int(nn_part.replace("nn", ""))  # properly removes "nn" prefix
                md_part = parts[2]  
                md = float(md_part.replace("md", "")) / 10
                nn_values.add(nn)
                md_values.add(md)
            except (IndexError, ValueError) as e:
                print(f"Skipping filename: {filename}")
                continue

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
            print(f"Error parsing {filename}: {e}")
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

process GENERATE_COORDINATES {
    publishDir "${params.outdir}",
        mode: 'copy'
        
    label 'gpu'  
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    memory { 50.GB }
    
    input:
        path(embeddings)   
        
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
    import pandas as pd
    from pathlib import Path

    # Create output directories
    os.makedirs("first_coordinates/coords", exist_ok=True)
   

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

        groups = {}
        for i, eid in enumerate(embedding_ids):
            contig_id = '_'.join(eid.split('_')[:-3])  # Parse embedding ID
            if contig_id not in groups:
                groups[contig_id] = []
            groups[contig_id].append(i)
        
        # Generate UMAP for each parameter combination
        for nn in ${params.nn}:
            for md in ${params.md}:
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
                coord_file = f"first_coordinates/coords/coordinates_nn{nn}_md{md_int}.tsv"

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
        connections_file = f"first_coordinates/connections.tsv"
        connections_df.to_csv(connections_file, sep='\t', index=False)

    """
}




workflow {
    // Create input channels
    ch_filtered_tsv = Channel.fromPath(params.second_run)
    ch_metadata = Channel.fromPath(params.genofeature_metadata)

    ch_embedding_datasets = Channel.value(params.embedding_datasets)

    // ch_embeddings = EMBEDDINGS(
    //     ch_embedding_datasets
    // )

    // ch_coordinates = GENERATE_COORDINATES(
    //     // ch_embeddings
    //     "/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/full_ena_output/embeddings"
    // )

    // remember to fix the input from coordinates the same way you did to pasv output.
    ch_umap = UMAP_PROJECTION(
        "/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/figures_folder/coordinates",
        // ch_coordinates,
        ch_filtered_tsv,
        ch_metadata
    )

    ch_hbd = HDBSCAN(
        "/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/figures_folder/coordinates",
        // ch_coordinates,
        ch_filtered_tsv,
        ch_metadata,
        ch_umap.plots
        // "/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/work/0b/eab48116d13496090e556b588f6dfa/plots"
        
    )
}