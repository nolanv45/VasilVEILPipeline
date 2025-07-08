nextflow.enable.dsl=2

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
    def selected_list = params.selected_genofeatures.collect { "\'${it}\'" }.join(', ')
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    mkdir -p embeddings

    extract_script="/mnt/VEIL/tools/embeddings/models/extract.py"
    model_path="/mnt/VEIL/tools/embeddings/models/esm2_t36_3B_UR50D.pt"
    ls -R "${split_dir}"

    # Get version info
    python_version=\$(python -V | sed 's/Python //')
    pytorch_version=\$(python -c "import torch; print(torch.__version__)")
    model_version=\$(python -c "import torch; m=torch.load('\${model_path}', map_location='cpu'); print(m['model_data']['model_version'] if 'model_data' in m and 'model_version' in m['model_data'] else 'unknown')")

    # Process each protein directory
    for protein_dir in "${split_dir}"/*; do
        if [ -d "\$protein_dir" ]; then
            protein=\$(basename "\$protein_dir")
            echo "Processing protein directory: \$protein"
            mkdir -p "embeddings/\$protein"
            
            # Process genofeature directories, but only if they're in selected_list
            for genofeature_dir in "\$protein_dir"/*; do
                if [ -d "\$genofeature_dir" ]; then
                    genofeature=\$(basename "\$genofeature_dir")
                    # Check if genofeature is in selected_list (case-insensitive)
                    if echo "${selected_list}" | tr '[:upper:]' '[:lower:]' | grep -iq "\'\$(echo \$genofeature | tr '[:upper:]' '[:lower:]')\'"; then
                        echo "Processing genofeature directory: \$genofeature"
                        output_dir="embeddings/\$protein/\$genofeature"                 
                        mkdir -p "\$output_dir"
               
                        # Process FASTA files
                        for fasta in "\$genofeature_dir"/*.fasta; do
                            if [ -f "\$fasta" ]; then
                                echo "Processing FASTA: \$fasta"


                                # Extract base name of the FASTA file
                                fasta_base=\$(basename "\$fasta" .fasta)

                                # Check if embeddings already exist
                                if [ -d "${params.outdir}/embeddings/\$protein/\$genofeature/batch_0" ] && [ -n "\$(find "${params.outdir}/embeddings/\$protein/\$genofeature/batch_0" -name '*.pt' 2>/dev/null)" ]; then
                                    echo "Embeddings already exist in output directory for \$fasta_base, skipping..."
                                    # Copy existing embeddings to work directory to ensure output is complete
                                    cp -r "${params.outdir}/embeddings/\$protein/\$genofeature" "embeddings/\$protein/"
                                    continue
                                fi
                                echo "Generating embeddings for \$fasta_base..."
                                /home/nolanv/.conda/envs/esm-umap/bin/python \\
                                    "\$extract_script" \\
                                    "\$model_path" \\
                                    "\$fasta" \\
                                    "\$output_dir" \\
                                    --repr_layers 36 \\
                                    --include mean || exit 1
                            fi
                        done
                    else
                        echo "Skipping genofeature \$genofeature (not in selected list)"
                    fi
                fi
            done
        fi
    done
    """
}


process HDBSCAN {
    publishDir "${params.outdir}/hdbscan",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.startsWith("plots/tiled_image_md")) filename
            else if (filename.endsWith(".png")) "plots/${filename}"
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
def load_and_resize_image(filepath, max_width=1200):
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
        coord_file = os.path.join("${coordinates_dir}", f"coords/coordinates_nn{nn}_md{md_int}.npy")
        id_file = os.path.join("${coordinates_dir}", f"ids/embedding_ids_nn{nn}_md{md_int}.txt")
        
        # Generate HDBSCAN plots with different min_cluster_size values
        for mc in [10, 20, 30, 40]:
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

    nn_values = sorted(list(nn_values))  # [75, 100, 125]
    md_values = sorted(list(md_values))
    mc_values = sorted([10, 20, 30, 40])

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
                if img.size[0] > 1200:
                    ratio = 1200 / img.size[0]
                    new_size = (1200, int(img.size[1] * ratio))
                    img = img.resize(new_size, Image.Resampling.LANCZOS)
                images_by_row[nn] = [img]  # Start row with UMAP image
        
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
        title = f"UMAP and HDBSCAN Clustering (md={md:.1f})"
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

        # Add Y-axis title (Number of Neighbors)
        y_title = "Number of Neighbors (nn)"
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
        canvas.save(output_file, optimize=True, quality=85)
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

process GENERATE_COORDINATES {
    publishDir "${params.outdir}",
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
    os.makedirs("coordinates/coords", exist_ok=True)
    os.makedirs("coordinates/ids", exist_ok=True)

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
                coord_file = f"coordinates/coords/coordinates_nn{nn}_md{md_int}.npy"
                id_file = f"coordinates/ids/embedding_ids_nn{nn}_md{md_int}.txt"
                
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
    
    ch_split_dir = Channel.fromPath(params.split_dir)

    // ch_embeddings = EMBEDDINGS(
    //     ch_split_dir
    // )

    // ch_coordinates = GENERATE_COORDINATES(
    //     ch_embeddings,
    //     ch_filtered_tsv
    // )

    // remember to fix the input from coordinates the same way you did to pasv output.
    // ch_umap = UMAP_PROJECTION(
    //     "/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/output_test/coordinates",
    //     // ch_coordinates,
    //     ch_filtered_tsv,
    //     ch_metadata
    // )

    ch_hbd = HDBSCAN(
        "/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/output_test/coordinates",
        // ch_coordinates,
        ch_filtered_tsv,
        ch_metadata,
        // ch_umap.plots
        "/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/work/da/3a0637061590d10d1b61f71e51a609/plots"
        
    )
}