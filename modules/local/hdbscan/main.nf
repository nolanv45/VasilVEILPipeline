process HDBSCAN {
    publishDir "${params.outdir}/04_parameter_selection/hdbscan",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.startsWith("plots/tiled_image_md")) filename
            else if (filename.startsWith("clusters_csv/")) filename
            else null
        }

    label "process_medium"
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/umap/umap.yml" 
    // container "containers/umap/umap.sif"
    
    input:
        path embeddings_dirs
        path coordinates_dir  // Directory containing the pre-generated coordinates
        path filtered_tsv    // TSV file with metadata
        path metadata_file   // Metadata file with colors/markers
        path plots
        
    output:
        path "plots/*.png", emit: plots
        path "plots/tiled_image_md*", emit: tiled_image
        path "clusters_csv/*.csv", emit: clusters_csv
        path "versions.yml", emit: versions

    script:
    """
#!/usr/bin/env python3
import os
# Ensure matplotlib has a writable config dir in container/work environments
os.environ.setdefault('MPLCONFIGDIR', '/tmp/matplotlib')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import hdbscan
import random 
import torch
from pathlib import Path
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

def find_pt_files(base_dir):
    pt_files = []
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith('.pt'):
                pt_files.append(os.path.join(root, file))
    return pt_files

def load_raw_embeddings(embeddings_dirs_str):
    embeddings = []
    embedding_ids = []
    base_dirs = embeddings_dirs_str.split() if isinstance(embeddings_dirs_str, str) else embeddings_dirs_str
    pt_files = []
    for base_dir in base_dirs:
        pt_files.extend(find_pt_files(base_dir))
    print(f"Found {len(pt_files)} .pt files")
    for file_path in pt_files:
        try:
            data = torch.load(file_path, map_location='cpu')
            emb = data["mean_representations"][36].numpy()
            eid = data.get("label")
            if emb is not None and eid is not None:
                embeddings.append(emb)
                embedding_ids.append(eid)
        except Exception as e:
            print(f"Warning: could not load {file_path}: {e}")
    return np.stack(embeddings), embedding_ids

def plot_umap_hdbscan(module_df, metadata_df, nn, md, mc, output_path, iteration_output, raw_embeddings, raw_embedding_ids):
    try:
        md_int = int(md * 10)

        # Load 2D coords for plotting only
        coord_file = os.path.join("${coordinates_dir}", f"coordinates_nn{nn}_md{md_int}.tsv")
        coord_df = pd.read_csv(coord_file, sep='\t')
        embedding_ids_2d = coord_df['embedding_id'].tolist()
        embedding_2d = coord_df[['x', 'y']].values

        # Align raw embeddings to the order of the 2D coordinate IDs
        id_to_raw = {eid: emb for eid, emb in zip(raw_embedding_ids, raw_embeddings)}
        aligned_embeddings = np.stack([id_to_raw[eid] for eid in embedding_ids_2d if eid in id_to_raw])
        aligned_ids = [eid for eid in embedding_ids_2d if eid in id_to_raw]
        aligned_2d = np.stack([embedding_2d[i] for i, eid in enumerate(embedding_ids_2d) if eid in id_to_raw])

        print(f"Clustering on {aligned_embeddings.shape[1]}-dimensional embeddings, plotting on 2D coords")

        # Cluster on HIGH-DIMENSIONAL embeddings
        labels = hdbscan.HDBSCAN(min_samples=2, min_cluster_size=mc).fit_predict(aligned_embeddings)
        clustered = (labels >= 0)

        module_df["normalized_orf_id"] = module_df["orf_id"].str.replace(".", "-", regex=False)

        cluster_df = pd.DataFrame({
            'embedding_id': aligned_ids,
            'cluster_label': labels,
            'genofeature': [module_df.loc[module_df["normalized_orf_id"] == eid, "genofeature"].iloc[0]
                if eid in module_df["normalized_orf_id"].values else "unknown"
                for eid in aligned_ids]
        })

        cluster_output_path = os.path.join(iteration_output,
                                           f"hdbscan_nn{nn}_md{md_int}_minclust{mc}_clusters.csv")
        cluster_df.to_csv(cluster_output_path, index=False)

        # Plot using 2D coords but colour by high-dim cluster labels
        fig, ax = plt.subplots(figsize=(10, 10))
        scatter = ax.scatter(aligned_2d[clustered, 0], aligned_2d[clustered, 1],
                             c=labels[clustered], s=10, alpha=0.7, cmap="Spectral")
        ax.scatter(aligned_2d[~clustered, 0], aligned_2d[~clustered, 1],
                   color="gray", s=10, alpha=0.5, label="not clustered")

        unique_labels = np.unique(labels[labels >= 0])
        for cluster in unique_labels:
            cluster_points = aligned_2d[labels == cluster]
            cx, cy = cluster_points[:, 0].mean(), cluster_points[:, 1].mean()
            plt.text(cx, cy, str(cluster), fontsize=10, ha='center', va='center')

        ax.set_xticks([])
        ax.set_yticks([])
        plt.grid(False)
        plt.title(f"UMAP (nn={nn}, md={md}) with HDBSCAN (minclust={mc})")
        plt.colorbar(scatter, label="Cluster Labels")
        plt.tight_layout()
        plt.savefig(output_path, dpi=600, bbox_inches="tight")
        plt.close()
        return True

    except Exception as e:
        print(f"Error processing plot: {e}")
        return False


# Main execution
print("Loading metadata...")
module_df = pd.read_csv("${filtered_tsv}", sep='\\t')
metadata_df = pd.read_csv("${metadata_file}", sep='\\t')


# Create required directories
os.makedirs("plots", exist_ok=True)
os.makedirs("clusters_csv", exist_ok=True)

# Copy input UMAP plots to working directory
os.system("cp ${plots} plots/")
# Load raw embeddings once before the loop
raw_embeddings, raw_embedding_ids = load_raw_embeddings("${embeddings_dirs}")

for umap_file in os.listdir("plots"):
    if umap_file.startswith("umap_nn"):
        parts = umap_file.replace(".png", "").split("_")
        nn = int(parts[1].replace("nn", ""))
        md = float(parts[2].replace("md", "")) / 10
        md_int = int(md * 10)
        for mc in ${params.mc}:
            output_path = f"plots/hdbscan_nn{nn}_md{md_int}_minclust{mc}.png"
            plot_umap_hdbscan(module_df, metadata_df, nn, md, mc, output_path, "clusters_csv",
                              raw_embeddings, raw_embedding_ids)

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
        
        # First, collect UMAP images for each nn value (these are the reference sizes)
        for nn in nn_values:
            umap_file = f"umap_nn{nn}_md{md_int}.png"
            umap_path = os.path.join(image_dir, umap_file)
            if os.path.exists(umap_path):
                print(f"Loading UMAP: {umap_file}")
                img = Image.open(umap_path)
                images_by_row[nn] = [img]

        # Then add HDBSCAN images for each nn value and resize them to match the UMAP
        for nn in list(images_by_row.keys()):
            umap_w, umap_h = images_by_row[nn][0].size
            for mc in mc_values:
                hdbscan_file = f"hdbscan_nn{nn}_md{md_int}_minclust{mc}.png"
                hdbscan_path = os.path.join(image_dir, hdbscan_file)
                if os.path.exists(hdbscan_path):
                    print(f"Loading HDBSCAN: {hdbscan_file}")
                    img = Image.open(hdbscan_path)
                    # Resize HDBSCAN images to match UMAP size for visual parity
                    if img.size != (umap_w, umap_h):
                        img = img.convert('RGB')
                        img = img.resize((umap_w, umap_h), Image.Resampling.LANCZOS)
                    images_by_row[nn].append(img)

        # Calculate grid dimensions using only rows that have images
        ordered_rows = sorted(images_by_row.keys())
        rows = len(ordered_rows)
        cols = 1 + len(mc_values)  # UMAP + HDBSCAN plots

        # Since HDBSCAN panels were resized to match each row's UMAP,
        # use the UMAP (first image) sizes to compute max cell size.
        max_width = max(row_images[0].size[0] for row_images in images_by_row.values())
        max_height = max(row_images[0].size[1] for row_images in images_by_row.values())

        grid_width = (max_width * cols) + (padding * (cols - 1))
        grid_height = (max_height * rows) + (padding * (rows - 1))
        
        label_height = 350
        label_width = 350
        axis_title_height = 100

        # Create canvas
        total_width = grid_width + label_width + padding * 2
        total_height = grid_height + label_height + axis_title_height * 2
        canvas = Image.new('RGBA', (total_width, total_height), bg_color)

        # Create a drawing object
        draw = ImageDraw.Draw(canvas)
        # Use a smaller fixed font so top labels don't crowd
        font = ImageFont.truetype("${projectDir}/fonts/DejaVuSans-Bold.ttf", 100)

        # Add main title (placed within the top label band, centered over the grid)
        title = f"UMAP and HDBSCAN Clustering (md={md/10:.1f})"
        title_x = label_width + grid_width // 2
        title_y = int(label_height * 0.18)
        draw.text((title_x, title_y), title, fill='black', font=font, anchor='mm')

        # Add X-axis title within the top label area (above column labels)
        x_title = "Minimum Cluster Size"
        x_title_img = Image.new('RGBA', (grid_width, axis_title_height), bg_color)
        x_draw = ImageDraw.Draw(x_title_img)
        x_draw.text((grid_width//2, axis_title_height//2), x_title, fill='black', font=font, anchor='mm')
        x_title_y = int(label_height * 0.32)
        canvas.paste(x_title_img, (label_width, x_title_y))

        # Add column labels below the title (moved further down to avoid overlap)
        col_label_y = int(label_height * 0.78)
        draw.text((label_width + max_width//2, col_label_y), "UMAP",
                  fill='black', font=font, anchor='mm')
        for col, mc in enumerate(mc_values, 1):
            x = label_width + col * (max_width + padding) + max_width//2
            draw.text((x, col_label_y), f"MC={mc}", fill='black', font=font, anchor='mm')

        # Add Y-axis title (rotated) at left of the image grid
        y_title = "Nearest Neighbors Value (nn)"
        y_title_img = Image.new('RGBA', (grid_height, axis_title_height), bg_color)
        y_draw = ImageDraw.Draw(y_title_img)
        y_draw.text((grid_height//2, axis_title_height//2), y_title,
                fill='black', font=font, anchor='mm')
        y_title_img = y_title_img.rotate(90, expand=True)
        canvas.paste(y_title_img, (padding * 2, label_height + (grid_height // 2) - (y_title_img.height // 2)))

        # (X-axis title already placed in the top label band)

        # Add row labels (only for rows we actually have)
        for row, nn in enumerate(ordered_rows):
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

        # Paste images (create blank placeholders for missing panels)
        placeholder = Image.new('RGBA', (max_width, max_height), bg_color)
        for row_idx, nn in enumerate(ordered_rows):
            row_imgs = images_by_row.get(nn, [])
            for col_idx in range(cols):
                if col_idx < len(row_imgs):
                    img = row_imgs[col_idx]
                else:
                    img = placeholder
                x = label_width + col_idx * (max_width + padding)
                y = label_height + row_idx * (max_height + padding)
                canvas.paste(img, (x, y))
        
        # Save the tiled image for this md value
        output_file = f"{output_prefix}_md{md_int}.png"
        canvas.save(output_file, optimize=True, quality=95, dpi=(600, 600))
        print(f"Saved tiled image for md={md:.1f}: {output_file}")

# Call the function
tile_images("plots", "plots/tiled_image")

import platform
from importlib import import_module

versions = [
    f'"${task.process}":',
    f'    python: {platform.python_version()}'
]

for package_name in ['numpy', 'pandas', 'hdbscan', 'torch']:
    try:
        module = import_module(package_name)
        version = getattr(module, '__version__', None)
        if version:
            versions.append(f'    {package_name}: {version}')
    except Exception:
        pass

with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write('\\n'.join(versions) + '\\n')
    """
}