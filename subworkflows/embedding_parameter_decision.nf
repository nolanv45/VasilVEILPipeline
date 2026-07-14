nextflow.enable.dsl=2

process EMBEDDING_PLAN {

    label "process_single"

    input:
    path(combined_tsv)

    output:
    path("**/*.fasta"), emit: planned_fastas, optional: true
    path("versions.yml"), emit: versions

    script:
    def datasetsJson = groovy.json.JsonOutput.toJson(params.datasets)
    """
    #!/usr/bin/env bash
    set -euo pipefail

    python3 - <<'PY'
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

combined_tsv = "${combined_tsv}"
datasets = json.loads('''${datasetsJson}''')

groups = defaultdict(set)

def normalize_id(value: str) -> str:
    token = (value or '').strip().split()[0]
    token = token.replace('.', '-')
    token = re.sub(r'[^A-Za-z0-9_-]', '', token)
    return token

with open(combined_tsv, 'r', encoding='utf-8', newline='') as combined_handle:
    reader = csv.DictReader(combined_handle, delimiter='\\t')
    fieldnames = set(reader.fieldnames or [])
    required_cols = {'dataset', 'protein', 'genofeature', 'orf_id'}
    missing_cols = required_cols - fieldnames
    if missing_cols:
        raise KeyError(f"Missing required columns in combined TSV: {sorted(missing_cols)}")

    for row in reader:
        dataset = (row.get('dataset') or '').strip()
        protein = (row.get('protein') or '').strip()
        genofeature = (row.get('genofeature') or '').strip()
        orf_id = normalize_id(row.get('orf_id') or '')
        if dataset and protein and genofeature and orf_id:
            groups[(dataset, protein, genofeature)].add(orf_id)

for (dataset, protein, genofeature), expected_ids in groups.items():

    fasta_path = Path(datasets.get(dataset, ''))
    if not str(fasta_path):
        raise KeyError(f"No fasta configured for dataset {dataset}")
    if not fasta_path.exists():
        raise FileNotFoundError(f"Configured fasta does not exist for {dataset}: {fasta_path}")

    output_root = Path("${params.outdir}") / "embeddings" / dataset / protein / genofeature

    existing_ids = {
        normalize_id(pt.stem)
        for pt in output_root.rglob('*.pt')
    } if output_root.exists() else set()

    # Fallback guard for slight naming differences between TSV IDs and pt stems.
    if len(existing_ids) >= len(expected_ids):
        continue

    missing_ids = sorted(expected_ids.difference(existing_ids))

    if not missing_ids:
        continue

    missing_set = set(missing_ids)
    out_fasta = Path(dataset) / protein / genofeature / f"{dataset}_{protein}_{genofeature}.fasta"
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    with fasta_path.open('r', encoding='utf-8') as fasta_handle, out_fasta.open('w', encoding='utf-8') as out_handle:
        header = None
        sequence_lines = []
        written = {'count': 0}

        def flush_record(record_header, record_lines):
            if record_header is None:
                return
            record_id = normalize_id(record_header)
            if record_id in missing_set:
                out_handle.write(f">{record_header}\\n")
                for seq_line in record_lines:
                    out_handle.write(seq_line)
                    if not seq_line.endswith('\\n'):
                        out_handle.write('\\n')
                written['count'] += 1

        for raw_line in fasta_handle:
            if raw_line.startswith('>'):
                flush_record(header, sequence_lines)
                header = raw_line[1:].strip()
                sequence_lines = []
            else:
                sequence_lines.append(raw_line)

        flush_record(header, sequence_lines)

    # Do not emit empty FASTA files; that would trigger unnecessary EMBEDDINGS tasks.
    if out_fasta.stat().st_size == 0 or written['count'] == 0:
        out_fasta.unlink(missing_ok=True)
PY

    python3 - <<'PY'
import platform

process_name = "${task.process}"
python_version = platform.python_version()

try:
    import pandas as pd
    pandas_version = pd.__version__
except Exception:
    pandas_version = 'not_installed'

with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"{process_name}":\\n')
    handle.write(f'    python: {python_version}\\n')
    handle.write(f'    pandas: {pandas_version}\\n')
PY
    """
}

process EMBEDDINGS {
    publishDir "${params.outdir}",
        mode: 'copy'
        
    label 'process_gpu' 
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/umap/umap.yml" 
    // container "containers/umap/umap.sif"
    
    input:
        tuple val(dataset), val(protein), val(genofeature), path(fasta)
        path model_folder
        
    output:
        path "embeddings/${dataset}/${protein}/${genofeature}", emit: embeddings_dirs
        path "versions.yml", emit: versions
        

    script:
    """
    set -euo pipefail

    out_dir="embeddings/${dataset}/${protein}/${genofeature}"
    mkdir -p "\$out_dir"

    python3 ${model_folder}/extract.py \
        ${model_folder}/esm2_t36_3B_UR50D.pt \
        "${fasta}" \
        "\$out_dir" \
        --repr_layers 36 \
        --include mean

    for batch_dir in "\$out_dir"/batch_*; do
        [ -d "\$batch_dir" ] || continue
        mv "\$batch_dir"/*.pt "\$out_dir"/
        rmdir "\$batch_dir" || true
    done

    python3 - <<'PY'
import platform
from importlib import import_module

process_name = "${task.process}"
versions = [
    f'"{process_name}":',
    f'    python: {platform.python_version()}'
]

for package_name in ['numpy', 'pandas', 'torch', 'umap']:
    try:
        module = import_module(package_name)
        version = getattr(module, '__version__', None)
        if version:
            versions.append(f'    {package_name}: {version}')
    except Exception:
        pass

with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write('\\n'.join(versions) + '\\n')
PY

    """
}

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



process UMAP_PROJECTION {
    publishDir "${params.outdir}/04_parameter_selection/umap",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith("tiled_image.png")) {
                return "umap_parameter_tiled.png"
            }
            else null
        }
    label "process_medium"
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/umap/umap.yml" 
    // container "containers/umap/umap.sif"
    
    input:
        path coordinates_dir  // Directory containing the pre-generated coordinates
        path filtered_tsv    // TSV file with metadata
        path metadata_file   // Metadata file with colors/markers
        
    output:
        path "plots/tiled_image.png", emit: tiled_image
        path "plots/*.png", emit: plots
        path "versions.yml", emit: versions

    script:
    """
#!/usr/bin/env python3
import os
# Ensure matplotlib has a writable config dir in container/work environments
os.environ.setdefault('MPLCONFIGDIR', '/tmp/matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set non-interactive backend
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from PIL import Image, ImageDraw, ImageFont, __version__ as PIL_VERSION

# Create output directory
os.makedirs("plots", exist_ok=True)

def plot_umap(module_df, metadata_df, nn, md, output_file):
    # Load pre-generated UMAP coordinates and IDs
    md_int = int(md * 10)
    coord_file = os.path.join("${coordinates_dir}", f"coordinates_nn{nn}_md{md_int}.tsv")
    conns_file = os.path.join("${coordinates_dir}", "connections.tsv")

    try:
        # Load coordinates and their corresponding IDs
        coord_df = pd.read_csv(coord_file, sep='\\t')
        embedding_ids = coord_df['embedding_id'].tolist()
        embedding_2d = coord_df[['x', 'y']].values
        
        print(f"Loaded coordinates from {coord_file}")
        print(f"Shape: {embedding_2d.shape}, IDs: {len(embedding_ids)}")
        
        # Prepare metadata mappings
        module_df["normalized_orf_id"] = module_df["orf_id"].str.replace(".", "-", regex=False)

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
            connections_df = pd.read_csv(conns_file, sep='\\t')
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
        fig, ax = plt.subplots(figsize=(10, 10))

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
    # try:
    font = ImageFont.truetype("${projectDir}/fonts/DejaVuSans-Bold.ttf", 120)
    #except Exception:
     #   font = ImageFont.load_default()


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

import platform

process_name = "${task.process}"
versions = [
    f'"{process_name}":',
    f'    python: {platform.python_version()}',
    f'    numpy: {np.__version__}',
    f'    pandas: {pd.__version__}',
    f'    matplotlib: {matplotlib.__version__}',
    f'    PIL: {PIL_VERSION}'
]

with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write('\\n'.join(versions) + '\\n')
    """
}

process GENERATE_COORDINATES {
    publishDir "${params.outdir}/04_parameter_selection/coordinates",
        mode: 'copy'
        
    label "process_high_memory"
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/umap/umap.yml" 
    // container "containers/umap/umap.sif"

    input:
        path embeddings_dirs
        
    output:
        path "coords", emit: coordinates_files
        path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Ensure a writable Numba cache directory in the work folder
    mkdir -p \$PWD/.numba_cache
    chmod 1777 \$PWD/.numba_cache || true
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache

    # Run the Python script using the conda environment's python binary
#!/usr/bin/env python3
import os
import torch
import numpy as np
import umap
import random
import pandas as pd
import json
from pathlib import Path

# Create output directories
os.makedirs("coords", exist_ok=True)

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
        for file in files:
            if file.endswith('.pt'):
                pt_files.append(os.path.join(root, file))
    return pt_files

# Load embeddings
embeddings = []
embedding_ids = []

base_dirs = "${embeddings_dirs}".split()

# Find all .pt files
pt_files = []
for base_dir in base_dirs:
    if os.path.isdir(base_dir):
        pt_files.extend(find_pt_files(base_dir))
print(f"Found {len(pt_files)} .pt files")

# Process each .pt file
for file_path in pt_files:
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
            coord_file = f"coords/coordinates_nn{nn}_md{md_int}.tsv"

            coord_df = pd.DataFrame({
                'embedding_id': embedding_ids,
                'x': coordinates[:, 0],
                'y': coordinates[:, 1]
            })
            coord_df.to_csv(coord_file, sep='\\t', index=False)
            
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
    connections_file = f"coords/connections.tsv"
    connections_df.to_csv(connections_file, sep='\\t', index=False)

print("Writing versions.yml")
import platform
from importlib import import_module

versions = [
    f'"${task.process}":',
    f'    python: {platform.python_version()}'
]

for package_name in ['numpy', 'pandas', 'torch', 'umap']:
    try:
        module = import_module(package_name)
        version = getattr(module, '__version__', None)
        if version:
            versions.append(f'    {package_name}: {version}')
    except Exception:
        pass

with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write('\\n'.join(versions) + '\\n')

PY
    """
}


workflow EMBEDDING_PARAMETER_DECISION {
    take:
        ch_combined_tsv
        
    main:
    ch_filtered_tsv = ch_combined_tsv
    ch_metadata = Channel.fromPath(params.genofeature_metadata)

    EMBEDDING_PLAN(ch_combined_tsv)

    ch_embeddings = EMBEDDINGS(
        EMBEDDING_PLAN.out.planned_fastas.flatten().map { fasta ->
            def genofeature = fasta.parent.name
            def protein = fasta.parent.parent.name
            def dataset = fasta.parent.parent.parent.name
            tuple(dataset, protein, genofeature, fasta)
        },
        "${baseDir}/tools/"
    )

    ch_existing_embedding_dirs = Channel
        .fromPath("${params.outdir}/embeddings/*", type: 'dir')

    ch_embedding_dirs = ch_existing_embedding_dirs
        .mix(ch_embeddings.embeddings_dirs)
        .unique()

    GENERATE_COORDINATES(
        ch_embedding_dirs.collect()
    )

    ch_umap = UMAP_PROJECTION(
        GENERATE_COORDINATES.out.coordinates_files,
        ch_filtered_tsv,
        ch_metadata
    )

    ch_hbd = HDBSCAN(
        ch_embedding_dirs.collect(),
        GENERATE_COORDINATES.out.coordinates_files,
        ch_filtered_tsv,
        ch_metadata,
        ch_umap.plots       
    )

    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_umap.tiled_image)
    ch_multiqc_files = ch_multiqc_files.mix(ch_hbd.tiled_image)
    ch_multiqc_files = ch_multiqc_files.mix(ch_hbd.plots)

    ch_versions = channel.empty()
    ch_versions = ch_versions.mix(EMBEDDING_PLAN.out.versions)
    ch_versions = ch_versions.mix(EMBEDDINGS.out.versions)
    ch_versions = ch_versions.mix(HDBSCAN.out.versions)
    ch_versions = ch_versions.mix(UMAP_PROJECTION.out.versions)
    ch_versions = ch_versions.mix(GENERATE_COORDINATES.out.versions)

    emit:
        ch_combined_tsv = ch_combined_tsv
        versions = ch_versions
        multiqc_files = ch_multiqc_files
}