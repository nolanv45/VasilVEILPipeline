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
    conda "${moduleDir}/environment.yml"
    
    input:
        val coordinates_dir  // Space-separated list of coordinate directories
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
import glob
import re
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

def normalize_input_dirs(raw_value):
    raw_value = str(raw_value).strip()
    if not raw_value:
        return []
    cleaned = raw_value.strip('[]')
    return [item for item in re.split(r'[\\s,]+', cleaned) if item]

def plot_umap(module_df, metadata_df, coord_file, output_file):
    # Load pre-generated UMAP coordinates and IDs
    conns_file = os.path.join(os.path.dirname(coord_file), "connections.tsv")

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

# Generate plots for all coordinate files in all coordinate directories
coordinate_files = []
for base_dir in normalize_input_dirs("${coordinates_dir}"):
    coordinate_files.extend(sorted(glob.glob(os.path.join(base_dir, "coordinates_nn*_md*.tsv"))))

if not coordinate_files:
    raise ValueError("No coordinate files found in the supplied coordinate directories")

legend_handles = None
for coord_file in coordinate_files:
    base_name = os.path.basename(coord_file)
    parts = base_name.replace(".tsv", "").split("_")
    nn = int(parts[1].replace("nn", ""))
    md = float(parts[2].replace("md", "")) / 10
    output_file = f"plots/umap_nn{nn}_md{int(md*10)}.png"
    handles = plot_umap(module_df, metadata_df, coord_file, output_file)
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