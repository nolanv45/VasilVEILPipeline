process TILE_UMAPS {
    publishDir "${params.outdir}/04_parameter_selection/umap", mode: 'copy'
    label "process_medium"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_pandas_matplotlib_pillow:bc19f6df598f92e8' :
        'community.wave.seqera.io/library/python_numpy_pandas_matplotlib_pillow:352bfae9d89bced3' }"
    
    input:
        path 'plots/*'
        path metadata_file
        
    output:
        path "tiled_image.png", emit: tiled_image
        path "versions.yml", emit: versions

     script:
    """
#!/usr/bin/env python3
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/matplotlib')
import io
import re
import platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, __version__ as PIL_VERSION

# ---------------------------------------------------------------------------
# Layout constants (all values in pixels, matched to dpi=600 plot outputs)
# ---------------------------------------------------------------------------
PADDING = 10
BG_COLOR = (255, 255, 255, 255)
LABEL_HEIGHT = 500       # vertical space reserved for column (md) labels
LABEL_WIDTH = 500        # horizontal space reserved for row (nn) labels
AXIS_TITLE_HEIGHT = 200  # space reserved for "Minimum Distance" / "Nearest Neighbors" titles
NN_LABEL_HEIGHT = 500    # working canvas size for a single rotated nn label, pre-rotation
NN_LABEL_WIDTH = 200
Y_TITLE_X = 50           # left margin before the rotated y-axis title
ROW_LABEL_GAP = 100      # gap between y-axis title and row labels
FONT_PATH = "${projectDir}/fonts/DejaVuSans-Bold.ttf"
FONT_SIZE = 120
LEGEND_FIGSIZE = (16, 36)
LEGEND_WIDTH_FRACTION = 0.25   # legend width as a fraction of grid height

FILENAME_RE = re.compile(r"umap_nn(\\d+)_md(\\d+p?\\d*)\\.png\$")

def parse_params(filename):
    match = FILENAME_RE.search(os.path.basename(filename))
    if not match:
        return None
    nn = int(match.group(1))
    md = float(match.group(2).replace('p', '.'))
    return nn, md


def build_grid_layout(image_dir):
    \"\"\"Discover, parse, and sort plot PNGs into a grid layout.

    Returns a dict with: png_files (sorted), images, nn_values, md_values,
    col_widths, row_heights, grid_width, grid_height.
    \"\"\"
    all_pngs = [
        os.path.join(image_dir, f)
        for f in os.listdir(image_dir)
        if f.lower().endswith(".png")
    ]

    parsed = {}
    for f in all_pngs:
        params = parse_params(f)
        if params is None:
            print(f"Skipping filename (does not match umap_nn<N>_md<N>.png): {f}")
            continue
        parsed[f] = params

    if not parsed:
        raise ValueError(f"No valid umap plot files found in {image_dir}")

    nn_values = sorted({nn for nn, _ in parsed.values()})
    md_values = sorted({md for _, md in parsed.values()})

    def sort_key(f):
        nn, md = parsed[f]
        return (nn_values.index(nn), md_values.index(md))

    png_files = sorted(parsed.keys(), key=sort_key)
    images = [Image.open(f) for f in png_files]
    sizes = [img.size for img in images]

    cols = len(md_values)
    rows = len(nn_values)

    col_widths = [0] * cols
    row_heights = [0] * rows
    for idx, (w, h) in enumerate(sizes):
        row, col = idx // cols, idx % cols
        col_widths[col] = max(col_widths[col], w)
        row_heights[row] = max(row_heights[row], h)

    grid_width = sum(col_widths) + PADDING * (cols - 1)
    grid_height = sum(row_heights) + PADDING * (rows - 1)

    return {
        "png_files": png_files,
        "images": images,
        "nn_values": nn_values,
        "md_values": md_values,
        "col_widths": col_widths,
        "row_heights": row_heights,
        "grid_width": grid_width,
        "grid_height": grid_height,
    }


def render_legend(legend_handles, height):
    \"\"\"Build the matplotlib legend and return it as a PIL Image sized to `height`.
    Uses an in-memory buffer instead of a temp file on disk.\"\"\"
    fig = plt.figure(figsize=LEGEND_FIGSIZE)
    ax = fig.add_subplot(111)
    ax.axis('off')

    if legend_handles:
        ax.legend(
            handles=legend_handles,
            title="Genofeature",
            loc="center left",
            bbox_to_anchor=(0, 0.5),
            frameon=False,
            fontsize=30,
            title_fontsize=36,
            borderaxespad=2,
            labelspacing=2.0,
            handlelength=4,
            handleheight=4,
            markerscale=4,
        )

    buf = io.BytesIO()
    fig.savefig(buf, bbox_inches='tight', dpi=150, facecolor='white',
                edgecolor='none', pad_inches=0.5)
    plt.close(fig)
    buf.seek(0)

    legend_img = Image.open(buf)
    width = int(height * LEGEND_WIDTH_FRACTION)
    legend_img = legend_img.resize((width, height), Image.Resampling.LANCZOS)
    return legend_img


def render_axis_titles(canvas, layout, font):
    \"\"\"Draw the 'Minimum Distance' and 'Nearest Neighbors' axis titles.\"\"\"
    grid_width = layout["grid_width"]
    grid_height = layout["grid_height"]

    x_title_img = Image.new('RGBA', (grid_width, AXIS_TITLE_HEIGHT), BG_COLOR)
    x_draw = ImageDraw.Draw(x_title_img)
    x_draw.text((grid_width // 2, AXIS_TITLE_HEIGHT // 2), "Minimum Distance (md)",
                fill='black', font=font, anchor='mm')
    canvas.paste(x_title_img, (LABEL_WIDTH, AXIS_TITLE_HEIGHT // 2))

    y_title_img = Image.new('RGBA', (grid_height, AXIS_TITLE_HEIGHT), BG_COLOR)
    y_draw = ImageDraw.Draw(y_title_img)
    y_draw.text((grid_height // 2, AXIS_TITLE_HEIGHT // 2), "Nearest Neighbors (nn)",
                fill='black', font=font, anchor='mm')
    y_title_img = y_title_img.rotate(90, expand=True)
    y_title_y = LABEL_HEIGHT + (grid_height // 2) - (y_title_img.height // 2)
    canvas.paste(y_title_img, (Y_TITLE_X, y_title_y), y_title_img)


def render_column_labels(draw, layout, font):
    \"\"\"Draw md-value labels along the top of the grid.\"\"\"
    col_widths = layout["col_widths"]
    for col, md in enumerate(layout["md_values"]):
        x = LABEL_WIDTH + sum(col_widths[:col]) + PADDING * col + col_widths[col] // 2
        y = AXIS_TITLE_HEIGHT + LABEL_HEIGHT // 2
        draw.text((x, y), f"md={md:.1f}", fill='black', font=font, anchor='mm')


def render_row_labels(canvas, layout, font):
    \"\"\"Draw nn-value labels along the left of the grid.\"\"\"
    row_heights = layout["row_heights"]
    x = Y_TITLE_X + AXIS_TITLE_HEIGHT + ROW_LABEL_GAP

    for row, nn in enumerate(layout["nn_values"]):
        label_img = Image.new('RGBA', (NN_LABEL_HEIGHT, NN_LABEL_WIDTH), BG_COLOR)
        label_draw = ImageDraw.Draw(label_img)
        label_draw.text((NN_LABEL_HEIGHT // 2, NN_LABEL_WIDTH // 2), f"nn={nn}",
                         fill='black', font=font, anchor='mm')
        label_img = label_img.rotate(90, expand=True)

        y = (LABEL_HEIGHT + sum(row_heights[:row]) + PADDING * row
             + row_heights[row] // 2 - label_img.height // 2)
        canvas.paste(label_img, (x, y), label_img)


def paste_grid_images(canvas, layout):
    \"\"\"Paste the individual umap plots into their grid positions.\"\"\"
    cols = len(layout["md_values"])
    col_widths = layout["col_widths"]
    row_heights = layout["row_heights"]

    for idx, img in enumerate(layout["images"]):
        row, col = idx // cols, idx % cols
        x = LABEL_WIDTH + sum(col_widths[:col]) + PADDING * col
        y = LABEL_HEIGHT + sum(row_heights[:row]) + PADDING * row
        canvas.paste(img, (x, y))


def tile_images(image_dir, output_file, legend_handles):
    layout = build_grid_layout(image_dir)
    grid_width = layout["grid_width"]
    grid_height = layout["grid_height"]

    total_width = grid_width + LABEL_WIDTH + PADDING * 2 + int(grid_height * 0.3)
    total_height = grid_height + LABEL_HEIGHT * 2 + AXIS_TITLE_HEIGHT * 2
    canvas = Image.new('RGBA', (total_width, total_height), BG_COLOR)
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype(FONT_PATH, FONT_SIZE)

    render_axis_titles(canvas, layout, font)
    render_column_labels(draw, layout, font)
    render_row_labels(canvas, layout, font)
    paste_grid_images(canvas, layout)

    legend_img = render_legend(legend_handles, grid_height)
    legend_x = LABEL_WIDTH + grid_width + PADDING
    legend_y = LABEL_HEIGHT
    canvas.paste(legend_img, (legend_x, legend_y))

    canvas.save(output_file)
    print(f"Tiled image with legend saved as {output_file}")


def build_legend_handles(metadata_df):
    \"\"\"Build legend handles directly from metadata (independent of any plot task).\"\"\"
    handles = []
    for genofeature in sorted(metadata_df["genofeature"].unique()):
        row = metadata_df.loc[metadata_df["genofeature"] == genofeature].iloc[0]
        handles.append(mlines.Line2D(
            [], [],
            color=row["color"],
            marker=row["marker"],
            linestyle='None',
            markersize=8,
            label=row["display_name"],
        ))
    return handles


def write_versions(output_path):
    process_name = "${task.process}"
    versions = [
        f'"{process_name}":',
        f'    python: {platform.python_version()}',
        f'    numpy: {np.__version__}',
        f'    pandas: {pd.__version__}',
        f'    matplotlib: {matplotlib.__version__}',
        f'    PIL: {PIL_VERSION}',
    ]
    with open(output_path, 'w', encoding='utf-8') as handle:
        handle.write('\\n'.join(versions) + '\\n')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
metadata_df = pd.read_csv("${metadata_file}", sep='\\t')
legend_handles = build_legend_handles(metadata_df)

tile_images("plots", "tiled_image.png", legend_handles)
write_versions("versions.yml")
    """
}