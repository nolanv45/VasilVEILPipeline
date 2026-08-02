process TILE_HDBSCAN {
    label "process_medium"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pytorch_numpy_pandas_pruned:81ab5dec161369ec' :
        'community.wave.seqera.io/library/python_pytorch_numpy_pandas_pruned:6c2364c205cfbeb9' }"
    
    publishDir "${params.outdir}/04_parameter_selection/hdbscan",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.startsWith("plots/tiled_image_md")) filename
            else if (filename.startsWith("clusters_csv/")) filename
            else null
        }
    input:
        path 'plots/*'

    output:
        path "plots/tiled_image_md*.png", emit: tiled_image
        path "versions.yml", emit: versions
        
    script:
    """
#!/usr/bin/env python3
import os
import re
os.environ.setdefault('MPLCONFIGDIR', '/tmp/matplotlib')
from PIL import Image, ImageDraw, ImageFont
Image.MAX_IMAGE_PIXELS = None

PADDING = 10
BG_COLOR = (255, 255, 255, 255)
LABEL_HEIGHT = 350
LABEL_WIDTH = 350
AXIS_TITLE_HEIGHT = 100
FONT_PATH = "${projectDir}/fonts/DejaVuSans-Bold.ttf"
FONT_SIZE = 100

# md tag may be plain digits (legacy) or "<int>p<digits>" (current convention)
FILENAME_RE = re.compile(r"nn(\\d+)_md(\\d+p?\\d*)")


def parse_params(filename):
    \"\"\"Extract (nn, md_tag) from a plot filename. Returns None if unparseable.\"\"\"
    match = FILENAME_RE.search(os.path.basename(filename))
    if not match:
        return None
    nn = int(match.group(1))
    md_tag = match.group(2)
    return nn, md_tag


def md_tag_to_float(md_tag):
    return float(md_tag.replace('p', '.'))


def tile_images(image_dir, output_prefix, mc_values, padding=PADDING, bg_color=BG_COLOR):
    png_files = []
    for f in os.listdir(image_dir):
        if f.endswith(".png") and (f.startswith("umap_nn") or f.startswith("hdbscan_nn")):
            png_files.append(os.path.join(image_dir, f))
    print(f"Found {len(png_files)} PNG files")

    nn_values = set()
    md_tags = set()
    for filename in png_files:
        params = parse_params(filename)
        if params is None:
            print(f"Skipping unparseable filename: {filename}")
            continue
        nn, md_tag = params
        nn_values.add(nn)
        md_tags.add(md_tag)

    nn_values = sorted(nn_values)
    md_tags = sorted(md_tags, key=md_tag_to_float)
    mc_values = sorted(mc_values)

    font = ImageFont.truetype(FONT_PATH, FONT_SIZE)

    for md_tag in md_tags:
        md = md_tag_to_float(md_tag)
        images_by_row = {}

        for nn in nn_values:
            umap_path = os.path.join(image_dir, f"umap_nn{nn}_md{md_tag}.png")
            if os.path.exists(umap_path):
                images_by_row[nn] = [Image.open(umap_path)]

        for nn in list(images_by_row.keys()):
            umap_w, umap_h = images_by_row[nn][0].size
            for mc in mc_values:
                hdbscan_path = os.path.join(image_dir, f"hdbscan_nn{nn}_md{md_tag}_minclust{mc}.png")
                if os.path.exists(hdbscan_path):
                    img = Image.open(hdbscan_path)
                    if img.size != (umap_w, umap_h):
                        img = img.convert('RGB')
                        img = img.resize((umap_w, umap_h), Image.Resampling.LANCZOS)
                    images_by_row[nn].append(img)

        if not images_by_row:
            print(f"No images found for md_tag={md_tag}, skipping")
            continue

        ordered_rows = sorted(images_by_row.keys())
        rows = len(ordered_rows)
        cols = 1 + len(mc_values)

        max_width = max(row_images[0].size[0] for row_images in images_by_row.values())
        max_height = max(row_images[0].size[1] for row_images in images_by_row.values())

        grid_width = (max_width * cols) + (padding * (cols - 1))
        grid_height = (max_height * rows) + (padding * (rows - 1))

        total_width = grid_width + LABEL_WIDTH + padding * 2
        total_height = grid_height + LABEL_HEIGHT + AXIS_TITLE_HEIGHT * 2
        canvas = Image.new('RGBA', (total_width, total_height), bg_color)
        draw = ImageDraw.Draw(canvas)

        title = f"UMAP and HDBSCAN Clustering (md={md:.1f})"
        title_x = LABEL_WIDTH + grid_width // 2
        title_y = int(LABEL_HEIGHT * 0.18)
        draw.text((title_x, title_y), title, fill='black', font=font, anchor='mm')

        x_title_img = Image.new('RGBA', (grid_width, AXIS_TITLE_HEIGHT), bg_color)
        x_draw = ImageDraw.Draw(x_title_img)
        x_draw.text((grid_width // 2, AXIS_TITLE_HEIGHT // 2), "Minimum Cluster Size",
                     fill='black', font=font, anchor='mm')
        x_title_y = int(LABEL_HEIGHT * 0.32)
        canvas.paste(x_title_img, (LABEL_WIDTH, x_title_y))

        col_label_y = int(LABEL_HEIGHT * 0.78)
        draw.text((LABEL_WIDTH + max_width // 2, col_label_y), "UMAP",
                   fill='black', font=font, anchor='mm')
        for col, mc in enumerate(mc_values, 1):
            x = LABEL_WIDTH + col * (max_width + padding) + max_width // 2
            draw.text((x, col_label_y), f"MC={mc}", fill='black', font=font, anchor='mm')

        y_title_img = Image.new('RGBA', (grid_height, AXIS_TITLE_HEIGHT), bg_color)
        y_draw = ImageDraw.Draw(y_title_img)
        y_draw.text((grid_height // 2, AXIS_TITLE_HEIGHT // 2), "Nearest Neighbors Value (nn)",
                     fill='black', font=font, anchor='mm')
        y_title_img = y_title_img.rotate(90, expand=True)
        canvas.paste(y_title_img, (padding * 2, LABEL_HEIGHT + (grid_height // 2) - (y_title_img.height // 2)))

        for row, nn in enumerate(ordered_rows):
            label_img = Image.new('RGBA', (max_height, LABEL_WIDTH // 2), bg_color)
            label_draw = ImageDraw.Draw(label_img)
            label_draw.text((max_height // 2, LABEL_WIDTH // 4), f"nn={nn}",
                             fill='black', font=font, anchor='mm')
            label_img = label_img.rotate(90, expand=True)
            x = LABEL_WIDTH // 3
            y = LABEL_HEIGHT + row * (max_height + padding) + max_height // 2 - label_img.height // 2
            canvas.paste(label_img, (x, y), label_img)

        placeholder = Image.new('RGBA', (max_width, max_height), bg_color)
        for row_idx, nn in enumerate(ordered_rows):
            row_imgs = images_by_row.get(nn, [])
            for col_idx in range(cols):
                img = row_imgs[col_idx] if col_idx < len(row_imgs) else placeholder
                x = LABEL_WIDTH + col_idx * (max_width + padding)
                y = LABEL_HEIGHT + row_idx * (max_height + padding)
                canvas.paste(img, (x, y))

        output_file = f"{output_prefix}_md{md_tag}.png"
        canvas.save(output_file, optimize=True, quality=95, dpi=(600, 600))
        print(f"Saved tiled image for md={md:.1f}: {output_file}")


tile_images("plots", "plots/tiled_image", mc_values=${params.mc})

import platform
process_name = "${task.process}"
versions = [f'"{process_name}":', f'    python: {platform.python_version()}']
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write('\\n'.join(versions) + '\\n')
    """
}
