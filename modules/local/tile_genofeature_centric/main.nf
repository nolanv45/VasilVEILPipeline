process TILE_GENOFEATURE_CENTRIC {
    label "process_medium"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib_pillow:3f1442a23078f4ef' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib_pillow:e3f9f64a94d792de' }"

    publishDir "${params.outdir}/05_final_analysis/umap",
        mode: 'copy'
    
    input:
        path 'genofeature_plots/*'

    output:
        path "genofeature_highlighting_tiled_all_mqc.png", emit: tiled
        path "versions.yml", emit: versions

    script:
    """
#!/usr/bin/env python3
import os
import math
import platform
from PIL import Image
import PIL

def tile_images(image_dir, output_file, padding=10):
    png_files = [os.path.join(image_dir, f) for f in os.listdir(image_dir)
                 if f.lower().endswith("_mqc.png")]
    png_files.sort()
    if not png_files:
        print("No plot PNG files found for tiling.")
        return
    print(f"Tiling {len(png_files)} plot images...")
    images = [Image.open(f) for f in png_files]
    heights = [img.height for img in images]
    target_height = max(heights)
    resized_images = []
    resized_sizes = []
    for img in images:
        w, h = img.size
        new_w = int(w * (target_height / h))
        resized_img = img.resize((new_w, target_height), Image.LANCZOS)
        resized_images.append(resized_img)
        resized_sizes.append((new_w, target_height))
    n = len(resized_images)
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)
    col_widths = [0] * cols
    row_heights = [target_height] * rows
    for idx, (w, h) in enumerate(resized_sizes):
        row = idx // cols
        col = idx % cols
        if w > col_widths[col]:
            col_widths[col] = w
    x_offsets = [0]
    for w in col_widths[:-1]:
        x_offsets.append(x_offsets[-1] + w + padding)
    y_offsets = [0]
    for h in row_heights[:-1]:
        y_offsets.append(y_offsets[-1] + h + padding)
    canvas_width = sum(col_widths) + padding * (cols - 1)
    canvas_height = sum(row_heights) + padding * (rows - 1)
    canvas = Image.new('RGBA', (canvas_width, canvas_height), (255, 255, 255, 255))
    for idx, img in enumerate(resized_images):
        row = idx // cols
        col = idx % cols
        x = x_offsets[col]
        y = y_offsets[row]
        canvas.paste(img, (x, y))
    canvas.save(output_file)
    print(f"Tiled image saved as {output_file}")

tiled_output_file = "genofeature_highlighting_tiled_all_mqc.png"
tile_images("genofeature_plots", tiled_output_file)
print("All plots completed!")

process_name = "${task.process}"
with open("versions.yml", "w", encoding="utf-8") as versions_handle:
    versions_handle.write(f'"{process_name}":\\n')
    versions_handle.write(f'    python: {platform.python_version()}\\n')
    versions_handle.write(f'    pillow: {PIL.__version__}\\n')
    """
}