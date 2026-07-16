process ANALYZE_AND_PLOT {
    tag "${meta.id}:${meta.protein}"
    // container "containers/phidra/phidra.sif"
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/phidra/phidra.yml" 
    label "process_single"
    publishDir "${params.outdir}/03_annotation_analysis", 
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.tsv')) "phidra_annotation/${meta.id}/${meta.protein}/${filename}"
            else if (filename.endsWith('.png')) "phidra_annotation/${meta.id}/${meta.protein}/${filename}"
            else null
        }

    input:
        tuple val(meta), path(input_tsv), path(metadata_file), path(cleaned_fasta)

    output:
        tuple val(meta),
            path("${meta.id}_${meta.protein}_stats.tsv"),
            path("${meta.id}_${meta.protein}_length_distribution_mqc.png"),
            emit: results
            path("${meta.id}_${meta.protein}_length_distribution_mqc.png"),
            emit: multiqc_plot
        path("versions.yml"), emit: versions
    script:
    """
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import os

def load_sequence_lengths(fasta_path):
    lengths = {}
    try:
        for record in SeqIO.parse(fasta_path, 'fasta'):
            lengths[record.id] = len(record.seq)
    except Exception as e:
        print(f"WARNING: unable to parse FASTA for ORF lengths ({fasta_path}): {e}")
    return lengths

def load_genofeature_colors(metadata_path):
    if not metadata_path or not os.path.exists(metadata_path):
        print(f"WARNING: metadata file not found for colors: {metadata_path}")
        return {}, {}
    try:
        meta = pd.read_csv(metadata_path, sep='\\t', dtype=str).fillna('')
        if 'genofeature' not in meta.columns or 'color' not in meta.columns:
            print("WARNING: metadata missing 'genofeature' or 'color' columns")
            return {}, {}
        color_map = dict(zip(meta['genofeature'], meta['color']))
        color_map_lower = {k.lower(): v for k, v in color_map.items()}
        return color_map, color_map_lower
    except Exception as e:
        print(f"WARNING: unable to load metadata colors from {metadata_path}: {e}")
        return {}, {}

# Read and process input data
df = pd.read_csv("${input_tsv}", sep='\\t', dtype=str, keep_default_na=False)
df = df[df["protein"] == "${meta.protein}"].copy()
print(f"Processing data with {len(df)} entries")

# Always compute ORF length from sequence to avoid header naming inconsistencies.
seq_lengths = load_sequence_lengths("${cleaned_fasta}")
df['ORF_length'] = df['orf_id'].map(seq_lengths)

# Keep only valid numeric lengths for statistics and plotting.
df['ORF_length'] = pd.to_numeric(df['ORF_length'], errors='coerce')
plot_df = df.dropna(subset=['ORF_length']).copy()
color_map, color_map_lower = load_genofeature_colors("${metadata_file}")

# Generate basic statistics
stats = plot_df.groupby('genofeature').agg({
    'ORF_length': ['count', 'mean', 'min', 'max']
}).reset_index()

# Flatten column names
stats.columns = ['genofeature', 'count', 'mean_length', 'min_length', 'max_length']
stats = stats.sort_values('genofeature')

# Save statistics
stats.to_csv("${meta.id}_${meta.protein}_stats.tsv", sep='\\t', index=False)

# Create a violin plot with the same overall treatment used by PASV_POST.
fig, ax_plot = plt.subplots(
    1, 1,
    figsize=(14, max(6, len(stats) * 0.55))
)

if not plot_df.empty:
    ordered_features = stats['genofeature'].tolist()
    feature_palette = {
        feat: color_map.get(feat, color_map_lower.get(str(feat).lower(), '#808080'))
        for feat in ordered_features
    }
    sns.violinplot(
        data=plot_df,
        x='ORF_length',
        y='genofeature',
        order=ordered_features,
        hue='genofeature',
        palette=feature_palette,
        orient='h',
        dodge=False,
        inner='quartile',
        cut=0,
        linewidth=1,
        ax=ax_plot
    )
    legend = ax_plot.get_legend()
    if legend is not None:
        legend.remove()

    # Add a single boxed count/mean label per genofeature to mirror PASV_POST.
    for idx, (_, row) in enumerate(stats.iterrows()):
        count = int(row['count'])
        mean = row['mean_length']
        ax_plot.text(
            1.01,
            idx,
            f"N={count}\\nμ={mean:.1f}",
            transform=ax_plot.get_yaxis_transform(),
            ha='left',
            va='center',
            fontsize=10,
            fontweight='bold',
            clip_on=False,
            bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='0.45', alpha=0.95)
        )
else:
    ax_plot.text(
        0.5,
        0.5,
        'No valid ORF lengths available for plotting',
        ha='center',
        va='center',
        transform=ax_plot.transAxes
    )

ax_plot.set_title("Length Distribution by genofeature")
ax_plot.set_xlabel("ORF Length (aa)")
ax_plot.set_ylabel("genofeature")
ax_plot.grid(True, axis='x', linestyle='--', alpha=0.7)
ax_plot.tick_params(axis='y', labelsize=10, pad=28)

fig.tight_layout()
fig.savefig("${meta.id}_${meta.protein}_length_distribution_mqc.png", bbox_inches='tight', dpi=300)
plt.close()

import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    handle.write(f'    matplotlib: {plt.matplotlib.__version__}\\n')
    handle.write(f'    seaborn: {sns.__version__}\\n')
    import Bio
    handle.write(f'    biopython: {Bio.__version__}\\n')
    """
}
