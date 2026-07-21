process PASV_POST {
    tag "${meta.id}:${meta.protein}"
    label "process_single"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/03_annotation_analysis",
        mode: 'copy',
        pattern: "*.{tsv,png}",
        saveAs: { filename ->
            if (filename.endsWith('.tsv')) "pasv_annotation/${meta.id}/${meta.protein}/${filename}"
            else if (filename.endsWith('.png')) "pasv_annotation/${meta.id}/${meta.protein}/${filename}"
            else null
        }

    input:
        tuple val(meta), 
              path(pasv_file), 
              path(metadata_file),
              path(cleaned_fasta),
              path(filtered_tsv)

    output:
        tuple val(meta),
              path("${meta.id}_${meta.protein}_signature_stats.tsv"),
              path("${meta.id}_${meta.protein}_signature_distribution_mqc.png"),
              emit: results
              path("${meta.id}_${meta.protein}_signature_distribution_mqc.png"),
              emit: multiqc_plot
        path("versions.yml"), emit: versions
    script:
    """
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import math
import sys
import os
import re

def load_sequence_lengths(fasta_path):
    lengths = {}
    try:
        for record in SeqIO.parse(fasta_path, 'fasta'):
            lengths[record.id] = len(record.seq)
    except Exception as e:
        print(f"WARNING: unable to parse FASTA for ORF lengths ({fasta_path}): {e}", file=sys.stderr)
    return lengths

def load_genofeature_colors(metadata_path):
    if not metadata_path or not os.path.exists(metadata_path):
        print(f"WARNING: metadata file not found for colors: {metadata_path}", file=sys.stderr)
        return {}, {}
    try:
        meta = pd.read_csv(metadata_path, sep='\\t', dtype=str).fillna('')
        if 'genofeature' not in meta.columns or 'color' not in meta.columns:
            print("WARNING: metadata missing 'genofeature' or 'color' columns", file=sys.stderr)
            return {}, {}
        color_map = dict(zip(meta['genofeature'], meta['color']))
        color_map_lower = {k.lower(): v for k, v in color_map.items()}
        return color_map, color_map_lower
    except Exception as e:
        print(f"WARNING: unable to load metadata colors from {metadata_path}: {e}", file=sys.stderr)
        return {}, {}

raw_expected = '${meta.expected_sigs}'
tokens = re.findall(r"[A-Za-z0-9_]+", raw_expected or '')
expected_sigs = set(t.upper() for t in tokens if t)

try:
    threshold = float(${params.pasv_threshold})
except Exception:
    threshold = 0.01

# Load filtered orf_ids that passed APPLY_CRITERIA rules for this protein
filtered = pd.read_csv("${filtered_tsv}", sep='\\t', dtype=str, keep_default_na=False)
passed_ids = set(filtered[filtered["protein"] == "${meta.protein}"]["orf_id"])

# Read and process the PASV file, filtering to passed IDs only
print(f"Processing PASV file for ${meta.protein}")
df = pd.read_csv("${pasv_file}", sep='\\t')
df = df[df['name'].isin(passed_ids)].copy()
print(f"After criteria filter: {df.shape}")

if df.empty:
    print("WARNING: no sequences passed criteria for ${meta.id} ${meta.protein}", file=sys.stderr)
    pd.DataFrame(columns=['signature', 'span_class', 'count']).to_csv("${meta.id}_${meta.protein}_processed.tsv", sep='\\t', index=False)
    pd.DataFrame(columns=['signature', 'span_class', 'count']).to_csv("${meta.id}_${meta.protein}_signature_stats.tsv", sep='\\t', index=False)
    fig, ax = plt.subplots()
    ax.text(0.5, 0.5, 'No data passed criteria', ha='center', va='center')
    plt.savefig("${meta.id}_${meta.protein}_signature_distribution_mqc.png", bbox_inches='tight', dpi=300)
    plt.close()
    exit(0)

# Remove signatures with gaps
processed_df = df[~df['signature'].str.contains('-')].copy()

# Derive ORF lengths from actual amino-acid sequences
seq_lengths = load_sequence_lengths("${cleaned_fasta}")
processed_df['orf_length'] = processed_df['name'].map(seq_lengths)
processed_df = processed_df.dropna(subset=['orf_length']).copy()

color_map, color_map_lower = load_genofeature_colors("${metadata_file}")

# Create span classes
processed_df['span_class'] = processed_df.apply(
    lambda row: 'both'    if row['spans_start'] == 'Yes' and row['spans_end'] == 'Yes'
               else 'start'   if row['spans_start'] == 'Yes'
               else 'end'     if row['spans_end'] == 'Yes'
               else 'neither',
    axis=1
)

total_features = len(processed_df)
min_count = max(1, math.ceil(total_features * threshold))
print(f"Total features: {total_features}, min_count threshold: {min_count}", file=sys.stderr)

sig_counts = processed_df['signature'].value_counts().to_dict()
keep_sigs = set([s for s, c in sig_counts.items() if c >= min_count]) | expected_sigs
processed_df = processed_df[processed_df['signature'].isin(keep_sigs)].copy()

# Generate statistics
stats = processed_df.groupby(['signature', 'span_class']).agg({
    'orf_length': ['count', 'mean', 'std', 'min', 'max']
}).reset_index()
stats.columns = ['signature', 'span_class', 'count', 'mean_length', 'std_length', 'min_length', 'max_length']

signature_order = stats.groupby('signature')['count'].sum().sort_values(ascending=True).index

span_classes = sorted(processed_df['span_class'].unique())
n_signatures = len(signature_order)

fig_height = max(8, n_signatures * 0.4)
fig_width = len(span_classes) * 12

fig, axes = plt.subplots(1, len(span_classes),
                         figsize=(fig_width, fig_height),
                         sharey=True)

if len(span_classes) == 1:
    axes = [axes]

for i, span in enumerate(span_classes):
    ax = axes[i]
    span_data = processed_df[processed_df['span_class'] == span]

    if not span_data.empty:
        signature_palette = {
            sig: color_map.get(sig, color_map_lower.get(str(sig).lower(), '#808080'))
            for sig in signature_order
        }
        sns.violinplot(data=span_data,
                       x='orf_length',
                       y='signature',
                       order=signature_order,
                       hue='signature',
                       palette=signature_palette,
                       orient='h',
                       dodge=False,
                       inner='quartile',
                       linewidth=1,
                       cut=0,
                       ax=ax)

        legend = ax.get_legend()
        if legend is not None:
            legend.remove()

        span_stats = stats[stats['span_class'] == span]
        for idx, sig in enumerate(signature_order):
            sig_stats = span_stats[span_stats['signature'] == sig]
            if not sig_stats.empty:
                count = int(sig_stats['count'].iloc[0])
                mean = sig_stats['mean_length'].iloc[0]
                ax.text(
                    1.01, idx,
                    f"N={count}\\nμ={mean:.1f}",
                    transform=ax.get_yaxis_transform(),
                    ha='left', va='center',
                    fontsize=10, fontweight='bold',
                    clip_on=False,
                    bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='0.45', alpha=0.95)
                )

    ax.set_title(f"Span: {span}")
    ax.set_xlabel("ORF Length (aa)")
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='y', labelsize=10, pad=40)

plt.subplots_adjust(left=0.25, right=0.95, bottom=0.1, top=0.9, wspace=0.2)
fig.suptitle("${meta.id} ${meta.protein} Signature Distribution", fontsize=16, y=1.02)

stats.to_csv("${meta.id}_${meta.protein}_signature_stats.tsv", sep='\\t', index=False)
plt.savefig("${meta.id}_${meta.protein}_signature_distribution_mqc.png", bbox_inches='tight', dpi=300)
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