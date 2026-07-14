nextflow.enable.dsl=2

process CLEAN_FASTA_HEADERS {
    tag "${meta.id}"
    container "containers/phidra/phidra.sif"
    label "process_single"

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("${meta.id}_cleaned.fasta"), emit: fasta
        path("versions.yml"), emit: versions
    script:
    """
    awk 'BEGIN{RS=">"; ORS=""} NR>1{n=index(\$0, "\\n"); header=substr(\$0,1,n-1); sub(/[[:space:]].*/, "", header); gsub(/\\./, "-", header); gsub(/[^A-Za-z0-9_-]/, "", header); if (header == "") header = "seq_" NR; seq=substr(\$0,n+1); gsub(/\\n/, "", seq); gsub(/\\*/, "", seq); print ">"header"\\n"seq"\\n"}' ${fasta} > "${meta.id}_cleaned.fasta"

    python3 - <<'PY'
import subprocess

awk_version = subprocess.check_output("awk --version | head -n1 | sed 's/.*Awk //; s/,.*//'", shell=True, text=True).strip()

with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    awk: {awk_version}\\n')
PY
    """
}

process PHIDRA {
    tag "${meta.id}:${meta.protein}"
    container "containers/phidra/phidra.sif"
    label "process_medium"
    publishDir "${params.outdir}/01_phidra",
        mode: 'copy',
        saveAs: { filename -> "${meta.id}/${filename}" }

    input:
        tuple val(meta), path(fasta), path(subject_db), path(pfamIDA)

    output:
        tuple val(meta), 
              path("${meta.protein}/${meta.protein}_validated_full_proteins.fa"), 
              path("${meta.protein}/${meta.protein}_pfam_validated_report.tsv"),
              path("${meta.protein}/${meta.protein}_unvalidated_full_proteins.fa"), 
              path("${meta.protein}/${meta.protein}_pfam_unvalidated_report.tsv"),
              path("*/mmseqs/initial/hits.tsv"),
              path("*/mmseqs/recursive/hits.tsv"),
              emit: results
        path("versions.yml"), emit: versions

    script:
    """
    WORK_DIR=\$PWD

    python ${workflow.projectDir}/bin/phidra/phidra_run.py \\
        -i ${fasta} \\
        -db ${subject_db} \\
        -pfam ${params.pfamDB} \\
        -ida ${pfamIDA} \\
        -f ${meta.protein} \\
        -o \$WORK_DIR \\
        -t ${task.cpus}
    
    if [ ! -s "\$WORK_DIR/${meta.protein}/mmseqs/recursive/hits.tsv" ]; then
        echo "recursive hits.tsv does not exist, creating it."
        mkdir -p "\$WORK_DIR/${meta.protein}/mmseqs/recursive"
        printf "Query_ID\tTarget_ID\tSequence_Identity\tAlignment_Length\tMismatches\tGap_Openings\tQuery_Start\tQuery_End\tTarget_Start\tTarget_End\tEvalue\tBit_Score\n" \
            > "\$WORK_DIR/${meta.protein}/mmseqs/recursive/hits.tsv"
    fi

    cp ${meta.protein}/final_results/validated_ida_pfams/full_proteins.fa \
    ${meta.protein}/${meta.protein}_validated_full_proteins.fa
    cp ${meta.protein}/final_results/validated_ida_pfams/pfam_validated_merged_report.tsv \
    ${meta.protein}/${meta.protein}_pfam_validated_report.tsv
    cp ${meta.protein}/final_results/unvalidated_ida_pfams/full_proteins.fa \
    ${meta.protein}/${meta.protein}_unvalidated_full_proteins.fa
    cp ${meta.protein}/final_results/unvalidated_ida_pfams/pfam_unvalidated_merged_report.tsv \
    ${meta.protein}/${meta.protein}_pfam_unvalidated_report.tsv

    python3 - <<'PY'
import platform

with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
PY
    """
}


process CRITERIA_TSV {
    tag "${meta.id}"
    container "containers/phidra/phidra.sif"
    label "process_single"

    input:
        // Changed path to val to allow safe passing of empty lists/placeholders
        tuple val(meta), path(all_files)

    output:
        tuple val(meta), path("${meta.id}_annotation_criteria.tsv"), emit: results
        path("versions.yml"), emit: versions

    script:
    def file_list_str = all_files.collect { "\"${it}\"" }.join(", ")
    """
#!/usr/bin/env python3
import pandas as pd
import os

file_list = [${file_list_str}]

validated_files   = [f for f in file_list if f.endswith("_pfam_validated_report.tsv")]
unvalidated_files = [f for f in file_list if f.endswith("_pfam_unvalidated_report.tsv")]

all_phidra_dfs = []

for f in validated_files:
    if os.path.exists(f):
        df = pd.read_csv(f, sep='\\t')
        df['PHIDRA'] = "Yes"
        df['PASV'] = "No"
        df['PASV_Spans'] = ""
        df['Dataset'] = "${meta.id}"
        df['Protein'] = os.path.basename(f).split("_", 1)[0]
        df['Genofeature'] = os.path.basename(f).split("_", 1)[0]
        all_phidra_dfs.append(df)

for f in unvalidated_files:
    if os.path.exists(f):
        df = pd.read_csv(f, sep='\\t')
        df['PHIDRA'] = "No"
        df['PASV'] = "No"
        df['PASV_Spans'] = ""
        df['Dataset'] = "${meta.id}"
        df['Protein'] = os.path.basename(f).split("_", 1)[0]
        df['Genofeature'] = os.path.basename(f).split("_", 1)[0]
        all_phidra_dfs.append(df)

combined_df = pd.concat(all_phidra_dfs, ignore_index=True)

combined_df.to_csv('${meta.id}_annotation_criteria.tsv', sep='\\t', index=False)
import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
"""
}

process DOMAIN_MATCH_ANNOTATION {
    tag "${meta.id}:${meta.protein}"
    container "containers/phidra/phidra.sif"
    label "process_single"

    input:
        tuple val(meta), path(criteria_tsv)

    output:
        tuple val(meta), path("${meta.protein}_annotation_criteria.tsv"), emit: results
        path("versions.yml"), emit: versions

    script:
    def map_json = groovy.json.JsonOutput.toJson(params.pfam_annotation_map)

    """
#!/usr/bin/env python3
import json
import pandas as pd

df = pd.read_csv("${criteria_tsv}", sep="\\t", dtype=str, keep_default_na=False)

pfam_map = json.loads('${map_json}')

mask = df["Protein"] == "${meta.protein}"

for idx in df[mask].index:
    pfams = df.at[idx, "Pfam_IDs"]
    if not pfams:
        continue
    pf_list = [p.strip().upper() for p in pfams.split("|") if p.strip()]
    for pf in pf_list:
        if pf in pfam_map:
            df.at[idx, "Genofeature"] = pfam_map[pf]
            break

df.to_csv("${meta.protein}_annotation_criteria.tsv", sep="\\t", index=False)
import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    """
}

process PASV {
    tag "${meta.id}:${meta.protein}"
    container "containers/phidra/phidra.sif"
    label "process_medium"
    publishDir "${params.outdir}/02_pasv",
        mode: 'copy',
        pattern: "pasv/output/*_putative.pasv_signatures.tsv",
        saveAs: { filename -> "${meta.id}/${meta.protein}/${filename}" }

    input:
        tuple val(meta), path(val_fasta, stageAs: 'val_full.fa'), path(unval_fasta, stageAs: 'unval_full.fa'), path(align_refs)

    output:
        tuple val(meta), path("pasv/output/${meta.protein}_putative.pasv_signatures.tsv"), emit: results
        path("versions.yml"), emit: versions

    script:
    """
    
    mkdir -p pasv/{input,output,pasv_tool}

    # Prepare input files
    cat "${val_fasta}" "${unval_fasta}" > "pasv/input/${meta.mapped_name}.fasta"
    cp "${align_refs}" "pasv/input/align_refs.fa"

    # Setup PASV
    PASV_DIR="pasv/pasv_tool"
    if [ ! -f "\${PASV_DIR}/pasv" ]; then
        mkdir -p "\${PASV_DIR}"
        cd "\${PASV_DIR}"
        wget -q https://github.com/mooreryan/pasv/releases/download/2.0.2/pasv-2.0.2-alpine-static.zip
        unzip pasv-2.0.2-alpine-static.zip
        chmod 755 pasv
        cd -
    fi

    # Run PASV
    \${PASV_DIR}/pasv msa \\
        --outdir=pasv/output \\
        --force \\
        --roi-start=${meta.roi_start} \\
        --roi-end=${meta.roi_end} \\
        --jobs=${task.cpus} \\
        --aligner=mafft \\
        pasv/input/${meta.mapped_name}.fasta \\
        pasv/input/align_refs.fa \\
        ${meta.cat_sites}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pasv: \$(\${PASV_DIR}/pasv --version | sed 's/pasv //')
    END_VERSIONS
    """
}

process PASV_ANNOTATION {
    tag "${meta.id}:${meta.protein}"
    container "containers/phidra/phidra.sif"
    label "process_single"

    input:
        tuple val(meta),
              path(criteria_tsv),
              path(pasv_signatures, stageAs: 'pasv_signatures.tsv')

    output:
        tuple val(meta), path("${meta.protein}_annotation_criteria.tsv"), emit: results
        path("versions.yml"), emit: versions
    script:
    """
#!/usr/bin/env python3
import pandas as pd

df = pd.read_csv("${criteria_tsv}", sep="\\t", dtype=str, keep_default_na=False)
df_pasv = pd.read_csv("pasv_signatures.tsv", sep="\\t", dtype=str, keep_default_na=False)

pasv_map = df_pasv.set_index("name")[["signature", "spans"]].to_dict(orient="index")

mask = df["Protein"] == "${meta.protein}"

for idx in df[mask].index:
    query_id = df.at[idx, "Query_ID"]
    if query_id in pasv_map:
        df.at[idx, "PASV"] = "Yes"
        df.at[idx, "Genofeature"] = pasv_map[query_id]["signature"]
        df.at[idx, "PASV_Spans"] = pasv_map[query_id]["spans"]

df.to_csv("${meta.protein}_annotation_criteria.tsv", sep="\\t", index=False)
import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    """
}

process ANALYZE_AND_PLOT {
    tag "${meta.id}:${meta.protein}"
    container "containers/phidra/phidra.sif"
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



process PASV_POST {
    tag "${meta.id}:${meta.protein}"
    container "containers/phidra/phidra.sif"
    label "process_single"
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


process TOP_HIT_ANNOTATION {
    tag "${meta.id}:${meta.protein}"
    container "containers/phidra/phidra.sif"
    label "process_single"

    input:
    tuple val(meta), 
          path(criteria_tsv),
          path(initial_search, stageAs: 'initial_search.tsv'),
          path(recursive_search, stageAs: 'recursive_search.tsv')

    output:
    tuple val(meta), path("${meta.protein}_annotation_criteria.tsv"), emit: results
    path("versions.yml"), emit: versions
    script:
    """
#!/usr/bin/env python3
import os
import pandas as pd

df = pd.read_csv("${criteria_tsv}", sep="\\t", dtype=str, keep_default_na=False)

def build_hit_map(filepath):
    if not os.path.exists(filepath):
        return {}
    hits = pd.read_csv(filepath, sep="\\t", dtype=str, keep_default_na=False)
    return {
        row.iloc[0]: row.iloc[1].split("_")[0]
        for _, row in hits.iterrows()
    }

def build_raw_map(filepath):
    if not os.path.exists(filepath):
        return {}
    hits = pd.read_csv(filepath, sep="\\t", dtype=str, keep_default_na=False)
    return {
        row.iloc[0]: row.iloc[1]
        for _, row in hits.iterrows()
    }

# initial: query_id -> genofeature (first part of target_id)
initial_map = build_hit_map("initial_search.tsv")

# recursive raw: query_id -> full target_id (which is a query_id from initial)
recursive_raw = build_raw_map("recursive_search.tsv")

# recursive final: query_id -> genofeature (via initial lookup)
recursive_map = {
    query_id: initial_map[target_id]
    for query_id, target_id in recursive_raw.items()
    if target_id in initial_map
}

mask = df["Protein"] == "${meta.protein}"

for idx in df[mask].index:
    query_id = df.at[idx, "Query_ID"]
    if query_id in initial_map:
        df.at[idx, "Genofeature"] = initial_map[query_id]
    elif query_id in recursive_map:
        df.at[idx, "Genofeature"] = recursive_map[query_id]

df.to_csv("${meta.protein}_annotation_criteria.tsv", sep="\\t", index=False)
import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    """
}

process DUPLICATE_HANDLE {
    tag "${meta.id}"
    debug true
    container "containers/phidra/phidra.sif"
    label "process_single"

    publishDir "${params.outdir}/03_annotation_analysis",
        mode: 'copy',
        pattern: "*.tsv",
        saveAs: { filename -> "${meta.id}/${filename}" }

    input:
        tuple val(meta), path(input_file)    // Properly declare the input file

    output:
        tuple val(meta), path("duplicate_orfs.tsv")   // Declare the output file
        path("versions.yml"), emit: versions
    script:
    """
#!/usr/bin/env python3
import pandas as pd

# Read the input file with proper quoting and correct separator
df = pd.read_csv("${input_file}", sep='\\t')

# Find duplicates based on both contig_id and orf_id
duplicates = df[df.duplicated(subset=['contig_id', 'orf_id'], keep=False)]

# Sort duplicates for better readability
duplicates = duplicates.sort_values(['contig_id', 'orf_id'])

# Save duplicates with tab separator
duplicates.to_csv('duplicate_orfs.tsv', sep='\\t', index=False)

# Print summary for debugging
print(f"Found {len(duplicates)} duplicate entries")

import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    """
}


process COMBINE_DATASETS {
    container "containers/phidra/phidra.sif"
    label "process_single"
    publishDir "${params.outdir}",
        mode: 'copy',
        pattern: "*.tsv"

    input:
        path(tsv_files)

    output:
        path "combined_datasets.tsv", emit: results
        path("versions.yml"), emit: versions
    script:
    """
#!/usr/bin/env python3
import pandas as pd
import os

dfs = []
input_files = "${tsv_files}".split()

for tsv_file in input_files:
    if os.path.exists(tsv_file):
        print(f"Reading file: {tsv_file}")
        try:
            df = pd.read_csv(tsv_file, sep='\\t')
            if not df.empty:
                dfs.append(df)
                print(f"Added {len(df)} rows from {tsv_file}")
            else:
                print(f"Warning: Empty DataFrame from {tsv_file}")
        except Exception as e:
            print(f"Error reading {tsv_file}: {str(e)}")
    else:
        print(f"Warning: File not found: {tsv_file}")

if dfs:
    combined = pd.concat(dfs, ignore_index=True)
    print(f"Combined DataFrame has {len(combined)} rows")
    combined.to_csv("combined_datasets.tsv", sep='\\t', index=False)
else:
    print("No data to combine, creating empty output")
    pd.DataFrame(columns=['contig_id', 'orf_id', 'genofeature', 'protein', 'dataset']).to_csv("combined_datasets.tsv", sep='\\t', index=False)

import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    """
}

process MERGE_TSV {
    tag "${meta.id}"
    container "containers/phidra/phidra.sif"
    label "process_single"
    publishDir "${params.outdir}/03_annotation_analysis",
        mode: 'copy',
        saveAs: { filename -> "${meta.id}/${filename}" }

    input:
        tuple val(meta), path(tsvs)

    output:
        tuple val(meta), path("${meta.id}_annotation_criteria.tsv"), emit: results
        path("versions.yml"), emit: versions
    script:
    def tsv_list = tsvs.collect { "\"${it}\"" }.join(", ")
    """
#!/usr/bin/env python3
import pandas as pd

tsv_files = [${tsv_list}]

dfs = [pd.read_csv(f, sep="\\t", dtype=str, keep_default_na=False) for f in tsv_files]

# Build a map of protein -> df that edited it
# Each file is named {protein}_annotation_criteria.tsv
import os
protein_df_map = {}
for f, df in zip(tsv_files, dfs):
    protein = os.path.basename(f).split("_")[0]
    protein_df_map[protein] = df

# Use first df as base structure
base = dfs[0].set_index("Query_ID")

# For each protein, take only that protein's rows from the authoritative df
for protein, df in protein_df_map.items():
    df = df.set_index("Query_ID")
    protein_rows = df[df["Protein"] == protein]
    base.update(protein_rows)

base.reset_index().to_csv("${meta.id}_annotation_criteria.tsv", sep="\\t", index=False)

import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    """
}

process PHIDRA_ONLY {
    tag "${meta.id}:${meta.protein}"
    container "containers/phidra/phidra.sif"
    label "process_single"

    input:
        tuple val(meta), path(criteria_tsv)

    output:
        tuple val(meta), path("${meta.protein}_annotation_criteria.tsv"), emit: results
        path("versions.yml"), emit: versions
    script:
    """
    cp ${criteria_tsv} ${meta.protein}_annotation_criteria.tsv

import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    """
}


process APPLY_CRITERIA {
    tag "${meta.id}"
    container "containers/phidra/phidra.sif"
    label "process_single"
    publishDir "${params.outdir}/03_annotation_analysis",
        mode: 'copy',
        saveAs: { filename -> "${meta.id}/${filename}" }

    input:
        tuple val(meta), path(criteria_tsv)

    output:
        tuple val(meta), path("${meta.id}_genofeatures.tsv"), emit: results
        path("versions.yml"), emit: versions
    script:
    def pasv_proteins_json = groovy.json.JsonOutput.toJson(params.pasv_proteins ?: [])
    def orf_drop = params.orf_coord_fields[meta.id] ?: -1

    """
#!/usr/bin/env python3
import json
import pandas as pd

df = pd.read_csv("${criteria_tsv}", sep="\\t", dtype=str, keep_default_na=False)

pasv_proteins = set(json.loads('${pasv_proteins_json}'))
orf_drop = ${orf_drop}  # negative index, e.g. -3 means drop last 3 fields

# Rule 1: filter out anything PHIDRA == No
df = df[df["PHIDRA"] == "Yes"]

# Rule 2: for PASV proteins, filter out PASV == No
pasv_mask = df["Protein"].isin(pasv_proteins)
df = df[~pasv_mask | (df["PASV"] == "Yes")]

# Parse contig_id and orf_id from Query_ID using orf_coord_fields
df["orf_id"] = df["Query_ID"]
df["contig_id"] = df["Query_ID"].apply(
    lambda x: "_".join(x.split("_")[:orf_drop])
)

# Output with required columns
df[["contig_id", "orf_id", "Genofeature", "Protein", "Dataset"]].rename(columns={
    "Protein":    "protein",
    "Dataset":    "dataset",
    "Genofeature": "genofeature"
}).to_csv("${meta.id}_genofeatures.tsv", sep="\\t", index=False)

import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    """
}



workflow ANNOTATE_PROTEINS {
    take:
        ch_dataset

    main:
        ch_metadata = channel.fromPath(params.genofeature_metadata)

        // Normalize fasta input file
        CLEAN_FASTA_HEADERS(ch_dataset)
         
        // Prepare PHIDRA input
        ch_dataset_proteins = CLEAN_FASTA_HEADERS.out.fasta
            .combine(channel.fromList(params.proteins))
            .map { meta, fasta, protein_config -> 
                def new_meta = meta + [
                    protein: protein_config.protein,
                    pfamDomain: protein_config.pfamDomain,
                    subjectDB: protein_config.subjectDB,
                    pasv_align_refs: protein_config.containsKey('pasv_align_refs') ? 
                        protein_config.pasv_align_refs : null,
                    cleaned_fasta: fasta
                ]
                // stage the subject DB file so Nextflow copies it into the task workdir
                tuple(new_meta, fasta, file(protein_config.subjectDB), file(protein_config.pfamDomain))
            }
        PHIDRA(ch_dataset_proteins)

        // ----------------------------------------------------------------------
        // CRITERIA TSV
        // ----------------------------------------------------------------------
        ch_phidra_crit = PHIDRA.out.results

        ch_unval_crit = ch_phidra_crit.map { meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs -> [meta.id, unval_pfam] }
        ch_val_crit   = ch_phidra_crit.map { meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs -> [meta.id, val_pfam] }

        ch_criteria_input = ch_unval_crit
            .mix(ch_val_crit)
            .groupTuple(by: 0)
            .map { id, files -> [ [id: id], files.flatten() ] }

        CRITERIA_TSV(ch_criteria_input)

        // ----------------------------------------------------------------------
        // JOIN CRITERIA TSV WITH PHIDRA RESULTS
        // ----------------------------------------------------------------------
        ch_tsv = CRITERIA_TSV.out.results.map { meta, tsv -> [meta.id, tsv] }

        ch_phidra_keyed = PHIDRA.out.results
            .map { meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                [meta.id, meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs]
            }
 
        ch_phidra_with_tsv = PHIDRA.out.results
            .map { meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                [meta.id, meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs]
            }
            .combine(ch_tsv, by: 0)
            .map { id, meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs, tsv ->
                [meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs]
            }

        ch_branched = ch_phidra_with_tsv
            .branch { meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                tophit: meta.protein.toLowerCase() in params.tophit_proteins.collect { it.toLowerCase() }
                pasv:   meta.protein.toLowerCase() in params.pasv_proteins.collect { it.toLowerCase() }
                domain: meta.protein.toLowerCase() in params.domain_proteins.collect { it.toLowerCase() }
                other: true
            }


        ch_branched.other
            .map { meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                [meta, tsv]
            }
            | PHIDRA_ONLY

        ch_branched.tophit
            .map { meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                [meta, tsv, init, recurs]
            }
            | TOP_HIT_ANNOTATION

        ch_pasv_input = ch_branched.pasv
            .map { meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                def settings      = (params.pasv_settings ?: [:]).get(meta.protein, [:])
                def mapped_name   = settings.mapped_name   ?: "${meta.protein}_putative"
                def roi_start     = settings.roi_start     ?: ''
                def roi_end       = settings.roi_end       ?: ''
                def cat_sites     = settings.cat_sites     ?: ''
                def expected_sigs = settings.expected_sigs ?: ''
                def new_meta = meta + [
                    mapped_name:   mapped_name,
                    roi_start:     roi_start,
                    roi_end:       roi_end,
                    cat_sites:     cat_sites,
                    expected_sigs: expected_sigs
                ]
                tuple(new_meta, val_fasta, unval_fasta, file(meta.pasv_align_refs))
            }

        ch_pasv_input | PASV

        ch_pasv_tsv = ch_branched.pasv
            .map { meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                ["${meta.id}:${meta.protein}", tsv]
            }

        PASV.out.results
            .map { meta, pasv_sig -> ["${meta.id}:${meta.protein}", meta, pasv_sig] }
            .join(ch_pasv_tsv, by: 0)
            .map { key, meta, pasv_sig, tsv -> [meta, tsv, pasv_sig] }
            | PASV_ANNOTATION

        ch_branched.domain
            .map { meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                [meta, tsv]
            }
            | DOMAIN_MATCH_ANNOTATION



        // ----------------------------------------------------------------------
        // MERGE ALL PER-PROTEIN TSVs BACK INTO ONE PER ID
        // ----------------------------------------------------------------------
        TOP_HIT_ANNOTATION.out.results
            .mix(PASV_ANNOTATION.out.results, DOMAIN_MATCH_ANNOTATION.out.results, PHIDRA_ONLY.out.results)
            .map { meta, tsv -> [meta.id, tsv] }
            .groupTuple(by: 0)
            .map { id, tsvs -> [ [id: id], tsvs.flatten() ] }
            | MERGE_TSV


        // APPLIES RULESET BELOW
        MERGE_TSV.out.results | APPLY_CRITERIA

        ch_cleaned_fasta = CLEAN_FASTA_HEADERS.out.fasta
            .map { meta, fasta -> [meta.id, fasta] }

        ch_analyze_input = PHIDRA.out.results
            .filter { meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                !(meta.protein.toLowerCase() in params.pasv_proteins.collect { it.toLowerCase() })
            }
            .map { meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs -> [meta.id, meta] }
            .combine(
                APPLY_CRITERIA.out.results.map { meta, tsv -> [meta.id, tsv] },
                by: 0
            )
            .combine(ch_cleaned_fasta, by: 0)
            .combine(ch_metadata)
            .map { id, meta, tsv, cleaned_fasta, metadata_file ->
                tuple(meta, tsv, metadata_file, cleaned_fasta)
            }

        ch_analyze_input
            | ANALYZE_AND_PLOT



        PASV.out.results
            .map { meta, pasv_sig -> [meta.id, meta, pasv_sig] }
            .combine(
                APPLY_CRITERIA.out.results.map { meta, tsv -> [meta.id, tsv] },
                by: 0
            )
            .combine(ch_cleaned_fasta, by: 0)
            .combine(ch_metadata)
            .map { id, meta, pasv_sig, filtered_tsv, cleaned_fasta, metadata_file ->
                tuple(meta, pasv_sig, metadata_file, cleaned_fasta, filtered_tsv)
            }
            | PASV_POST



        ch_combine_datasets = APPLY_CRITERIA.out.results
            .map { meta, file -> file.toAbsolutePath() }
            .collect()

        COMBINE_DATASETS(ch_combine_datasets)

        DUPLICATE_HANDLE(APPLY_CRITERIA.out.results)

        ch_annotation_plots = ANALYZE_AND_PLOT.out.multiqc_plot
            .mix(PASV_POST.out.multiqc_plot)

        ch_versions = PHIDRA.out.versions
            .mix(CRITERIA_TSV.out.versions)
            .mix(TOP_HIT_ANNOTATION.out.versions)
            .mix(PASV_ANNOTATION.out.versions)
            .mix(DOMAIN_MATCH_ANNOTATION.out.versions)
            .mix(PHIDRA_ONLY.out.versions)
            .mix(MERGE_TSV.out.versions)
            .mix(APPLY_CRITERIA.out.versions)
            .mix(ANALYZE_AND_PLOT.out.versions)
            .mix(PASV.out.versions)
            .mix(PASV_POST.out.versions)
            .mix(COMBINE_DATASETS.out.versions)


    emit:
        ch_combined_tsv = COMBINE_DATASETS.out.results
        multiqc_files = ch_annotation_plots
        versions = ch_versions
}



