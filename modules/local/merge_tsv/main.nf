process MERGE_TSV {
    tag "${meta.id}"
    label "process_single"
    conda "${moduleDir}/environment.yml"
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
