process CRITERIA_TSV {
    tag "${meta.id}"
    // container "containers/phidra/phidra.sif"
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/phidra/phidra.yml"
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