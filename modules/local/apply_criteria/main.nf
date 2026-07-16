process APPLY_CRITERIA {
    tag "${meta.id}"
    // container "containers/phidra/phidra.sif"
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/phidra/phidra.yml"
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
