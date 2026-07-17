process PASV_ANNOTATION {
    tag "${meta.id}:${meta.protein}"
    label "process_single"
    conda "${moduleDir}/environment.yml"

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