process PASV_ANNOTATION {
    tag "${meta.id}:${meta.protein}"
    label "process_single"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas:2a5ca2e7dd4ced9c' :
        'community.wave.seqera.io/library/python_pandas:eb68d9296e3f036a' }"

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
import math
import re

df = pd.read_csv("${criteria_tsv}", sep="\\t", dtype=str, keep_default_na=False)
df_pasv = pd.read_csv("pasv_signatures.tsv", sep="\\t", dtype=str, keep_default_na=False)

pasv_map = df_pasv.set_index("name")[["signature", "spans"]].to_dict(orient="index")

raw_expected = '${meta.expected_sigs}'
tokens = re.findall(r"[A-Za-z0-9_]+", raw_expected or '')
expected_sigs = set(t.upper() for t in tokens if t)

try:
    threshold = float(${params.pasv_threshold})
except Exception:
    threshold = 0.01

filtered_pasv = df_pasv[~df_pasv['signature'].astype(str).str.contains('-')].copy()
total_features = len(filtered_pasv)
min_count = max(1, math.ceil(total_features * threshold)) if total_features else 0
sig_counts = filtered_pasv['signature'].value_counts().to_dict()
keep_sigs = set([s for s, c in sig_counts.items() if c >= min_count]) | expected_sigs

mask = df["Protein"] == "${meta.protein}"

for idx in df[mask].index:
    query_id = df.at[idx, "Query_ID"]
    if query_id in pasv_map:
        signature = pasv_map[query_id]["signature"]
        if signature not in keep_sigs:
            continue
        df.at[idx, "PASV"] = "Yes"
        df.at[idx, "Genofeature"] = signature
        df.at[idx, "PASV_Spans"] = pasv_map[query_id]["spans"]

df.to_csv("${meta.protein}_annotation_criteria.tsv", sep="\\t", index=False)
import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    """
}