process TOP_HIT_ANNOTATION {
    tag "${meta.id}:${meta.protein}"
    // container "containers/phidra/phidra.sif"
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/phidra/phidra.yml"
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