process DOMAIN_MATCH_ANNOTATION {
    tag "${meta.id}:${meta.protein}"
    // container "containers/phidra/phidra.sif"
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/phidra/phidra.yml"
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