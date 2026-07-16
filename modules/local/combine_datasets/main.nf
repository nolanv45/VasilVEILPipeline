process COMBINE_DATASETS {
    // container "containers/phidra/phidra.sif"
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/phidra/phidra.yml"
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