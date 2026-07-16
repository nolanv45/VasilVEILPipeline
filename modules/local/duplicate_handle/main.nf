process DUPLICATE_HANDLE {
    tag "${meta.id}"
    debug true
    // container "containers/phidra/phidra.sif"
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/phidra/phidra.yml"
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
