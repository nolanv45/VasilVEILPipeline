process PHIDRA_ONLY {
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
    """
    cp ${criteria_tsv} ${meta.protein}_annotation_criteria.tsv

import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    """
}
