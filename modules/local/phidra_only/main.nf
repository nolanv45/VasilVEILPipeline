process PHIDRA_ONLY {
    tag "${meta.id}:${meta.protein}"
    label "process_single"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas:2a5ca2e7dd4ced9c' :
        'community.wave.seqera.io/library/python_pandas:eb68d9296e3f036a' }"

    input:
        tuple val(meta), path(criteria_tsv)

    output:
        tuple val(meta), path("${meta.protein}_annotation_criteria.tsv"), emit: results
        path("versions.yml"), emit: versions
    script:
    """
    cp ${criteria_tsv} ${meta.protein}_annotation_criteria.tsv

    python3 - <<'PY'
    import platform
    with open('versions.yml', 'w', encoding='utf-8') as handle:
        handle.write(f'"${task.process}":\\n')
        handle.write(f'    python: {platform.python_version()}\\n')
    PY
    """
}