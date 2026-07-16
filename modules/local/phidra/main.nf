process PHIDRA {
    tag "${meta.id}:${meta.protein}"
    // container "containers/phidra/phidra.sif"
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/phidra/phidra.yml" 
    label "process_medium"
    publishDir "${params.outdir}/01_phidra",
        mode: 'copy',
        saveAs: { filename -> "${meta.id}/${filename}" }

    input:
        tuple val(meta), path(fasta), path(subject_db), path(pfamIDA)

    output:
        tuple val(meta), 
              path("${meta.protein}/${meta.protein}_validated_full_proteins.fa"), 
              path("${meta.protein}/${meta.protein}_pfam_validated_report.tsv"),
              path("${meta.protein}/${meta.protein}_unvalidated_full_proteins.fa"), 
              path("${meta.protein}/${meta.protein}_pfam_unvalidated_report.tsv"),
              path("*/mmseqs/initial/hits.tsv"),
              path("*/mmseqs/recursive/hits.tsv"),
              emit: results
        path("versions.yml"), emit: versions

    script:
    """
    WORK_DIR=\$PWD

    python ${workflow.projectDir}/bin/phidra/phidra_run.py \\
        -i ${fasta} \\
        -db ${subject_db} \\
        -pfam ${params.pfamDB} \\
        -ida ${pfamIDA} \\
        -f ${meta.protein} \\
        -o \$WORK_DIR \\
        -t ${task.cpus}
    
    if [ ! -s "\$WORK_DIR/${meta.protein}/mmseqs/recursive/hits.tsv" ]; then
        echo "recursive hits.tsv does not exist, creating it."
        mkdir -p "\$WORK_DIR/${meta.protein}/mmseqs/recursive"
        printf "Query_ID\tTarget_ID\tSequence_Identity\tAlignment_Length\tMismatches\tGap_Openings\tQuery_Start\tQuery_End\tTarget_Start\tTarget_End\tEvalue\tBit_Score\n" \
            > "\$WORK_DIR/${meta.protein}/mmseqs/recursive/hits.tsv"
    fi

    cp ${meta.protein}/final_results/validated_ida_pfams/full_proteins.fa \
    ${meta.protein}/${meta.protein}_validated_full_proteins.fa
    cp ${meta.protein}/final_results/validated_ida_pfams/pfam_validated_merged_report.tsv \
    ${meta.protein}/${meta.protein}_pfam_validated_report.tsv
    cp ${meta.protein}/final_results/unvalidated_ida_pfams/full_proteins.fa \
    ${meta.protein}/${meta.protein}_unvalidated_full_proteins.fa
    cp ${meta.protein}/final_results/unvalidated_ida_pfams/pfam_unvalidated_merged_report.tsv \
    ${meta.protein}/${meta.protein}_pfam_unvalidated_report.tsv

    python3 - <<'PY'
import platform

with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
PY
    """
}
