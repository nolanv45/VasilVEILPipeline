process PASV {
    tag "${meta.id}:${meta.protein}"
    label "process_medium"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/02_pasv",
        mode: 'copy',
        pattern: "pasv/output/*_putative.pasv_signatures.tsv",
        saveAs: { filename -> "${meta.id}/${meta.protein}/${filename}" }

    input:
        tuple val(meta), path(val_fasta, stageAs: 'val_full.fa'), path(unval_fasta, stageAs: 'unval_full.fa'), path(align_refs)

    output:
        tuple val(meta), path("pasv/output/${meta.protein}_putative.pasv_signatures.tsv"), emit: results
        path("versions.yml"), emit: versions

    script:
    """
    
    mkdir -p pasv/{input,output,pasv_tool}

    # Prepare input files
    cat "${val_fasta}" "${unval_fasta}" > "pasv/input/${meta.mapped_name}.fasta"
    cp "${align_refs}" "pasv/input/align_refs.fa"

    # Setup PASV
    PASV_DIR="pasv/pasv_tool"
    if [ ! -f "\${PASV_DIR}/pasv" ]; then
        mkdir -p "\${PASV_DIR}"
        cd "\${PASV_DIR}"
        wget -q https://github.com/mooreryan/pasv/releases/download/2.0.2/pasv-2.0.2-alpine-static.zip
        unzip pasv-2.0.2-alpine-static.zip
        chmod 755 pasv
        cd -
    fi

    # Run PASV
    \${PASV_DIR}/pasv msa \\
        --outdir=pasv/output \\
        --force \\
        --roi-start=${meta.roi_start} \\
        --roi-end=${meta.roi_end} \\
        --jobs=${task.cpus} \\
        --aligner=mafft \\
        pasv/input/${meta.mapped_name}.fasta \\
        pasv/input/align_refs.fa \\
        ${meta.cat_sites}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pasv: \$(\${PASV_DIR}/pasv --version | sed 's/pasv //')
    END_VERSIONS
    """
}