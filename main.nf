nextflow.enable.dsl=2

process PHIDRA {
    conda "/home/nolanv/.conda/envs/phidra"

    output:
    path "${params.outdir}/"

    script:
    """
    source /etc/profile.d/conda.sh
    conda activate /home/nolanv/.conda/envs/phidra

    mkdir -p phidra

    # Copy the phidra directory contents to the task directory
    cp -r ${params.phidra_dir}/* phidra/

    # Debugging: List the contents of the phidra directory in the task's working directory
    ls -al phidra

    cd phidra
    python phidra_run.py \
        -i ${params.input_fasta} \
        -db ${params.subjectDB} \
        -pfam ${params.pfamDB} \
        -ida ${params.pfamDomain} \
        -f ${params.protein} \
        -o ${params.outdir}
    """
}

workflow {
    PHIDRA()
}