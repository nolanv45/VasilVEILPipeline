nextflow.enable.dsl=2

process PHIDRA {
    conda "/home/nolanv/.conda/envs/phidra"


    input:
    val p

    output:
    path "${params.outdir}/"

    script:
    """
    source /etc/profile.d/conda.sh
    conda activate /home/nolanv/.conda/envs/phidra
    python /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/phidra/phidra_run.py \
    -i ${params.input_fasta} \
    -db ${params.subjectDB} \
    -pfam ${params.pfamDB} \
    -ida ${params.pfamDomain} \
    -f ${p} \
    -o ${params.outdir} 

    """
}

workflow {
    Channel.from(params.proteins)
        | PHIDRA // Send each protein to PHIDRA process
}
