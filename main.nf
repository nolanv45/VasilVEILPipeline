nextflow.enable.dsl=2

process PHIDRA {
    input:
    val

    output:
    stdout

    script:
    """
    conda activate

    for p in ${PROTEINS}
    do
    DB=${p}_db
    IDA=${p}_ida

    python phidra_run.py \
    -i ${QUERY_FILE} \
    -db ${!DB} \
    -pfam ${PFAM} \
    -ida ${IDA_DIR}/${!IDA} \
    -f ${p} \
    -o ${OUTDIR} \
    -t 18

    done

    conda deactivate
    """
}

workflow {
    sayHello()
}
