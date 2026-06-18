nextflow.enable.dsl=2

include { VEILREPMOD } from './workflows/veilrepmod'
include { PIPELINE_INITIALIZATION } from './subworkflows/local/pipeline_initialization'
include { PIPELINE_COMPLETION } from './subworkflows/local/pipeline_completion'

// WORKFLOW: Run main analysis pipeline
workflow {
    main:

    PIPELINE_INITIALIZATION(


    )



    VEILREPMOD()

    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,

    )


}