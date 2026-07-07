nextflow.enable.dsl=2

// SUBWORKFLOW: Three subworkflows consisting of local modules
include { REP_MODULE_ANALYSIS } from './subworkflows/rep_module_analysis'
include { ANNOTATE_PROTEINS } from './subworkflows/annotate_proteins'
include { EMBEDDING_PARAMETER_DECISION } from './subworkflows/embedding_parameter_decision'

workflow {
  if (!params.final_analysis) {
    Channel
        .fromList(params.datasets.entrySet())
        .map { entry ->
            def meta = [ id: entry.key, path: entry.value ]
            tuple(meta, meta.path)
        }
        .set { ch_datasets }

    def combinedDatasetsExists = file("${params.outdir}/combined_datasets.tsv").exists()
    println "combinedDatasetsExists = ${combinedDatasetsExists}"
    if ( combinedDatasetsExists ) {
        EMBEDDING_PARAMETER_DECISION(Channel.fromPath("${params.outdir}/combined_datasets.tsv"))
    }
    else {
        // re-run FIRST_RUN for all datasets, then pass FIRST_RUN emitted channels to SECOND_RUN
        def firstRes = ANNOTATE_PROTEINS(ch_datasets)

        EMBEDDING_PARAMETER_DECISION(firstRes.ch_combined_tsv)
    }
  }
  if (params.final_analysis) {
        REP_MODULE_ANALYSIS("${params.outdir}/combined_datasets.tsv")
  }
}