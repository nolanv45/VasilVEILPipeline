nextflow.enable.dsl=2

// SUBWORKFLOW: Three subworkflows consisting of local modules
include { REP_MODULE_ANALYSIS } from '../subworkflows/rep_module_analysis'
include { ANNOTATE_PROTEINS } from '../subworkflows/annotate_proteins'
include { EMBEDDING_PARAMETER_DECISION } from '../subworkflows/embedding_parameter_decision'
include { MULTIQC } from '../modules/nf-core/multiqc'
include { paramsSummaryMap } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_veilpipeline_pipeline'

workflow VEILPIPELINE {
    take:
    multiqc_config
    multiqc_logo
    multiqc_methods_description



    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    if (!params.final_analysis) {
        channel
            .fromList(params.datasets.entrySet())
            .map { entry ->
                def meta = [ id: entry.key, path: entry.value ]
                tuple(meta, meta.path)
            }
            .set { ch_datasets }

        def combinedDatasetsExists = file("${params.outdir}/combined_datasets.tsv").exists()
        println "combinedDatasetsExists = ${combinedDatasetsExists}"
        if ( combinedDatasetsExists ) {
            def embedRes = EMBEDDING_PARAMETER_DECISION(channel.fromPath("${params.outdir}/combined_datasets.tsv"))
            ch_versions = ch_versions.mix(embedRes.versions)
            ch_multiqc_files = ch_multiqc_files.mix(embedRes.multiqc_files)
        }
        else {
            def annotationRes = ANNOTATE_PROTEINS(ch_datasets)
            ch_versions = ch_versions.mix(annotationRes.versions)
            ch_multiqc_files = ch_multiqc_files.mix(annotationRes.multiqc_files)

            def embedRes = EMBEDDING_PARAMETER_DECISION(annotationRes.ch_combined_tsv)
            ch_versions = ch_versions.mix(embedRes.versions)
            ch_multiqc_files = ch_multiqc_files.mix(embedRes.multiqc_files)
        }
    }
    if (params.final_analysis) {
            def repRes = REP_MODULE_ANALYSIS("${params.outdir}/combined_datasets.tsv")
            ch_versions = ch_versions.mix(repRes.versions)
     }

    ch_multiqc_files = ch_multiqc_files.mix(ch_versions)

    ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    ch_software_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(
            name: "${workflow.manifest.name}_software_mqc_versions.yml",
            storeDir: "${params.outdir}/pipeline_info",
            newLine: true
        )

    ch_multiqc_files = ch_multiqc_files.mix(ch_software_versions)

    ch_multiqc_input = ch_multiqc_files.flatten().collect().map { files ->
        tuple(
            [ id: workflow.manifest.name ],
            files,
            multiqc_config ? file(multiqc_config, checkIfExists: true) : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
            multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
            [],
            []
        )
    }

    MULTIQC(ch_multiqc_input)

    ch_multiqc_report = MULTIQC.out.report.map { _meta, report -> report }

    emit:
    multiqc_report = ch_multiqc_report
    versions = ch_software_versions
}