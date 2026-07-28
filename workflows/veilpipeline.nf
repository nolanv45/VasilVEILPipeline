nextflow.enable.dsl=2

// SUBWORKFLOW: Three subworkflows consisting of local modules
include { REP_MODULE_ANALYSIS } from '../subworkflows/local/rep_module_analysis'
include { ANNOTATE_PROTEINS } from '../subworkflows/local/annotate_proteins'
include { EMBEDDING_PARAMETER_DECISION } from '../subworkflows/local/embedding_parameter_decision'
include { MULTIQC as MULTIQC_SWEEP } from '../modules/nf-core/multiqc'
include { MULTIQC as MULTIQC_FINAL } from '../modules/nf-core/multiqc'
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
    ch_multiqc_files_sweep = channel.empty()
    ch_multiqc_files_final = channel.empty()

    // shared prep items, built once, mixed 
    ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
        .collectFile(name: 'workflow_summary_mqc.yaml')

    ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        .collectFile(name: 'methods_description_mqc.yaml', sort: true)

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
            ch_versions = ch_versions.mix(embedRes.ch_versions)
            ch_multiqc_files_sweep = ch_multiqc_files_sweep.mix(embedRes.ch_multiqc_files)
            
            // ANNOTATE_PROTEINS didn't run this time (combined_datasets.tsv already existed),
            // so pull in whatever it published from a prior run to keep those plots in the report
            ch_existing_annotation_multiqc = channel
                .fromPath("${params.outdir}/03_annotation_analysis/**/*_mqc.png", checkIfExists: false)
            ch_multiqc_files_sweep = ch_multiqc_files_sweep.mix(ch_existing_annotation_multiqc)
        }
        else {
            def annotationRes = ANNOTATE_PROTEINS(ch_datasets)
            ch_versions = ch_versions.mix(annotationRes.ch_versions)
            ch_multiqc_files_sweep = ch_multiqc_files_sweep.mix(annotationRes.ch_multiqc_files)

            def embedRes = EMBEDDING_PARAMETER_DECISION(annotationRes.ch_combined_tsv)
            ch_versions = ch_versions.mix(embedRes.ch_versions)
            ch_multiqc_files_sweep = ch_multiqc_files_sweep.mix(embedRes.ch_multiqc_files)
        }
    }
    if (params.final_analysis) {
        def repRes = REP_MODULE_ANALYSIS("${params.outdir}/combined_datasets.tsv")
        ch_versions = ch_versions.mix(repRes.ch_versions)
        ch_multiqc_files_final = ch_multiqc_files_final.mix(repRes.ch_multiqc_files)
    }
    // software versions
    ch_software_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(
            name: "${workflow.manifest.name}_software_mqc_versions.yml",
            storeDir: "${params.outdir}/pipeline_info",
            newLine: true
        )

    if (!params.final_analysis) {
        ch_multiqc_files_sweep = ch_multiqc_files_sweep
            .mix(ch_versions)
            .mix(ch_workflow_summary)
            .mix(ch_methods_description)
            .mix(ch_software_versions)

        ch_multiqc_input_sweep = ch_multiqc_files_sweep.flatten().collect().map { files ->
            tuple(
                [ id: "${workflow.manifest.name}_parameter_sweep" ],
                files,
                multiqc_config ? file(multiqc_config, checkIfExists: true) : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                []
            )
        }

        MULTIQC_SWEEP(ch_multiqc_input_sweep)
        ch_multiqc_report_sweep = MULTIQC_SWEEP.out.report.map { _meta, report -> report }
    }
    else {
        ch_multiqc_report_sweep = channel.empty()
    }

    if (params.final_analysis) {
        ch_multiqc_files_final = ch_multiqc_files_final
            .mix(ch_versions)
            .mix(ch_workflow_summary)
            .mix(ch_methods_description)
            .mix(ch_software_versions)

        ch_multiqc_input_final = ch_multiqc_files_final.flatten().collect().map { files ->
            tuple(
                [ id: "${workflow.manifest.name}_final_analysis" ],
                files,
                multiqc_config ? file(multiqc_config, checkIfExists: true) : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                []
            )
        }

        MULTIQC_FINAL(ch_multiqc_input_final)
        ch_multiqc_report_final = MULTIQC_FINAL.out.report.map { _meta, report -> report }
    }
    else {
        ch_multiqc_report_final = channel.empty()
    }

    emit:
    multiqc_report_sweep = ch_multiqc_report_sweep
    multiqc_report_final = ch_multiqc_report_final
    versions             = ch_software_versions
}