include { EMBEDDING_PLAN } from '../../../modules/local/embedding_plan'
include { EMBEDDINGS } from '../../../modules/local/embeddings'
include { GENERATE_COORDINATES } from '../../../modules/local/generate_coordinates'
include { UMAP_PROJECTION } from '../../../modules/local/umap_projection'
include { HDBSCAN } from '../../../modules/local/hdbscan'

workflow EMBEDDING_PARAMETER_DECISION {
    take:
        ch_combined_tsv
        
    main:
    ch_filtered_tsv = ch_combined_tsv
    ch_metadata = channel.fromPath(params.genofeature_metadata)

    EMBEDDING_PLAN(ch_combined_tsv)

    ch_embeddings = EMBEDDINGS(
        EMBEDDING_PLAN.out.planned_fastas.flatten().map { fasta ->
            def genofeature = fasta.parent.name
            def protein = fasta.parent.parent.name
            def dataset = fasta.parent.parent.parent.name
            tuple(dataset, protein, genofeature, fasta)
        },
        "${baseDir}/tools/"
    )

    ch_existing_embedding_dirs = channel
        .fromPath("${params.outdir}/embeddings/*", type: 'dir')

    ch_embedding_dirs = ch_existing_embedding_dirs
        .mix(ch_embeddings.embeddings_dirs)
        .unique()

    def excluded_genofeatures = null

    // create all (nn, md) combinations
    def parameter_combinations = []

    params.nn.each { nn ->
        params.md.each { md ->
            parameter_combinations << tuple(nn, md)
        }
    }

    // Create one input tuple for each parameter combination
    ch_generate_coordinates_inputs = ch_embedding_dirs
        .collect()
        .flatMap { dirs ->

            parameter_combinations.collect { combo ->

                def nn = combo[0]
                def md = combo[1]

                def md_tag = (md * 10).intValue()

                tuple(
                    dirs,
                    excluded_genofeatures,
                    nn,
                    md,
                    "04_parameter_selection/coordinates/nn${nn}_md${md_tag}"
                )
            }
        }

    ch_generate_coordinates_inputs.view()

    GENERATE_COORDINATES(ch_generate_coordinates_inputs)

    ch_coordinate_dirs = GENERATE_COORDINATES.out.coordinates_files.collect()

    ch_umap = UMAP_PROJECTION(
        ch_coordinate_dirs,
        ch_filtered_tsv,
        ch_metadata
    )

    ch_hbd = HDBSCAN(
        ch_embedding_dirs.collect(),
        ch_coordinate_dirs,
        ch_filtered_tsv,
        ch_metadata,
        ch_umap.plots       
    )

    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_umap.tiled_image)
    ch_multiqc_files = ch_multiqc_files.mix(ch_hbd.tiled_image)
    ch_multiqc_files = ch_multiqc_files.mix(ch_hbd.plots)

    ch_versions = channel.empty()
    ch_versions = ch_versions.mix(EMBEDDING_PLAN.out.versions)
    ch_versions = ch_versions.mix(EMBEDDINGS.out.versions)
    ch_versions = ch_versions.mix(HDBSCAN.out.versions)
    ch_versions = ch_versions.mix(UMAP_PROJECTION.out.versions)
    ch_versions = ch_versions.mix(GENERATE_COORDINATES.out.versions)

    emit:
        ch_combined_tsv = ch_combined_tsv
        versions = ch_versions
        multiqc_files = ch_multiqc_files
}