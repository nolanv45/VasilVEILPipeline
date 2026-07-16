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

    channel
        .from(parameter_combinations)
        .set { ch_parameter_combinations }

    // Create one input tuple for each parameter combination
    ch_embedding_dirs
        .collect()
        .combine(ch_parameter_combinations)
        .map { embedding_dirs, param_tuple ->
            def (nn, md) = param_tuple
            tuple(
                embedding_dirs,
                excluded_genofeatures,
                nn,
                md,
                "04_parameter_selection/coordinates"
            )
        }
        .set { ch_generate_coordinates_inputs }

    GENERATE_COORDINATES(ch_generate_coordinates_inputs)

    ch_umap = UMAP_PROJECTION(
        GENERATE_COORDINATES.out.coordinates_files,
        ch_filtered_tsv,
        ch_metadata
    )

    ch_hbd = HDBSCAN(
        ch_embedding_dirs.collect(),
        GENERATE_COORDINATES.out.coordinates_files,
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