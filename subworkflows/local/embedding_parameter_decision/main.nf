include { EMBEDDING_PLAN } from '../../../modules/local/embedding_plan'
include { EMBEDDINGS } from '../../../modules/local/embeddings'
include { GENERATE_COORDINATES } from '../../../modules/local/generate_coordinates'
include { PLOT_UMAP } from '../../../modules/local/plot_umap'
include { TILE_UMAPS } from '../../../modules/local/tile_umaps'
include { HDBSCAN } from '../../../modules/local/hdbscan'

workflow EMBEDDING_PARAMETER_DECISION {
    take:
        ch_combined_tsv
        ch_cleaned_fasta
        
    main:
    ch_filtered_tsv = ch_combined_tsv
    ch_metadata = channel.fromPath(params.genofeature_metadata)

    EMBEDDING_PLAN(
        ch_combined_tsv,
        ch_cleaned_fasta.map { meta, fasta -> fasta }.collect()
    )

    ch_embeddings = EMBEDDINGS(
        EMBEDDING_PLAN.out.planned_fastas.flatten().map { fasta ->
            def genofeature = fasta.parent.name
            def protein = fasta.parent.parent.name
            def dataset = fasta.parent.parent.parent.name
            tuple(dataset, protein, genofeature, fasta)
        },
        "${baseDir}/tools/"
    )

    ch_expected_keys = EMBEDDING_PLAN.out.expected_dirs
        .splitText()
        .map { it.trim() }
        .filter { it }
        .collect()
        .map { it as Set }

    ch_existing_embedding_dirs = channel
        .fromPath("${params.outdir}/embeddings/*/*/*", type: 'dir', checkIfExists: false)
        .map { dir -> tuple("${dir.parent.parent.name}/${dir.parent.name}/${dir.name}", dir) }
        .combine(ch_expected_keys)
        .filter { key, dir, expectedSet -> expectedSet.contains(key) }
        .map { key, dir, expectedSet -> dir }

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

    ch_coord_files = GENERATE_COORDINATES.out.coordinates_files
        .map { dir ->
            def m = (dir.name =~ /nn(\d+)_md(\d+)/)
            if (!m) {
                error "Could not parse nn/md from coordinates dir name: ${dir.name}"
            }
            def nn = m[0][1] as Integer
            def md_tag = m[0][2]

            def coordFile = dir.listFiles().find { it.name ==~ /coordinates_nn\d+_md\d+\.tsv/ }
            if (!coordFile) {
                error "No coordinates_nn*_md*.tsv file found in ${dir}"
            }

            def connsFile = dir.resolve('connections.tsv')
            if (!connsFile.exists()) {
                error "Expected connections.tsv in ${dir} but it was missing"
            }

            tuple(nn, md_tag, coordFile, connsFile)
        }

    ch_plot_umap = PLOT_UMAP(
        ch_coord_files,
        ch_filtered_tsv.first(),  
        ch_metadata.first()    
    )

    ch_umap = TILE_UMAPS(
        ch_plot_umap.plot.collect(),
        ch_metadata
    )

    ch_hbd = HDBSCAN(
        ch_embedding_dirs.collect(),
        ch_coordinate_dirs,
        ch_filtered_tsv,
        ch_metadata,
        ch_plot_umap.plot.collect()
    )

    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_umap.tiled_image)
    ch_multiqc_files = ch_multiqc_files.mix(ch_hbd.tiled_image)
    ch_multiqc_files = ch_multiqc_files.mix(ch_hbd.plots)

    ch_versions = channel.empty()
    ch_versions = ch_versions.mix(EMBEDDING_PLAN.out.versions)
    ch_versions = ch_versions.mix(EMBEDDINGS.out.versions)
    ch_versions = ch_versions.mix(HDBSCAN.out.versions)
    ch_versions = ch_versions.mix(TILE_UMAPS.out.versions)
    ch_versions = ch_versions.mix(GENERATE_COORDINATES.out.versions)

    emit:
        ch_combined_tsv = ch_combined_tsv
        ch_versions = ch_versions
        ch_multiqc_files = ch_multiqc_files
}