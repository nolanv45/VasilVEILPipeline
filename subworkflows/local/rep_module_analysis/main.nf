include { GENERATE_COORDINATES } from '../../../modules/local/generate_coordinates'
include { MODIFY_CLUSTERS } from '../../../modules/local/modify_clusters'
include { HDBSCAN_TSV } from '../../../modules/local/hdbscan_tsv'
include { HDBSCAN_VISUALS } from '../../../modules/local/hdbscan_visuals'
include { PLOT_GENOFEATURE_CENTRIC } from '../../../modules/local/plot_genofeature_centric'
include { TILE_GENOFEATURE_CENTRIC } from '../../../modules/local/tile_genofeature_centric'
include { MODULE_FILE } from '../../../modules/local/module_file'
workflow REP_MODULE_ANALYSIS {
    take:
        ch_combined_tsv

    main:
    ch_metadata = channel.fromPath(params.genofeature_metadata)

    MODULE_FILE(
        ch_combined_tsv,
        ch_metadata
    )
    ch_module_file = MODULE_FILE.out.module_file

    def excluded_genofeatures = params.excluded_genofeatures.collect { feature -> "\'${feature}\'" }.join(', ')
    def final_nn = "${params.final_nn}"
    def final_md = "${params.final_md}"
    def final_md_tag = String.format('%.1f', params.final_md as Double).replace('.', 'p')
    ch_generate_coordinates_inputs = channel.of(
        tuple(["${params.outdir}/embeddings"], excluded_genofeatures, final_nn, final_md, "05_final_analysis/coordinates/nn${final_nn}_md${final_md_tag}")
    )




    GENERATE_COORDINATES(ch_generate_coordinates_inputs)
    ch_coordinates_dir = GENERATE_COORDINATES.out.coordinates_files
    ch_coordinates = ch_coordinates_dir.map { dir ->
        file(dir).listFiles().find { it.name.startsWith("coordinates_nn") }
    }
    ch_connections = ch_coordinates_dir.map { dir ->
        file("${dir}/connections.tsv")
    }

    ch_coordinates.view { "coordinates: $it" }
    ch_connections.view { "connections: $it" }
    MODIFY_CLUSTERS("${params.outdir}/04_parameter_selection/hdbscan/clusters_csv")
    
    HDBSCAN_TSV(MODIFY_CLUSTERS.out.modified_clusters, ch_combined_tsv)
    HDBSCAN_VISUALS(HDBSCAN_TSV.out.cluster_info, ch_metadata)

    ch_features = ch_module_file
        .map { f ->
            def header = f.readLines()[0].split('\t')
            header.findAll { col -> !(col in ['dataset', 'contig_id', 'module']) }
        }
        .flatMap { it }

    ch_genofeature_centric_plots = PLOT_GENOFEATURE_CENTRIC(
        ch_features,
        ch_module_file.first(),
        ch_coordinates.first(),
        ch_connections.first(),
        params.genofeature_metadata
    )
    ch_genofeature_centric = TILE_GENOFEATURE_CENTRIC(
        ch_genofeature_centric_plots.overview.mix(ch_genofeature_centric_plots.by_dataset).collect()
    )

    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(HDBSCAN_VISUALS.out.composition_heatmap)
    ch_multiqc_files = ch_multiqc_files.mix(HDBSCAN_VISUALS.out.dataset_stacked)
    ch_multiqc_files = ch_multiqc_files.mix(HDBSCAN_VISUALS.out.genofeature_stacked)
    ch_multiqc_files = ch_multiqc_files.mix(HDBSCAN_VISUALS.out.genofeature_size_scaled)
    ch_multiqc_files = ch_multiqc_files.mix(ch_genofeature_centric.tiled)

    ch_versions = channel.empty()
    ch_versions = ch_versions.mix(MODULE_FILE.out.versions)
    ch_versions = ch_versions.mix(GENERATE_COORDINATES.out.versions)
    ch_versions = ch_versions.mix(MODIFY_CLUSTERS.out.versions)
    ch_versions = ch_versions.mix(HDBSCAN_TSV.out.versions)
    ch_versions = ch_versions.mix(HDBSCAN_VISUALS.out.versions)
    ch_versions = ch_versions.mix(TILE_GENOFEATURE_CENTRIC.out.versions)

    emit:
        ch_versions = ch_versions
        ch_multiqc_files = ch_multiqc_files

}