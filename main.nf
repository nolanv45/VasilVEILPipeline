nextflow.enable.dsl=2

include { EMBEDDINGS; HDBSCAN; GENERATE_COORDINATES; UMAP_PROJECTION } from './second.nf'
include {CLEAN_FASTA_HEADERS; 
PHIDRA;
ANNOTATE_HITS;
PASV;
PASV_POST;
PHIDRA_ONLY_SUMMARY;
DOMAIN_MATCH;
COMBINE_PHIDRA_TSV;
ANALYZE_AND_PLOT;
STANDARDIZE_OUTPUTS;
COMBINE_DATASETS;
DUPLICATE_HANDLE;
SPLIT_BY_GENOFEATURE;} from './first.nf'

include {MODULE_FILE; GENERATE_COORDINATES_2; GENOFEATURE_CENTRIC; MODIFY_CLUSTERS; ZEROFIVEC} from './third.nf'




workflow {
    Channel
        .fromList(params.datasets.entrySet())
        .map { entry ->
            def meta = [ id: entry.key, path: entry.value ]
            tuple(meta, meta.path)
        }
        .set { ch_datasets }

    def anyMissing = params.datasets.any { id, path ->
        def emb_dir = file("${params.outdir}/${id}/files_for_embeddings")
        !emb_dir.exists()
    }
    println "anyMissing = ${anyMissing}"
    if ( anyMissing ) {
        // re-run FIRST_RUN for all datasets, then pass FIRST_RUN emitted channels to SECOND_RUN
        def firstRes = FIRST_RUN(ch_datasets)

        SECOND_RUN(firstRes.ch_combined_tsv)
    }
    else {
        SECOND_RUN("${params.outdir}/combined_datasets.tsv")
    }
}

workflow FIRST_RUN {
    take:
        ch_dataset

    main:
        CLEAN_FASTA_HEADERS(ch_dataset)
         
        ch_dataset_proteins = CLEAN_FASTA_HEADERS.out.fasta
            .combine(Channel.fromList(params.proteins))
            .map { meta, fasta, protein_config -> 
                def new_meta = meta + [
                    protein: protein_config.protein,
                    pfamDomain: protein_config.pfamDomain,
                    subjectDB: protein_config.subjectDB,
                    pasv_align_refs: protein_config.containsKey('pasv_align_refs') ? 
                        protein_config.pasv_align_refs : null,
                    cleaned_fasta: fasta
                ]
                tuple(new_meta, fasta)
            }

        PHIDRA(ch_dataset_proteins)

        ch_branched = PHIDRA.out.results
            .branch { meta, fasta, init_search, pfam -> 
                annotation: meta.protein.toLowerCase() in (params.tophit_proteins.collect { it.toLowerCase() })
                pasv: meta.protein.toLowerCase() in (params.pasv_proteins.collect { it.toLowerCase() })
                pfam: meta.protein.toLowerCase() in (params.domain_proteins.collect { it.toLowerCase() })
                phidra_only: true
            }

        ch_pasv = ch_branched.pasv
            .map { meta, fasta, init_search, pfam ->
                def settings = (params.pasv_settings ?: [:]).get(meta.protein, [:])
                def mapped_name = settings.mapped_name ?: "${meta.protein}_putative"
                def roi_start   = settings.roi_start   ?: ''
                def roi_end     = settings.roi_end     ?: ''
                def cat_sites   = settings.cat_sites   ?: ''
                def expected_sigs  = settings.expected_sigs   ?: ''
    
                def new_meta = meta + [
                    mapped_name: mapped_name,
                    roi_start: roi_start,
                    roi_end: roi_end,
                    cat_sites: cat_sites,
                    expected_sigs: expected_sigs
                ]
                tuple(new_meta, fasta, meta.pasv_align_refs)
            }

        ANNOTATE_HITS(ch_branched.annotation)
        PASV(ch_pasv)

        PASV_POST(PASV.out.results)
        PHIDRA_ONLY_SUMMARY(ch_branched.phidra_only)
        DOMAIN_MATCH(ch_branched.pfam
            .map { meta, fasta, init_search, pfam -> 
                tuple(meta, fasta, pfam, params.pfam_annotation_map)
            }
        )

        ch_annotate = ANNOTATE_HITS.out.results ?: Channel.empty()
        ch_phidra_only = PHIDRA_ONLY_SUMMARY.out.results ?: Channel.empty()
        
        // Combine channels
        ch_phidra_results = ch_annotate
            .mix(ch_phidra_only)
            .map { meta, file -> 
                tuple(meta.id, meta, file)  // Add grouping key
            }
            .groupTuple(by: 0)  // Group by dataset ID
            .map { id, metas, files ->
                tuple(metas[0], files.flatten())  // Flatten the files array
            }

        COMBINE_PHIDRA_TSV(ch_phidra_results)

        ANALYZE_AND_PLOT(COMBINE_PHIDRA_TSV.out.combined)

        ch_pasv_map = PASV_POST.out
            .map { meta, processed_file, stats_file, plot_file ->
                tuple(meta, processed_file)
            }

        ch_domain_annotate_map = DOMAIN_MATCH.out.results
            .map { meta, processed_file ->
                tuple(meta, processed_file)
            }

        ch_combined_phidra = COMBINE_PHIDRA_TSV.out.combined ?: Channel.empty()
        ch_pasv_map = ch_pasv_map ?: Channel.empty()

        ch_files_to_standardize = ch_combined_phidra
            .mix(ch_pasv_map, ch_domain_annotate_map)
            .map { meta, file -> 
                tuple(meta.id, meta, file)
            }
            .groupTuple(by: 0)
            .map { id, metas, files ->
                def meta0 = metas[0]
                def files_list = files.flatten().collect { it.toString() }
                tuple(meta0, files_list)
            }

        STANDARDIZE_OUTPUTS(ch_files_to_standardize)

        ch_combine_datasets = STANDARDIZE_OUTPUTS.out.standardized
            .map { meta, file ->
                file.toAbsolutePath()
            }
            .collect()
        
        COMBINE_DATASETS(ch_combine_datasets)

        DUPLICATE_HANDLE(STANDARDIZE_OUTPUTS.out)

        ch_split_by_genofeature = SPLIT_BY_GENOFEATURE(
            STANDARDIZE_OUTPUTS.out
        )
        ch_combined_tsv = COMBINE_DATASETS.out.combined

    emit:
        ch_combined_tsv
}

workflow SECOND_RUN {
    take:
        ch_combined_tsv
        
    main:
    ch_filtered_tsv = ch_combined_tsv
    ch_metadata = Channel.fromPath(params.genofeature_metadata)

    ch_embedding_datasets = Channel.value(params.embedding_datasets)

    ch_embeddings = EMBEDDINGS(
        ch_embedding_datasets, ch_combined_tsv
    )
    // embeddings found, generate new embeddings for these bc didnt find in output file.

    ch_coordinates = GENERATE_COORDINATES(
        ch_embeddings
        // "/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/full_ena_output/embeddings"
    )

    // remember to fix the input from coordinates the same way you did to pasv output.
    ch_umap = UMAP_PROJECTION(
        // "/mnt/VEIL/users//nolanv/pipeline_project/VasilVEILPipeline/figures_folder/coordinates",
        ch_coordinates,
        ch_filtered_tsv,
        ch_metadata
    )

    ch_hbd = HDBSCAN(
        // "/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/figures_folder/coordinates",
        ch_coordinates,
        ch_filtered_tsv,
        ch_metadata,
        ch_umap.plots
        // "/mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/work/0b/eab48116d13496090e556b588f6dfa/plots"
        
    )
    emit:
        ch_combined_tsv
}

// workflow THIRD_RUN {
//     take:
//         ch_combined_tsv

//     main:
//     ch_filtered_tsv = Channel.fromPath(params.second_run)
//     ch_metadata = Channel.fromPath(params.genofeature_metadata)

//     ch_module_file = MODULE_FILE(
//         ch_filtered_tsv,
//         ch_metadata
//     )

//     GENERATE_COORDINATES_2(params.embeddings)
//     ch_coordinates = GENERATE_COORDINATES_2.out.coordinates_tsv
//     ch_connections = GENERATE_COORDINATES_2.out.connections_tsv
//     MODIFY_CLUSTERS(ch_cluster_dir)
//     // ZEROFIVEC(GENERATE_COORDINATES_2.coordinates_tsv, MODIFY_CLUSTERS.out)

//     GENOFEATURE_CENTRIC(ch_module_file, ch_coordinates, ch_connections)

// }