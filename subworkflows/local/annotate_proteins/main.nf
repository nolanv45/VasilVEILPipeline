include { CLEAN_FASTA_HEADERS } from '../../../modules/local/clean_fasta_headers'
include { PHIDRA } from '../../../modules/local/phidra'
include { CRITERIA_TSV } from '../../../modules/local/criteria_tsv'
include { TOP_HIT_ANNOTATION } from '../../../modules/local/top_hit_annotation'
include { PASV } from '../../../modules/local/pasv'
include { PASV_ANNOTATION } from '../../../modules/local/pasv_annotation'
include { DOMAIN_MATCH_ANNOTATION } from '../../../modules/local/domain_match_annotation'
include { PHIDRA_ONLY } from '../../../modules/local/phidra_only'
include { MERGE_TSV } from '../../../modules/local/merge_tsv'
include { APPLY_CRITERIA } from '../../../modules/local/apply_criteria'
include { ANALYZE_AND_PLOT } from '../../../modules/local/analyze_and_plot'
include { PASV_POST } from '../../../modules/local/pasv_post'
include { COMBINE_DATASETS } from '../../../modules/local/combine_datasets'
include { DUPLICATE_HANDLE } from '../../../modules/local/duplicate_handle'

workflow ANNOTATE_PROTEINS {
    take:
        ch_dataset

    main:
        ch_metadata = channel.fromPath(params.genofeature_metadata)

        // Normalize fasta input file
        CLEAN_FASTA_HEADERS(ch_dataset)
         
        // Prepare PHIDRA input
        ch_dataset_proteins = CLEAN_FASTA_HEADERS.out.fasta
            .combine(channel.fromList(params.proteins))
            .map { meta, fasta, protein_config -> 
                def new_meta = meta + [
                    protein: protein_config.protein,
                    pfamDomain: protein_config.pfamDomain,
                    subjectDB: protein_config.subjectDB,
                    pasv_align_refs: protein_config.containsKey('pasv_align_refs') ? 
                        protein_config.pasv_align_refs : null,
                    cleaned_fasta: fasta
                ]
                // stage the subject DB file so Nextflow copies it into the task workdir
                tuple(new_meta, fasta, file(protein_config.subjectDB), file(protein_config.pfamDomain))
            }
        PHIDRA(ch_dataset_proteins)

        // ----------------------------------------------------------------------
        // CRITERIA TSV
        // ----------------------------------------------------------------------
        ch_phidra_crit = PHIDRA.out.results

        ch_unval_crit = ch_phidra_crit.map { meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs -> [meta.id, unval_pfam] }
        ch_val_crit   = ch_phidra_crit.map { meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs -> [meta.id, val_pfam] }

        ch_criteria_input = ch_unval_crit
            .mix(ch_val_crit)
            .groupTuple(by: 0)
            .map { id, files -> [ [id: id], files.flatten() ] }

        CRITERIA_TSV(ch_criteria_input)

        // ----------------------------------------------------------------------
        // JOIN CRITERIA TSV WITH PHIDRA RESULTS
        // ----------------------------------------------------------------------
        ch_tsv = CRITERIA_TSV.out.results.map { meta, tsv -> [meta.id, tsv] }

        ch_phidra_keyed = PHIDRA.out.results
            .map { meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                [meta.id, meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs]
            }
 
        ch_phidra_with_tsv = PHIDRA.out.results
            .map { meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                [meta.id, meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs]
            }
            .combine(ch_tsv, by: 0)
            .map { id, meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs, tsv ->
                [meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs]
            }

        ch_branched = ch_phidra_with_tsv
            .branch { meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                tophit: meta.protein.toLowerCase() in params.tophit_proteins.collect { it.toLowerCase() }
                pasv:   meta.protein.toLowerCase() in params.pasv_proteins.collect { it.toLowerCase() }
                domain: meta.protein.toLowerCase() in params.domain_proteins.collect { it.toLowerCase() }
                other: true
            }


        ch_branched.other
            .map { meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                [meta, tsv]
            }
            | PHIDRA_ONLY

        ch_branched.tophit
            .map { meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                [meta, tsv, init, recurs]
            }
            | TOP_HIT_ANNOTATION

        ch_pasv_input = ch_branched.pasv
            .map { meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                def settings      = (params.pasv_settings ?: [:]).get(meta.protein, [:])
                def mapped_name   = settings.mapped_name   ?: "${meta.protein}_putative"
                def roi_start     = settings.roi_start     ?: ''
                def roi_end       = settings.roi_end       ?: ''
                def cat_sites     = settings.cat_sites     ?: ''
                def expected_sigs = settings.expected_sigs ?: ''
                def new_meta = meta + [
                    mapped_name:   mapped_name,
                    roi_start:     roi_start,
                    roi_end:       roi_end,
                    cat_sites:     cat_sites,
                    expected_sigs: expected_sigs
                ]
                tuple(new_meta, val_fasta, unval_fasta, file(meta.pasv_align_refs))
            }

        ch_pasv_input | PASV

        ch_pasv_tsv = ch_branched.pasv
            .map { meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                ["${meta.id}:${meta.protein}", tsv]
            }

        PASV.out.results
            .map { meta, pasv_sig -> ["${meta.id}:${meta.protein}", meta, pasv_sig] }
            .join(ch_pasv_tsv, by: 0)
            .map { key, meta, pasv_sig, tsv -> [meta, tsv, pasv_sig] }
            | PASV_ANNOTATION

        ch_branched.domain
            .map { meta, tsv, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                [meta, tsv]
            }
            | DOMAIN_MATCH_ANNOTATION



        // ----------------------------------------------------------------------
        // MERGE ALL PER-PROTEIN TSVs BACK INTO ONE PER ID
        // ----------------------------------------------------------------------
        TOP_HIT_ANNOTATION.out.results
            .mix(PASV_ANNOTATION.out.results, DOMAIN_MATCH_ANNOTATION.out.results, PHIDRA_ONLY.out.results)
            .map { meta, tsv -> [meta.id, tsv] }
            .groupTuple(by: 0)
            .map { id, tsvs -> [ [id: id], tsvs.flatten() ] }
            | MERGE_TSV


        // APPLIES RULESET BELOW
        MERGE_TSV.out.results | APPLY_CRITERIA

        ch_cleaned_fasta = CLEAN_FASTA_HEADERS.out.fasta
            .map { meta, fasta -> [meta.id, fasta] }

        ch_analyze_input = PHIDRA.out.results
            .filter { meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs ->
                !(meta.protein.toLowerCase() in params.pasv_proteins.collect { it.toLowerCase() })
            }
            .map { meta, val_fasta, val_pfam, unval_fasta, unval_pfam, init, recurs -> [meta.id, meta] }
            .combine(
                APPLY_CRITERIA.out.results.map { meta, tsv -> [meta.id, tsv] },
                by: 0
            )
            .combine(ch_cleaned_fasta, by: 0)
            .combine(ch_metadata)
            .map { id, meta, tsv, cleaned_fasta, metadata_file ->
                tuple(meta, tsv, metadata_file, cleaned_fasta)
            }

        ch_analyze_input
            | ANALYZE_AND_PLOT



        PASV.out.results
            .map { meta, pasv_sig -> [meta.id, meta, pasv_sig] }
            .combine(
                APPLY_CRITERIA.out.results.map { meta, tsv -> [meta.id, tsv] },
                by: 0
            )
            .combine(ch_cleaned_fasta, by: 0)
            .combine(ch_metadata)
            .map { id, meta, pasv_sig, filtered_tsv, cleaned_fasta, metadata_file ->
                tuple(meta, pasv_sig, metadata_file, cleaned_fasta, filtered_tsv)
            }
            | PASV_POST



        ch_combine_datasets = APPLY_CRITERIA.out.results
            .map { meta, file -> file.toAbsolutePath() }
            .collect()

        COMBINE_DATASETS(ch_combine_datasets)

        DUPLICATE_HANDLE(APPLY_CRITERIA.out.results)

        ch_annotation_plots = ANALYZE_AND_PLOT.out.multiqc_plot
            .mix(PASV_POST.out.multiqc_plot)

        ch_versions = PHIDRA.out.versions
            .mix(CRITERIA_TSV.out.versions)
            .mix(TOP_HIT_ANNOTATION.out.versions)
            .mix(PASV_ANNOTATION.out.versions)
            .mix(DOMAIN_MATCH_ANNOTATION.out.versions)
            .mix(PHIDRA_ONLY.out.versions)
            .mix(MERGE_TSV.out.versions)
            .mix(APPLY_CRITERIA.out.versions)
            .mix(ANALYZE_AND_PLOT.out.versions)
            .mix(PASV.out.versions)
            .mix(PASV_POST.out.versions)
            .mix(COMBINE_DATASETS.out.versions)


    emit:
        ch_combined_tsv = COMBINE_DATASETS.out.results
        multiqc_files = ch_annotation_plots
        versions = ch_versions
}



