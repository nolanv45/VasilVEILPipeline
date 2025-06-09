nextflow.enable.dsl=2

process PHIDRA {
    conda "/home/nolanv/.conda/envs/phidra"

    input:
    val protein_config


    output:
    path "phidra_output_test/*", emit: "phidra_output"

    script:
    """
    source /etc/profile.d/conda.sh
    conda activate /home/nolanv/.conda/envs/phidra

    WORK_DIR=\$PWD

    mkdir -p phidra_output_test/

    cd ${params.phidra_dir}

    python phidra_run.py \
        -i ${params.input_fasta} \
        -db ${protein_config.subjectDB} \
        -pfam ${params.pfamDB} \
        -ida ${protein_config.pfamDomain} \
        -f ${protein_config.protein} \
        -o \$WORK_DIR/phidra_output_test/ \
        -t 24
    """
}

// process preparePASVInput {
//     input:
//     path phidra_output

//     output:
//     path "pasv_ready"

//     script:
//     """
//     mkdir -p pasv_ready
//     """
// }

process PASV {
    conda "/home/nolanv/.conda/envs/phidra"

    input:
    path phidra_dirs

    output:
    path "pasv_output/", emit: "pasv_output"

    script:
    """
    # Use current working directory
    WORK_DIR=\$(pwd)

    mkdir -p \$WORK_DIR/pasv_output/input
    mkdir -p \$WORK_DIR/pasv_output/output
    mkdir -p \$WORK_DIR/pasv_output/pasv


    declare -A MAP=(
        [PolA]="POL_putative"
        [RNR]="RNR_putative"
    )

    for dir in ${phidra_dirs}; do
        protein=\$(basename "\$dir")

        fasta_file="\$dir/final_results/pfam_validated_full_protein.fa"

        mapped_name="\${MAP[\$protein]}"

        if [[ -z "\$mapped_name" ]]; then
            echo "Warning: Protein '\$protein' not recognized in mapping. Skipping."
            continue
        fi

        

        if [[ -f "\$fasta_file" ]]; then
            cp "\$fasta_file" "\$WORK_DIR/pasv_output/input/\$mapped_name.fasta"
        else
            echo "Warning: File not found: \$fasta_file"
        fi
    done

    /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/pasv.sh \$WORK_DIR/pasv_output
    """

}

process annotate_top_hits {
    input:
    path initial_search
    path recursive_search
    val output_file

    output:
    path "${output_file}"

    script:
    """
    touch annotated_hits.tsv
    python annotate_phidra_top_hit_type.py \\
        ${initial_search} \\
        ${recursive_search} \\
        ${output_file} \\
    
    """
}

workflow {
    Channel.fromList(params.proteins)
        .set { protein_configs }
    phidra_results = protein_configs | PHIDRA
    phidra_results.collect().set { all_phidra_outputs }

    // Run PASV only if PolA or RNR is involved
    if (params.proteins.any { it.toLowerCase() == "pola" || it.toLowerCase() == "rnr" }) {
        PASV(all_phidra_outputs)
    }



    // Run annotation if helicase or polB is used
    for (protein in params.proteins) {
        def lower = protein.toLowerCase()
        if (lower == "helicase" || lower == "polb") {
            def init_path = file("phidra_output_test/${protein}/mmseqs_results/initial_search/${protein}_TopHit_Evalue.tsv")
            def rec_path = file("phidra_output_test/${protein}/mmseqs_results/initial_search/${protein}_Recursive_TopHit_Evalue.tsv")
            def out_path = "phidra_output_test/${protein}_annotated_hits.tsv"
            annotate_top_hits(init_path, rec_path, out_path)
        }
    }
}