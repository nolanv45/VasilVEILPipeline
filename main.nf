nextflow.enable.dsl=2

process PHIDRA {
    conda "/home/nolanv/.conda/envs/phidra"

    input:
    val protein_config


    output:
    path "phidra_output_test/", emit: "phidra_output"

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
    input:
    path phidra_output_test

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
        [helicase]="HEL_putative"
        [RecA]="HEL_putative"
        [UvsX]="HEL_putative"
        [RecB]="HEL_putative"
        [RecD]="HEL_putative"
        [PcrA]="HEL_putative"
        [UvrD]="HEL_putative"
        [Dda]="HEL_putative"
        [UvsW]="HEL_putative"
        [UvrB]="HEL_putative"
        [RecG]="HEL_putative"
        [SNF2]="HEL_putative"
        [Gp4]="HEL_putative"
        [Gp41]="HEL_putative"
        [DnaB]="HEL_putative"
        [RNR]="RNR_putative"
    )

    for dir in ${phidra_output_test}/*; do 
        protein=\$(basename "\$dir")

        if [[ -z "\${MAP[\$protein]}" ]]; then
            echo "Warning: Protein '\$protein' not recognized in mapping. Skipping."
            continue
        fi
        src="\$dir/final_results/pfam_validated_full_protein.fa"
        dest="\$WORK_DIR/pasv_output/input/\${MAP[\$protein]}.fasta"

        if [[ -f "\$src" ]]; then
            cp "\$src" "\$dest"
        else
            echo "Warning: Source file '\$src' does not exist. Skipping."
        fi
    done

    /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/pasv.sh \$WORK_DIR/pasv_output
    """

}


workflow {
    Channel.fromList(params.proteins)
        .set { protein_configs }
    phidra_results = protein_configs | PHIDRA

    phidra_results | PASV
}
