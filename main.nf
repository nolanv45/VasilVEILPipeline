nextflow.enable.dsl=2

// process PHIDRA {
//     conda "/home/nolanv/.conda/envs/phidra"

//     input:
//     val protein_config


//     output:
//     path "*", emit: "phidra_output"

//     script:
//     """
//     source /etc/profile.d/conda.sh
//     conda activate /home/nolanv/.conda/envs/phidra

//     WORK_DIR=\$PWD

//     cd ${params.phidra_dir}

//     python phidra_run.py \
//         -i ${params.input_fasta} \
//         -db ${protein_config.subjectDB} \
//         -pfam ${params.pfamDB} \
//         -ida ${protein_config.pfamDomain} \
//         -f ${protein_config.protein} \
//         -o \$WORK_DIR \
//         -t 24
//     """
// }

process PHIDRA {
    conda "/home/nolanv/.conda/envs/phidra"

    input:
    val protein_config


    output:
    tuple val(protein_config.protein), 
          path("*/final_results/pfam_validated_full_protein.fa"), 
          path("*/mmseqs_results/initial_search/*_TopHit_Evalue.tsv"), 
          path("*/mmseqs_results/recursive_search/*_TopHit_Evalue.tsv")

    script:
    """
    source /etc/profile.d/conda.sh
    conda activate /home/nolanv/.conda/envs/phidra

    WORK_DIR=\$PWD

    cd ${params.phidra_dir}

    python phidra_run.py \
        -i ${params.input_fasta} \
        -db ${protein_config.subjectDB} \
        -pfam ${params.pfamDB} \
        -ida ${protein_config.pfamDomain} \
        -f ${protein_config.protein} \
        -o \$WORK_DIR \
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

// process PASV {
//     conda "/home/nolanv/.conda/envs/phidra"

//     input:
//     path phidra_dirs

//     output:
//     path "pasv_output/", emit: "pasv_output"

//     script:
//     """
//     # Use current working directory
//     WORK_DIR=\$(pwd)

//     mkdir -p \$WORK_DIR/pasv_output/input
//     mkdir -p \$WORK_DIR/pasv_output/output
//     mkdir -p \$WORK_DIR/pasv_output/pasv


//     declare -A MAP=(
//         [PolA]="POL_putative"
//         [RNR]="RNR_putative"
//     )

//     for dir in ${phidra_dirs}; do
//         protein=\$(basename "\$dir")

//         fasta_file="\$dir/final_results/pfam_validated_full_protein.fa"

//         mapped_name="\${MAP[\$protein]}"

//         if [[ -z "\$mapped_name" ]]; then
//             echo "Warning: Protein "\$protein" not recognized in mapping. Skipping."
//             continue
//         fi

        

//         if [[ -f "\$fasta_file" ]]; then
//             cp "\$fasta_file" "\$WORK_DIR/pasv_output/input/\$mapped_name.fasta"
//         else
//             echo "Warning: File not found: \$fasta_file"
//         fi
//     done

//     /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/pasv.sh \$WORK_DIR/pasv_output
//     """

// }

process PASV {
    conda "/home/nolanv/.conda/envs/phidra"

    input:
    tuple val(protein), path(fasta)

    output:
    path "pasv_output/", emit: "pasv_output"

    script:
    """
    # Use current working directory
    WORK_DIR=\$(pwd)

    # Create directories
    mkdir -p \$WORK_DIR/pasv_output/{input,output,pasv}

    # Define mapping
    if [ "${protein}" == "PolA" ]; then
        mapped_name="POL_putative"
    elif [ "${protein}" == "RNR" ]; then
        mapped_name="RNR_putative"
    else
        echo "Error: Unknown protein ${protein}"
        exit 1
    fi

    # Copy input file with mapped name
    cp ${fasta} "\$WORK_DIR/pasv_output/input/\${mapped_name}.fasta"

    # Run PASV
    /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/pasv.sh \$WORK_DIR/pasv_output
    """
}



process PASV_test {
    input:
    path pasv_dir

    output:
    path "pasv_output", emit: "pasv_test_output"

    script:
    """
    /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/pasv.sh ${pasv_dir}
    """
}

process cleanHeaders {
    script:
    """
    echo "cleanHeaders placeholder"
    """
}

process annotate_top_hits {
    debug true

    input:
    tuple val(protein), 
          path(initial_search, stageAs: { it[0] + "_initial.tsv" }), 
          path(recursive_search, stageAs: { it[0] + "_recursive.tsv" })

    output:
    path "annotated_hits_${protein}.tsv"

    script:
    """
    #!/usr/bin/env python3
import sys
import os

protein = "${protein}"
file1_path = f"{protein}_initial.tsv"
file2_path = f"{protein}_recursive.tsv"
output_file = f"annotated_hits_{protein}.tsv"

def parse_contig_orf(query_id):
    parts = query_id.split("_")
    contig_id = "_".join(parts[:-3])
    return contig_id, query_id

def get_signature(target_id):
    return target_id.split("_")[0]

def process_files(file1_path, file2_path):
    primary_mappings = {}
    results = []
    unique_file2_orfs = set()
    file1_queries = set()
    with open(file1_path, "r") as f:
        next(f)
        for line in f:
            parts = line.strip().split("\\t")
            query_id = parts[0]
            target_id = parts[1]
            file1_queries.add(query_id)
            contig_id, orf_id = parse_contig_orf(query_id)
            signature = get_signature(target_id)
            primary_mappings[query_id] = signature
            results.append((contig_id, orf_id, signature, "Initial"))

    with open(file2_path, "r") as f:
        next(f)
        for line in f:
            parts = line.strip().split("\\t")
            query_id = parts[0]
            target_id = parts[1]
            if query_id not in file1_queries:
                unique_file2_orfs.add(query_id)
            if query_id in primary_mappings:
                continue
            if target_id in primary_mappings:
                contig_id, orf_id = parse_contig_orf(query_id)
                signature = primary_mappings[target_id]
                results.append((contig_id, orf_id, signature, "Recursive"))
    return results, unique_file2_orfs

def main():
    print("DEBUG - Current directory:", os.getcwd())
    print("DEBUG - Files in directory:", os.listdir("."))
    print("DEBUG - Processing files:", file1_path, file2_path)
    results, unique_file2_orfs = process_files(file1_path, file2_path)
    with open(output_file, "w") as f:
        f.write("Genome_ID\\tORF_ID\\tIdentified\\tsignature\\n")
        for contig_id, orf_id, signature, identified in results:
            f.write(f"{contig_id}\\t{orf_id}\\t{identified}\\t{signature}\\n")
    print(f"Results have been saved to: {output_file}")

if __name__ == "__main__":
    main()
    """
}


// Input preparation workflow
workflow {
    // 1. Input channel from config
    Channel
        .fromList(params.proteins)
        .view { config -> "[DEBUG] Input config: ${config}" }
        .set { ch_input }

    // 2. Run PHIDRA
    PHIDRA(ch_input)
        .view { "[DEBUG] PHIDRA output: $it" }
        .branch { protein, fasta, init_search, rec_search ->
            def name = protein.toLowerCase()
            annotation: name in ['polb', 'helicase']
            pasv: name in ['pola', 'rnr']
        }
        .set { branched }

    // 3. Annotation path
    branched.annotation
        .filter { it != null }
        .map { protein, fasta, init_search, rec_search ->
            tuple(protein, init_search, rec_search)
        }
        .view { "[DEBUG] To annotate_top_hits: $it" }
        | annotate_top_hits

    // 4. PASV path
    branched.pasv
        .filter { it != null }
        .map { protein, fasta, init_search, rec_search ->
            tuple(protein, fasta)  // Only emit protein and fasta for PASV
        }
        .view { "[DEBUG] To PASV: $it" }
        | PASV
}