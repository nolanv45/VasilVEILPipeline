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

process cleanHeaders {

}

process annotate_top_hits {
    input:
    tuple path(initial_search, stageAs: { "initial_search_${protein_name}.tsv" }), \
          path(recursive_search, stageAs: { "recursive_search_${protein_name}.tsv" }), \
          val(output_file), \
          val(protein_name)

    output:
    path "${output_file}"

    script:
    """
    #!/usr/bin/env python
    import sys

    def parse_contig_orf(query_id):
        parts = query_id.split('_')
        contig_id = '_'.join(parts[:-3])
        return contig_id, query_id

    def get_signature(target_id):
        return target_id.split('_')[0]

    def process_files(file1_path, file2_path):
        primary_mappings = {}
        results = []
        unique_file2_orfs = set()
    
        # Process first file and store all its Query_IDs
        file1_queries = set()
        with open(file1_path, 'r') as f:
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                query_id = parts[0]
                target_id = parts[1]
                file1_queries.add(query_id)
                contig_id, orf_id = parse_contig_orf(query_id)
                signature = get_signature(target_id)
                primary_mappings[query_id] = signature
                # Add "Initial" for entries from first file
                results.append((contig_id, orf_id, signature, "Initial"))
    
        # Process second file
        with open(file2_path, 'r') as f:
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                query_id = parts[0]
                target_id = parts[1]
                # If this Query_ID isn't in file1, add it to unique ORFs
                if query_id not in file1_queries:
                    unique_file2_orfs.add(query_id)
                # Skip if query_id already exists in primary mappings
                if query_id in primary_mappings:
                    continue
                # Look up the target_id in primary mappings
                if target_id in primary_mappings:
                    contig_id, orf_id = parse_contig_orf(query_id)
                    signature = primary_mappings[target_id]
                    # Add "Recursive" for entries from second file
                    results.append((contig_id, orf_id, signature, "Recursive"))
    
        return results, unique_file2_orfs

    def main():   
        file1_path = "${initial_search}"
        file2_path = "${recursive_search}"
        output_file = "${output_file}"    
        # Process files
        results, unique_file2_orfs = process_files(file1_path, file2_path)
        # Write results to file
        with open(output_file, 'w') as f:
            # Updated header with new column
            f.write("Genome_ID\tORF_ID\tIdentified\tsignature\n")
            # Write data with identification source
            for contig_id, orf_id, signature, identified in results:
                f.write(f"{contig_id}\t{orf_id}\t{identified}\t{signature}\n")
        print(f"\nResults have been saved to: {output_file}")
    
    if __name__ == "__main__":
        main()
    """
}

workflow {
    // Channel.fromList(params.proteins)
    //     .set { protein_configs }
    // phidra_results = protein_configs | PHIDRA
    // phidra_results.collect().set { all_phidra_outputs }

    // // Run PASV only if PolA or RNR is involved
    // def proteins_lower = params.proteins.collect { it.toString().toLowerCase() }
    // if (proteins_lower.any { it == "pola" || it == "rnr" }) {
    //     PASV(all_phidra_outputs)
    // }







        annotation_jobs = params.proteins
        .findAll { p -> 
            def lower = p.protein.toString().toLowerCase()
            lower == "helicase" || lower == "polb"
        }
        .collect { p ->
            def init_path = "phidra_output_test/${protein}/mmseqs_results/initial_search/${protein}_TopHit_Evalue.tsv"
            def rec_path = "phidra_output_test/${protein}/mmseqs_results/initial_search/${protein}_Recursive_TopHit_Evalue.tsv"
            def out_path = "phidra_output_test/${protein}_annotated_hits.tsv"
            tuple(init_path, rec_path, out_path, p.protein)
        }

    Channel
        .from(annotation_jobs)
        .set { annotation_channel }

    annotation_channel
        | annotate_top_hits
}
