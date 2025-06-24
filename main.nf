nextflow.enable.dsl=2


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
    WORK_DIR=\$PWD
    # Clean input fasta headers - only replace dots in headers
    echo "[INFO] Cleaning input FASTA headers..."    
    # sed '/^>/ s/\\./-/g' ${params.input_fasta} > cleaned_input_tmp.fasta

    # Remove asterisks from sequences
    echo "[INFO] Removing asterisks from sequences..."
    # sed -i 's/\\*//g' cleaned_input_tmp.fasta

    # Filter sequences to minimum length 200 using seqkit
    # seqkit seq -m 200 cleaned_input_tmp.fasta -o cleaned_input.fasta

    # Clean subjectDB fasta headers - only replace dots in headers
    # sed '/^>/ s/\\./-/g' ${protein_config.subjectDB} > cleaned_subjectDB_tmp.fasta

    # Remove asterisks from subjectDB sequences
    # sed -i 's/\\*//g' cleaned_subjectDB_tmp.fasta

    # Filter subjectDB sequences to minimum length 200 (optional)
    # seqkit seq -m 200 cleaned_subjectDB_tmp.fasta -o cleaned_subjectDB.fasta


    source /etc/profile.d/conda.sh
    conda activate /home/nolanv/.conda/envs/phidra


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


process PASV {
    conda "/home/nolanv/.conda/envs/phidra"
    cpus 4

    input:
    tuple val(protein), path(fasta), path(align_refs)

    output:
    tuple val(protein), file("pasv_output/output/${protein}_putative.pasv_signatures.tsv")

    script:
    """
    #!/usr/bin/env bash
    set -e  # Exit on error

    # Create directories
    mkdir -p pasv_output/{input,output,pasv}

    # Convert protein name to lowercase for comparison
    protein_lower=\$(echo "${protein}" | tr '[:upper:]' '[:lower:]')

    # Map protein names to PASV expected names and settings
    if [ "\${protein_lower}" == "pola" ]; then
        mapped_name="POL_putative"
        roi_start=521
        roi_end=923
        cat_sites="668,705,758,762"
    elif [ "\${protein_lower}" == "rnr" ]; then
        mapped_name="RNR_putative"
        roi_start=437
        roi_end=625
        cat_sites="437,439,441,462,438"
    else
        echo "Error: Unknown protein ${protein}"
        exit 1
    fi

    # Copy and prepare input files
    cp "${fasta}" "pasv_output/input/\${mapped_name}.fasta"
    cp "${align_refs}" "pasv_output/input/align_refs.fa"

    # Download and prepare PASV if needed
    PASV_DIR="pasv_output/pasv"
    if [ ! -f "\${PASV_DIR}/pasv" ]; then
        mkdir -p "\${PASV_DIR}"
        cd "\${PASV_DIR}"
        wget -q https://github.com/mooreryan/pasv/releases/download/2.0.2/pasv-2.0.2-alpine-static.zip
        unzip pasv-2.0.2-alpine-static.zip
        chmod 755 pasv
        ./pasv --help
        cd -
    fi

    # Run PASV
    \${PASV_DIR}/pasv msa \\
        --outdir=pasv_output/output \\
        --force \\
        --roi-start=\${roi_start} \\
        --roi-end=\${roi_end} \\
        --jobs=${task.cpus} \\
        --aligner=mafft \\
        pasv_output/input/\${mapped_name}.fasta \\
        pasv_output/input/align_refs.fa \\
        \${cat_sites}

    # Rename output file if it exists
    PASV_FILE="pasv_output/output/\${mapped_name}.pasv_signatures.tsv"
    if [ -f "\${PASV_FILE}" ]; then
        mv "\${PASV_FILE}" "pasv_output/output/${protein}_putative.pasv_signatures.tsv"
    else
        echo "Error: PASV did not produce expected output file: \${PASV_FILE}"
        exit 1
    fi

    # Run post-processing (optional)
    # Uncomment and update paths if needed
    # python /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/pasv_post.py \\
    #     pasv_output/output/${protein}_putative.pasv_signatures.tsv \\
    #     pasv_output/output/${protein}_pasv_boxplots.png \\
    #     ${protein}
    """
}


process PASV_post {
    conda "/home/nolanv/.conda/envs/phidra"
    input:
    tuple val(protein), path(pasv_signatures)

    output:
    path "${protein}_pasv_boxplots.png", emit: "pasv_boxplots"

    script:
    """
#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import sys

pasv_signatures = "${pasv_signatures}"
pasv_boxplots = "pasv_output/output/${protein}_pasv_boxplots.png"
protein = "${protein}"

def pasv_post(file):
    df = pd.read_csv(file, sep='\t')

    gapless = df[~df['signature'].str.contains('-')].copy()


    gapless['orf_length'] = gapless['name'].apply(lambda x: (abs(int(x.split('_')[2]) - int(x.split('_')[1])) + 1) // 3)


    def classify_span(s):
        s = s.lower()
        if 'both' in s:
            return 'both'
        elif 'start' in s:
            return 'start'
        elif 'end' in s:
            return 'end'
        elif 'neither' in s:
            return 'neither'
        return 'unknown'
    gapless['span_class'] = gapless['spans'].apply(classify_span)
    tally = gapless.groupby(['signature', 'span_class']).size().reset_index(name='count')
    print("Tally of signature and span class combinations:")
    print(tally)
    return gapless, tally

def boxplots(df, tally, output_file, protein):
    span_classes = df['span_class'].unique()
    sig_counts = [tally[tally['span_class'] == span_class]['signature'].nunique() for span_class in span_classes]
    width_ratios = [max(1, n) *0.8 for n in sig_counts]
    fig_width = sum(width_ratios)

    orf_min = df['orf_length'].min()
    orf_max = df['orf_length'].max()

    fig, axes = plt.subplots(nrows=1, ncols=len(span_classes), figsize=(fig_width, 8), sharey=True, gridspec_kw={'width_ratios': width_ratios})

    if len(span_classes) == 1:
        axes = [axes]
    
    for i, span_class in enumerate(span_classes):
        ax = axes[i] 
        sub = df[df['span_class'] == span_class]
        sub_tally = (
            tally[tally['span_class'] == span_class]
            .sort_values(by='count', ascending=False)
        )
        signature_sort = sub_tally['signature'].tolist()
        data = [sub[sub['signature'] == sig]['orf_length'] for sig in signature_sort]
        ax.boxplot(data, labels=signature_sort)

        ax.yaxis.set_tick_params(labelleft=True)
        ax.set_yticks(ax.get_yticks())

        for j, sig in enumerate(signature_sort):
            values = sub[sub['signature'] == sig]['orf_length']
            count = sub_tally[sub_tally['signature'] == sig]['count'].values
            mean_val = round(values.mean())
            ax.annotate(f"N={count[0] if len(count) > 0 else 0}",
                xy=(j + 1, 0), xycoords=('data', 'axes fraction'),
                xytext=(0, -20), textcoords='offset points',
                ha='center', va='top', fontsize=10, color='blue')
            ax.annotate(f"Mean={mean_val}",
                        xy=(j + 1, values.max()),
                        xytext=(0, 5), textcoords='offset points',
                        ha='center', va='bottom', fontsize=9, color='green')
        ax.set_title(f"Span: {span_class}")
        ax.set_ylabel("ORF Length (aa)")

    fig.suptitle(f"Signature Frequency of {protein}", fontsize=16)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

# Run the functions directly
df, tally = pasv_post("${pasv_signatures}")
boxplots(df, tally, "${protein}_pasv_boxplots.png", "${protein}")
    """
}



process STANDARDIZE_OUTPUTS {
    debug true
    publishDir "${params.outdir}/combined", mode: 'copy'
    conda "/home/nolanv/.conda/envs/phidra"

    input:
    tuple val(proteins), file('*.tsv')  // Accept multiple TSV files as input

    output:
    file "combined_results.tsv"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob
    import os

    # Get list of input TSV files
    input_files = sorted(glob.glob('*.tsv'))
    proteins = "${proteins}".strip('[]').split(', ')
    
    print(f"Processing {len(input_files)} files for proteins: {proteins}")

    def calculate_orf_length(orf_id):
        try:
            parts = orf_id.split('_')
            return (abs(int(parts[2]) - int(parts[1])) + 1) // 3
        except (IndexError, ValueError) as e:
            print(f"Warning: Could not calculate length for {orf_id}: {str(e)}")
            return None

    all_results = []
    
    # Use zip to pair proteins with files directly
    for protein, input_file in zip(proteins, input_files):
        try:
            print(f"Processing {protein} file: {input_file}")
            df = pd.read_csv(input_file, sep='\\t')
            
            if 'Genome_ID' in df.columns and 'Identified' in df.columns:
                print(f"Detected PHIDRA format for {protein}")
                df['ORF_length'] = df['ORF_ID'].apply(calculate_orf_length)
                df['Source'] = 'PHIDRA'
                df['Protein'] = protein
                if 'signature' not in df.columns:
                    df['signature'] = ''
                
            elif 'name' in df.columns and 'spans' in df.columns:
                print(f"Detected PASV format for {protein}")
                df['Genome_ID'] = df['name'].apply(lambda x: '_'.join(x.split('_')[:-3]))
                df['ORF_ID'] = df['name']
                df['ORF_length'] = df['ORF_ID'].apply(calculate_orf_length)
                df['Identified'] = 'PASV_' + df['spans'].fillna('Unknown')
                df['Source'] = 'PASV'
                df['Protein'] = protein
                df['signature'] = df['signature'].fillna('')
            
            df = df[['Genome_ID', 'ORF_ID', 'ORF_length', 'Identified', 'signature', 'Source', 'Protein']]
            all_results.append(df)
            print(f"Successfully processed {protein} with {len(df)} entries")
            
        except Exception as e:
            print(f"Error processing {protein}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue

    if all_results:
        # Combine all results and sort
        final_df = pd.concat(all_results, ignore_index=True)
        final_df = final_df.sort_values(['Protein', 'Genome_ID', 'ORF_ID'])
        
        # Print summary statistics
        print("\\nSummary Statistics:")
        print(f"Total entries: {len(final_df)}")
        print("\\nEntries by Source and Protein:")
        print(final_df.groupby(['Source', 'Protein']).size())
        
        # Write combined results
        final_df.to_csv("combined_results.tsv", sep='\\t', index=False)
        print(f"Written {len(final_df)} total entries to combined_results.tsv")
    else:
        # Create empty file with headers if no results
        pd.DataFrame(columns=['Genome_ID', 'ORF_ID', 'ORF_length', 'Identified', 'signature', 'Source', 'Protein']).to_csv(
            "combined_results.tsv", sep='\\t', index=False)
        print("No results to process, created empty file with headers")
    """
}












process annotate_top_hits {
    debug true

    input:
    tuple val(protein), 
          path(pfam_validated_fasta),
          path(initial_search, stageAs: 'initial_search.tsv'), 
          path(recursive_search, stageAs: 'recursive_search.tsv')

    output:
    tuple val(protein), file("${protein}_annotated_hits.tsv")

    script:
    """
    #!/usr/bin/env python3
import sys
import os
from Bio import SeqIO

protein = "${protein}"
pfam_validated_fasta = "${pfam_validated_fasta}"
file1_path = "initial_search.tsv"
file2_path = "recursive_search.tsv"
output_file = "${protein}_annotated_hits.tsv"

def load_validated_ids():
    validated_ids = set()
    with open(pfam_validated_fasta, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            validated_ids.add(record.id)
    print(f"Loaded {len(validated_ids)} validated sequences from Pfam")
    return validated_ids

def parse_contig_orf(query_id):
    parts = query_id.split("_")
    contig_id = "_".join(parts[:-3])
    return contig_id, query_id

def get_signature(target_id):
    return target_id.split("_")[0]

def process_files(file1_path, file2_path, validated_ids):
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
            if query_id in validated_ids:
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
            if query_id in validated_ids:
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
    validated_ids = load_validated_ids()
    print("DEBUG - Files in directory:", os.listdir("."))
    print("DEBUG - Processing files:", file1_path, file2_path)
    results, unique_file2_orfs = process_files(file1_path, file2_path, validated_ids)
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
    def ch_input = Channel.fromList(params.proteins)


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
    def ch_annotation = branched.annotation
        .filter { it != null }
        .map { protein, fasta, init_search, rec_search ->
            tuple(protein, fasta, init_search, rec_search)
        }
        | annotate_top_hits

    // 4. PASV path
    def ch_pasv = branched.pasv
        .filter { it != null }
        .map { protein, fasta, init_search, rec_search ->
            def proteinConfig = params.proteins.find { it.protein.toLowerCase() == protein.toLowerCase() }
            def align_refs = file(proteinConfig.pasv_align_refs)
            tuple(protein, fasta, align_refs)
        }
        | PASV

    // Standardize results
    ch_annotation
        .mix(ch_pasv)
        .collect()
        .map { results ->
            def file_map = [:] // Create a map of protein -> file
            results.collate(2).each { protein, file ->
                file_map[protein] = file
            }
            tuple(file_map.keySet().toList(), file_map.values().toList())
        }
        .view { "[DEBUG] Files for STANDARDIZE_OUTPUTS: proteins=${it[0]}, files=${it[1]}" }
        | STANDARDIZE_OUTPUTS
}