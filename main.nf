nextflow.enable.dsl=2


process PHIDRA {
    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}/phidra", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith("pfam_validated_full_protein.fa")) "validated_sequences/${filename}"
        else if (filename.contains("_TopHit_Evalue.tsv")) "mmseqs_results/${filename}"
        else null
    }

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
    publishDir "${params.outdir}/pasv/${protein}", 
        mode: 'copy',
        pattern: "pasv_output/output/*.tsv"

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
        mapped_name="${protein}_putative"
        roi_start=521
        roi_end=923
        cat_sites="668,705,758,762"
    elif [ "\${protein_lower}" == "rnr" ]; then
        mapped_name="${protein}_putative"
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
    # PASV_FILE="pasv_output/output/\${mapped_name}.pasv_signatures.tsv"
    #if [ -f "\${PASV_FILE}" ]; then
     #   mv "\${PASV_FILE}" "pasv_output/output/${protein}_putative.pasv_signatures.tsv"
    #else
     #   echo "Error: PASV did not produce expected output file: \${PASV_FILE}"
      #  exit 1
    #fi
    """
}


process ANALYZE_AND_PLOT {
    debug true
    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}", 
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.tsv')) "phidra_analysis/stats/${filename}"
            else if (filename.endsWith('.png')) "phidra_analysis/plots/${filename}"
            else null
        }

    input:
    path(input_tsv)

    output:
    path "protein_stats.tsv"
    path "length_distribution.png"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Read and process input data
    df = pd.read_csv("${input_tsv}", sep='\\t')
    print(f"Processing data with {len(df)} entries")

    # Calculate ORF length if not present
    if 'ORF_length' not in df.columns:
        df['ORF_length'] = df['orf_id'].apply(
            lambda x: (abs(int(x.split('_')[2]) - int(x.split('_')[1])) + 1) // 3
        )

    # Generate basic statistics
    stats = df.groupby('genofeature').agg({
        'ORF_length': ['count', 'mean', 'min', 'max']
    }).reset_index()
    
    # Flatten column names
    stats.columns = ['genofeature', 'count', 'mean_length', 'min_length', 'max_length']
    stats = stats.sort_values('genofeature')
    
    # Save statistics
    stats.to_csv("protein_stats.tsv", sep='\\t', index=False)

    # Create plot
    plt.figure(figsize=(10, max(6, len(stats) * 0.5)))
    
    # Create boxplot using seaborn
    sns.boxplot(data=df, x='ORF_length', y='genofeature', 
               orient='h', color='skyblue')
    
    # Add annotations
    for i, row in stats.iterrows():
        # Count annotation (left)
        plt.annotate(f"N={int(row['count'])}",
            xy=(row['min_length'], i),
            xytext=(-10, 0),
            textcoords='offset points',
            ha='right', va='center',
            fontsize=9, color='blue')
        
        # Mean annotation (right)
        plt.annotate(f"Mean={row['mean_length']:.1f}",
            xy=(row['max_length'], i),
            xytext=(10, 0),
            textcoords='offset points',
            ha='left', va='center',
            fontsize=9, color='green')
    
    plt.title("Length Distribution by genofeature")
    plt.xlabel("ORF Length (aa)")
    plt.grid(True, axis='x', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig("length_distribution.png", bbox_inches='tight', dpi=300)
    plt.close()
    """
}


process PASV_POST {
    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}/pasv_analysis/${protein}",
        mode: 'copy',
        pattern: "*.{tsv,png}"

    input:
    tuple val(protein), file(pasv_file)

    output:
    tuple val(protein), 
        file("${protein}_processed.tsv"),
        file("${protein}_signature_stats.tsv"),
        file("${protein}_signature_distribution.png")

    script:
    """
#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import sys

def pasv_post(file, protein):
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

    # Generate signature statistics
    sig_stats = gapless.groupby('signature').agg({
        'orf_length': ['count', 'mean']
    }).reset_index()
    sig_stats.columns = ['signature', 'count', 'mean_length']
    sig_stats['protein'] = protein
        
    return gapless, tally, sig_stats


def boxplots(df, tally, output_file, protein):
    span_classes = sorted(df['span_class'].unique())
    sig_counts = [
        tally[tally['span_class'] == span_class]['signature'].nunique()
        for span_class in span_classes
    ]
    
    # Calculate statistics for TSV while making plot
    stats_data = []
    
    for span_class in span_classes:
        sub = df[df['span_class'] == span_class]
        sub_tally = (
            tally[tally['span_class'] == span_class]
            .sort_values(by='count', ascending=True)
        )
        
        for sig in sub_tally['signature'].tolist():
            values = sub[sub['signature'] == sig]['orf_length']
            count = sub_tally[sub_tally['signature'] == sig]['count'].values[0]
            mean_val = round(values.mean())
            min_val = values.min()
            max_val = values.max()
            
            stats_data.append({
                'protein': protein,
                'signature': sig,
                'span_class': span_class,
                'count': count,
                'mean_length': mean_val,
                'min_length': min_val,
                'max_length': max_val
            })
    
    # Create and save statistics DataFrame
    stats_df = pd.DataFrame(stats_data)
    stats_df.to_csv(output_file.replace('_distribution.png', '_stats.tsv'), 
                    sep='\t', index=False)
    
    # Create plot with existing code
    height_per_sig = 0.4
    fig_height = max(8, max(sig_counts) * height_per_sig)
    fig_width = len(span_classes) * 6
    
    fig, axes = plt.subplots(nrows=1, ncols=len(span_classes), 
                            figsize=(fig_width, fig_height), sharey=True)

    if len(span_classes) == 1:
        axes = [axes]
    
    for i, span_class in enumerate(span_classes):
        ax = axes[i] 
        sub = df[df['span_class'] == span_class]

        sub_tally = (
            tally[tally['span_class'] == span_class]
            .sort_values(by='count', ascending=True)  # Changed to True for bottom-to-top ordering
        )
        signature_sort = sub_tally['signature'].tolist()
        data = [sub[sub['signature'] == sig]['orf_length'] for sig in signature_sort]
        
        # Create horizontal boxplot
        bp = ax.boxplot(data, positions=range(len(signature_sort)), 
                       vert=False,  # Make boxplot horizontal
                       patch_artist=True)  # Fill boxes with color
        
        # Set y-axis labels (signatures)
        ax.set_yticks(range(len(signature_sort)))
        ax.set_yticklabels(signature_sort)
        
        # Add counts and means
        for j, sig in enumerate(signature_sort):
            values = sub[sub['signature'] == sig]['orf_length']
            count = sub_tally[sub_tally['signature'] == sig]['count'].values
            mean_val = round(values.mean())
            
            # Add count annotation
            ax.annotate(f"N={count[0] if len(count) > 0 else 0}",
                xy=(values.min(), j),  # Place at start of boxplot
                xytext=(-10, 0),  # Offset to the left
                textcoords='offset points',
                ha='right', va='center',
                fontsize=9, color='blue')
            
            # Add mean annotation
            ax.annotate(f"Mean={mean_val}",
                xy=(values.max(), j),  # Place at end of boxplot
                xytext=(10, 0),  # Offset to the right
                textcoords='offset points',
                ha='left', va='center',
                fontsize=9, color='green')
        
        # Customize appearance
        ax.set_title(f"Span: {span_class}")
        ax.set_xlabel("ORF Length (aa)")
        ax.grid(True, axis='x', linestyle='--', alpha=0.7)

    # Adjust layout
    fig.suptitle(f"Signature Distribution of {protein}", fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()


# Process the PASV file
protein = "${protein}"
processed_df, tally, sig_stats = pasv_post("${pasv_file}", protein)

# Save processed results
processed_df.to_csv(f"{protein}_processed.tsv", sep='\\t', index=False)
sig_stats.to_csv(f"{protein}_signature_stats.tsv", sep='\\t', index=False)

# Generate plots
boxplots(processed_df, tally, f"{protein}_signature_distribution.png", protein)

print(f"\\nSignature Statistics Summary for {protein}:")
print(sig_stats.sort_values('count', ascending=False))
    """
}

process STANDARDIZE_OUTPUTS {
    debug true
    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}", 
        mode: 'copy',
        pattern: "*.tsv"

    input:
    tuple val(tag), path(tsv_files)

    output:
    path "combined_results.tsv"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob
    import os

    print("Processing input files:", "${tsv_files}".split())
    
    # Initialize empty list to store DataFrames
    dfs = []
    
    # Process each input file
    for tsv_file in "${tsv_files}".split():
        if not os.path.exists(tsv_file):
            print(f"File not found: {tsv_file}")
            continue
            
        print(f"Reading file: {tsv_file}")
        df = pd.read_csv(tsv_file, sep='\\t')
        
        # Detect if this is a PASV format file by checking columns
        pasv_columns = ['name', 'signature', 'spans']
        is_pasv = all(col in df.columns for col in pasv_columns)
        
        if is_pasv:
            print(f"Converting PASV format file: {tsv_file}")
            protein = os.path.basename(tsv_file).split('_')[0]
            
            # Convert PASV format to standard format
            standardized = pd.DataFrame({
                'genome_id': df['name'].apply(lambda x: '_'.join(x.split('_')[:-3])),
                'orf_id': df['name'],
                'identified': 'PASV',
                'genofeature': df['signature'],
                'protein': protein
            })
            print(f"Converted {len(standardized)} PASV entries")
            dfs.append(standardized)
        else:
            # Verify standard format
            standard_columns = ['genome_id', 'orf_id', 'identified', 'genofeature', 'protein']
            if all(col in df.columns for col in standard_columns):
                print(f"Adding standard format file: {tsv_file}")
                dfs.append(df)
            else:
                print(f"Warning: File {tsv_file} has unexpected format. Columns: {df.columns.tolist()}")
        
        print(f"Processed {len(df)} rows from {tsv_file}")

    # Combine all DataFrames if any were read
    if dfs:
        print(f"Combining {len(dfs)} DataFrames")
        combined = pd.concat(dfs, ignore_index=True)
        print(f"Total rows: {len(combined)}")
        print("Final column names:", combined.columns.tolist())
        combined.to_csv("combined_results.tsv", sep='\\t', index=False)
        print(f"Created combined file with {len(combined)} total rows")
    else:
        print("No data to combine")
        pd.DataFrame(columns=['genome_id', 'orf_id', 'identified', 'genofeature', 'protein']).to_csv("combined_results.tsv", sep='\\t', index=False)
    """
}




process phidra_only_summary {
    debug true
    publishDir "${params.outdir}/phidra/phidra_only/${protein}", 
        mode: 'copy',
        pattern: "*.tsv"
    conda "/home/nolanv/.conda/envs/phidra"

    input:
    tuple val(protein), file(fasta)

    output:
    tuple val(protein), file("${protein}_phidra_only.tsv")

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    
    protein = "${protein}"

    # Process FASTA file to get ORF information
    orfs = []
    current_header = None
    with open("${fasta}") as f:
        for line in f:
            if line.startswith('>'):
                current_header = line.strip()[1:]
                orfs.append(current_header)
    
    # Create DataFrame with same structure as annotate_top_hits output
    results = []
    for orf in orfs:
        genome_id = '_'.join(orf.split('_')[:-3])
        results.append({
            'genome_id': genome_id,
            'orf_id': orf,
            'identified': 'phidra_only',
            'genofeature': protein,  # Empty as specified
            'protein': protein
        })
    
    # Create and save DataFrame
    df = pd.DataFrame(results)
    df.to_csv("${protein}_phidra_only.tsv", sep='\\t', index=False)

    print(f"Processed {len(results)} ORFs with no hits for ${protein}")
    """
}



process annotate_top_hits {
    debug true
    publishDir "${params.outdir}/phidra/annotate_hits/${protein}", 
        mode: 'copy',
        pattern: "*_annotated_hits.tsv"

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
        f.write("genome_id\\torf_id\\tidentified\\tgenofeature\\tprotein\\n")
        for contig_id, orf_id, signature, identified in results:
            f.write(f"{contig_id}\\t{orf_id}\\t{identified}\\t{signature}\\t{protein}\\n")
    print(f"Results have been saved to: {output_file}")

if __name__ == "__main__":
    main()
    """
}


process DUPLICATE_HANDLE {
    debug true
    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}", 
        mode: 'copy',
        pattern: "*.tsv"

    input:
        path input_file    // Properly declare the input file

    output:
        path "duplicate_orfs.tsv"   // Declare the output file

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    
    # Read the input file with proper quoting and correct separator
    df = pd.read_csv("${input_file}", sep='\\t')
    
    # Find duplicates based on both genome_id and orf_id
    duplicates = df[df.duplicated(subset=['genome_id', 'orf_id'], keep=False)]
    
    # Sort duplicates for better readability
    duplicates = duplicates.sort_values(['genome_id', 'orf_id'])
    
    # Save duplicates with tab separator
    duplicates.to_csv('duplicate_orfs.tsv', sep='\\t', index=False)
    
    # Print summary for debugging
    print(f"Found {len(duplicates)} duplicate entries")
    """
}



process COMBINE_PHIDRA_TSV {
    debug true
    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}/phidra_analysis", 
        mode: 'copy',
        pattern: "*.tsv"

    input:
    tuple val(tag), path(tsv_files)

    output:
    path 'combined_phidra_output.tsv'

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob
    import os

    print("Debug: Current directory contents:")
    print(os.listdir('.'))

    # Get list of input files
    tsv_files = "${tsv_files}".split()  
    print(f"Processing TSV files: {tsv_files}")

    dfs = []
    for tsv_file in tsv_files:
        if os.path.exists(tsv_file):
            print(f"Reading file: {tsv_file}")
            df = pd.read_csv(tsv_file, sep='\\t')
            dfs.append(df)
            print(f"Added {len(df)} rows from {tsv_file}")
        else:
            print(f"Warning: File not found: {tsv_file}")

    if dfs:
        # Combine all DataFrames into one
        combined = pd.concat(dfs, ignore_index=True)
        combined.to_csv("combined_phidra_output.tsv", sep='\\t', index=False)
        print(f"Debug: Written combined file with {len(combined)} total rows")
    else:
        print("Warning: No data to combine")
        pd.DataFrame().to_csv("combined_phidra_output.tsv", sep='\\t', index=False)
    """
}

workflow {
    def ch_input = Channel.fromList(params.proteins)

    PHIDRA(ch_input)
        .view { "[DEBUG] PHIDRA output: $it" }
        .branch { protein, fasta, init_search, rec_search ->
            def name = protein.toLowerCase()
            annotation: name in ['polb', 'helicase']
            pasv: name in ['pola', 'rnr']
            phidra_only: true
        }
        .set { branched }

    def ch_annotation = branched.annotation
        .filter { it != null }
        .map { protein, fasta, init_search, rec_search ->
            tuple(protein, fasta, init_search, rec_search)
        }
        | annotate_top_hits
        | map { protein, file -> file }

    def ch_pasv_results = branched.pasv
        .filter { it != null }
        .map { protein, fasta, init_search, rec_search ->
            def proteinConfig = params.proteins.find { it.protein.toLowerCase() == protein.toLowerCase() }
            def align_refs = file(proteinConfig.pasv_align_refs)
            tuple(protein, fasta, align_refs)
        }
        | PASV
        | PASV_POST

    def ch_phidra_only = branched.phidra_only
        .filter { it != null }
        .map { protein, fasta, init_search, rec_search ->
            tuple(protein, fasta)
        }
        | phidra_only_summary
        | map { protein, file -> file }

    def ch_pasv_processed = ch_pasv_results
        .map { protein, processed_file, stats, plot ->
            // Return only the processed file for combining with annotation results
            processed_file
        }

    def ch_combined_phidra_results = Channel.empty()
        .mix(ch_annotation, ch_phidra_only)
        .collect()
        .map { files -> tuple('combined', files) }  
        .view { "[DEBUG] Processing files: ${it[1]}" }
        | COMBINE_PHIDRA_TSV

    // Use the same channel output for both processes
    ch_combined_phidra_results | ANALYZE_AND_PLOT

    // Mix with PASV results and standardize
    def ch_standardized_output = Channel.empty()
        .mix(
            ch_combined_phidra_results,
            ch_pasv_processed
        )
        .collect()
        .map { files -> 
            def flattened = files.flatten()
            println "[DEBUG] Files to standardize: ${flattened}"
            tuple('combined', flattened) 
        }
        | STANDARDIZE_OUTPUTS

    def ch_duplicated_output = ch_standardized_output
        | DUPLICATE_HANDLE
}