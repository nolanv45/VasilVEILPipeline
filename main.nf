nextflow.enable.dsl=2

process CLEAN_FASTA_HEADERS {
    tag "${meta.id}"
    conda "/home/nolanv/.conda/envs/phidra"

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("${meta.id}_cleaned.fasta"), emit: fasta

    script:
    """
    # Clean headers
    seqkit replace -p '>' -r '>' ${fasta} | sed '/^>/ s/\\./-/g' > temp.fasta
    seqkit seq -w 0 temp.fasta | sed '/^>/! s/\\*//g' > "${meta.id}_cleaned.fasta"
    rm temp.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/^seqkit v//')
    END_VERSIONS
    """
}

process PHIDRA {
    tag "${meta.id}:${meta.protein}"
    conda "/home/nolanv/.conda/envs/phidra"
    label 'standard'
    publishDir "${params.outdir}/${meta.id}/phidra",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith("pfam_validated_full_protein.fa")) "validated_sequences/${filename}"
            else if (filename.contains("_TopHit_Evalue.tsv")) "mmseqs_results/${filename}"
            else null
        }

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), 
              path("*/final_results/pfam_validated_full_protein.fa"), 
              path("*/mmseqs_results/initial_search/*_TopHit_Evalue.tsv"),
              path("*/mmseqs_results/recursive_search/*_TopHit_Evalue.tsv"),
              emit: results

    script:
    """
    WORK_DIR=\$PWD
    source /etc/profile.d/conda.sh
    conda activate /home/nolanv/.conda/envs/phidra

    cd ${params.phidra_dir}

    INPUT_FASTA=\$WORK_DIR/${fasta}

    python phidra_run.py \\
        -i \$INPUT_FASTA \\
        -db ${meta.subjectDB} \\
        -pfam ${params.pfamDB} \\
        -ida ${meta.pfamDomain} \\
        -f ${meta.protein} \\
        -o \$WORK_DIR \\
        -t ${task.cpus}

    if [ ! -f "\$WORK_DIR/${meta.protein}/mmseqs_results/recursive_search/${meta.protein}_TopHit_Evalue.tsv" ]; then
        mkdir -p "\$WORK_DIR/output/mmseqs_results/recursive_search"
        echo -e "Query_ID\tTarget_ID\tSequence_Identity\tAlignment_Length\tMismatches\tGap_Openings\tQuery_Start\tQuery_End\tTarget_Start\tTarget_End\tEvalue\tBit_Score" > \
            "\$WORK_DIR/${meta.protein}/mmseqs_results/recursive_search/${meta.protein}_TopHit_Evalue.tsv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}

process PASV {
    tag "${meta.id}:${meta.protein}"
    conda "/home/nolanv/.conda/envs/phidra"
    label 'standard'
    publishDir "${params.outdir}/${meta.id}/pasv/${meta.protein}",
        mode: 'copy',
        pattern: "pasv_output/output/*.tsv"

    input:
        tuple val(meta), path(fasta), path(align_refs)

    output:
        tuple val(meta), path("pasv_output/output/${meta.protein}_putative.pasv_signatures.tsv"), emit: results   

    script:
    def protein_lower = meta.protein.toLowerCase()
    """
    mkdir -p pasv_output/{input,output,pasv}

    # Map protein names to settings
    if [ "${protein_lower}" == "pola" ]; then
        mapped_name="${meta.protein}_putative"
        roi_start=521
        roi_end=923
        cat_sites="668,705,758,762"
    elif [ "${protein_lower}" == "rnr" ]; then
        mapped_name="${meta.protein}_putative"
        roi_start=437
        roi_end=625
        cat_sites="437,439,441,462,438"
    else
        echo "Error: Unknown protein ${meta.protein}"
        exit 1
    fi

    # Prepare input files
    cp "${fasta}" "pasv_output/input/\${mapped_name}.fasta"
    cp "${align_refs}" "pasv_output/input/align_refs.fa"

    # Setup PASV
    PASV_DIR="pasv_output/pasv"
    if [ ! -f "\${PASV_DIR}/pasv" ]; then
        mkdir -p "\${PASV_DIR}"
        cd "\${PASV_DIR}"
        wget -q https://github.com/mooreryan/pasv/releases/download/2.0.2/pasv-2.0.2-alpine-static.zip
        unzip pasv-2.0.2-alpine-static.zip
        chmod 755 pasv
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pasv: \$(\${PASV_DIR}/pasv --version | sed 's/pasv //')
    END_VERSIONS
    """
}


process ANALYZE_AND_PLOT {
    tag "${meta.id}"
    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}/${meta.id}", 
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.tsv')) "phidra_analysis/stats/${filename}"
            else if (filename.endsWith('.png')) "phidra_analysis/plots/${filename}"
            else null
        }

    input:
        tuple val(meta), path(input_tsv)

    output:
        tuple val(meta), 
              path("protein_stats.tsv"),
              path("length_distribution.png"),
              emit: results

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
            lambda x: (abs(int(x.split('_')[-3]) - int(x.split('_')[-2])) + 1) // 3
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
    tag "${meta.id}:${meta.protein}"
    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}/${meta.id}/pasv_analysis/${meta.protein}",
        mode: 'copy',
        pattern: "*.{tsv,png}"

    input:
        tuple val(meta), path(pasv_file)

    output:
        tuple val(meta),
              path("${meta.protein}_processed.tsv"),
              path("${meta.protein}_signature_stats.tsv"),
              path("${meta.protein}_signature_distribution.png"),
              emit: results

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Read and process the PASV file
    print(f"Processing PASV file for ${meta.protein}")
    df = pd.read_csv("${pasv_file}", sep='\\t')
    print(f"Initial data shape: {df.shape}")

    # Remove signatures with gaps
    processed_df = df[~df['signature'].str.contains('-')].copy()
    print(f"After removing gaps: {processed_df.shape}")

    # Calculate ORF length
    processed_df['orf_length'] = processed_df['name'].apply(
        lambda x: (abs(int(x.split('_')[2]) - int(x.split('_')[1])) + 1) // 3
    )

    # Create span classes based on spans_start and spans_end columns
    processed_df['span_class'] = processed_df.apply(
        lambda row: 'both' if row['spans_start'] == 'Yes' and row['spans_end'] == 'Yes'
                   else 'start' if row['spans_start'] == 'Yes'
                   else 'end' if row['spans_end'] == 'Yes'
                   else 'neither',
        axis=1
    )

    # Generate statistics
    stats = processed_df.groupby(['signature', 'span_class']).agg({
        'orf_length': ['count', 'mean', 'std', 'min', 'max']
    }).reset_index()
    stats.columns = ['signature', 'span_class', 'count', 'mean_length', 'std_length', 'min_length', 'max_length']
    
    # Sort signatures by count for better visualization
    signature_order = stats.groupby('signature')['count'].sum().sort_values(ascending=True).index

    # Create visualization
    span_classes = sorted(processed_df['span_class'].unique())
    n_signatures = len(signature_order)
    
    # Set figure dimensions
    height_per_sig = 0.4
    fig_height = max(8, n_signatures * height_per_sig)
    fig_width = len(span_classes) * 12
    
    # Create figure
    fig, axes = plt.subplots(1, len(span_classes), 
                            figsize=(fig_width, fig_height),
                            sharey=True)
    
    if len(span_classes) == 1:
        axes = [axes]
    
    print(f"Creating plots for {len(span_classes)} span classes and {n_signatures} signatures")
    
    # Plot for each span class
    for i, span in enumerate(span_classes):
        ax = axes[i]
        span_data = processed_df[processed_df['span_class'] == span]
        
        if not span_data.empty:
            # Create boxplot
            sns.boxplot(data=span_data,
                       x='orf_length',
                       y='signature',
                       order=signature_order,
                       orient='h',
                       ax=ax,
                       color='skyblue')

            # Add statistics
            span_stats = stats[stats['span_class'] == span]
            for idx, sig in enumerate(signature_order):
                sig_stats = span_stats[span_stats['signature'] == sig]
                if not sig_stats.empty:
                    # Add count
                    ax.text(ax.get_xlim()[0], idx, 
                           f"n={int(sig_stats['count'].iloc[0])}",
                           ha='right', va='center',
                           fontsize=8, color='blue')
                    
                    # Add mean
                    if sig_stats['count'].iloc[0] > 0:
                        ax.text(ax.get_xlim()[1], idx,
                               f"μ={sig_stats['mean_length'].iloc[0]:.1f}",
                               ha='left', va='center',
                               fontsize=8, color='green')
        
        ax.set_title(f"Span: {span}")
        ax.set_xlabel("ORF Length (aa)")
        ax.grid(True, alpha=0.3)
        
        # Ensure y-labels are readable
        ax.tick_params(axis='y', labelsize=8, pad=40)
    
    # Adjust layout
    plt.subplots_adjust(left=0.25, right=0.95, bottom=0.1, top=0.9, wspace=0.2)
    fig.suptitle(f"${meta.protein} Signature Distribution", fontsize=16, y=1.02)
    
    # Save outputs
    processed_df.to_csv("${meta.protein}_processed.tsv", sep='\\t', index=False)
    stats.to_csv("${meta.protein}_signature_stats.tsv", sep='\\t', index=False)
    plt.savefig("${meta.protein}_signature_distribution.png", bbox_inches='tight', dpi=300)
    plt.close()
    """
}

process STANDARDIZE_OUTPUTS {
    tag "${meta.id}"

    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}/${meta.id}", 
        mode: 'copy',
        pattern: "*.tsv"

    input:
    tuple val(meta), path(tsv_files)

    output:
    tuple val(meta), path("${meta.id}_combined_results.tsv"), emit: standardized

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os

    required_columns = ['contig_id', 'orf_id', 'identified', 'genofeature', 'protein', 'dataset']
    dfs = []

    print(f"Processing files for ${meta.id}: {repr('${tsv_files}'.split())}")

    for tsv_file in "${tsv_files}".split():
        if not os.path.exists(tsv_file):
            print(f"File not found: {tsv_file}")
            continue

        print(f"Reading file: {tsv_file}")
        try:
            df = pd.read_csv(tsv_file, sep='\\t')
            print(f"Columns in {tsv_file}: {df.columns.tolist()}")

            # Detect file type and standardize
            if 'name' in df.columns and 'signature' in df.columns:  # PASV format
                print(f"Converting PASV format: {tsv_file}")
                protein = os.path.basename(tsv_file).split('_')[0]
                
                standardized = pd.DataFrame({
                    'contig_id': df['name'].apply(lambda x: '_'.join(x.split('_')[:-3])),
                    'orf_id': df['name'],
                    'identified': 'PASV',
                    'genofeature': df['signature'],
                    'protein': protein,
                    'dataset': "${meta.id}"
                })
                print(f"Converted {len(standardized)} PASV entries")
                dfs.append(standardized)
            
            elif 'genome_id' in df.columns:  # Standard format with old column name
                print(f"Processing standard format: {tsv_file}")
                df = df.rename(columns={'genome_id': 'contig_id'})  # Rename the column
                if 'dataset' not in df.columns:
                    df['dataset'] = "${meta.id}"
                dfs.append(df)
            
            else:
                print(f"Warning: Unrecognized format in {tsv_file}")
                print(f"Found columns: {df.columns.tolist()}")
                continue

        except Exception as e:
            print(f"Error processing {tsv_file}: {str(e)}")
            continue

    if dfs:
        print("Combining DataFrames...")
        combined = pd.concat(dfs, ignore_index=True)
        
        # Ensure all required columns exist
        for col in required_columns:
            if col not in combined:
                combined[col] = None
        
        # Reorder columns
        combined = combined[required_columns]
        
        print(f"Final combined shape: {combined.shape}")
        print(f"Final columns: {combined.columns.tolist()}")
        combined.to_csv("${meta.id}_combined_results.tsv", sep='\\t', index=False)
    else:
        print("No data to combine, creating empty output")
        pd.DataFrame(columns=required_columns).to_csv("${meta.id}_combined_results.tsv", sep='\\t', index=False)
    """
}


process PHIDRA_ONLY_SUMMARY {
    tag "${meta.id}:${meta.protein}"
    publishDir "${params.outdir}/${meta.id}/phidra/phidra_only/${protein}", 
        mode: 'copy',
        pattern: "*.tsv"
    conda "/home/nolanv/.conda/envs/phidra"

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("${meta.protein}_phidra_only.tsv"), emit: results

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    
    protein = "${meta.protein}"

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
    df.to_csv("${meta.protein}_phidra_only.tsv", sep='\\t', index=False)

    print(f"Processed {len(results)} ORFs with no hits for ${meta.protein}")
    """
}



process ANNOTATE_HITS {
    tag "${meta.id}:${meta.protein}"
    publishDir "${params.outdir}/${meta.id}/phidra/annotate_hits/${meta.protein}", 
        mode: 'copy',
        pattern: "*_annotated_hits.tsv"

    input:
    tuple val(meta), 
          path(pfam_validated_fasta),
          path(initial_search, stageAs: 'initial_search.tsv'), 
          path(recursive_search, stageAs: 'recursive_search.tsv')

    output:
    tuple val(meta), file("${meta.protein}_annotated_hits.tsv"), emit: results

    script:
    """
#!/usr/bin/env python3
import sys
import os
from Bio import SeqIO

protein = "${meta.protein}"
pfam_validated_fasta = "${pfam_validated_fasta}"
file1_path = "initial_search.tsv"
file2_path = "recursive_search.tsv"
output_file = "${meta.protein}_annotated_hits.tsv"

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
    tag "${meta.id}"
    debug true
    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}/${meta.id}", 
        mode: 'copy',
        pattern: "*.tsv"

    input:
        tuple val(meta), path(input_file)    // Properly declare the input file

    output:
        tuple val(meta), path("duplicate_orfs.tsv")   // Declare the output file

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    
    # Read the input file with proper quoting and correct separator
    df = pd.read_csv("${input_file}", sep='\\t')

    # Find duplicates based on both contig_id and orf_id
    duplicates = df[df.duplicated(subset=['contig_id', 'orf_id'], keep=False)]
    
    # Sort duplicates for better readability
    duplicates = duplicates.sort_values(['contig_id', 'orf_id'])

    # Save duplicates with tab separator
    duplicates.to_csv('duplicate_orfs.tsv', sep='\\t', index=False)
    
    # Print summary for debugging
    print(f"Found {len(duplicates)} duplicate entries")
    """
}



process COMBINE_PHIDRA_TSV {
    tag "${meta.id}"
    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}/${meta.id}/phidra_analysis", 
        mode: 'copy'

    input:
    tuple val(meta), path(tsv_files)

    output:
        tuple val(meta), path('combined_phidra_output.tsv'), emit: combined

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os


    dfs = []
    for tsv_file in "${tsv_files}".split():
        if os.path.exists(tsv_file):
            print(f"Reading file: {tsv_file}")
            try:
                df = pd.read_csv(tsv_file, sep='\\t')
                if not df.empty:
                    dfs.append(df)
                    print(f"Added {len(df)} rows from {tsv_file}")
                else:
                    print(f"Warning: Empty DataFrame from {tsv_file}")
            except Exception as e:
                print(f"Error reading {tsv_file}: {str(e)}")
        else:
            print(f"Warning: File not found: {tsv_file}")


    if dfs:
        combined = pd.concat(dfs, ignore_index=True)
        print(f"Combined DataFrame has {len(combined)} rows")
        combined.to_csv("combined_phidra_output.tsv", sep='\\t', index=False)
    else:
        print("No data to combine, creating empty output")
        pd.DataFrame(columns=['genome_id', 'orf_id', 'identified', 'genofeature', 'protein']).to_csv("combined_phidra_output.tsv", sep='\\t', index=False)
    """
}


process SPLIT_BY_GENOFEATURE {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}",
        mode: 'copy'

    input:
        tuple val(meta), path(filtered_tsv)

    output:
        tuple val(meta), path("files_for_embeddings"), emit: fastas_for_embeddings

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import pandas as pd
    import os

    # Create the base output directory
    os.makedirs("files_for_embeddings", exist_ok=True)
    df = pd.read_csv("${filtered_tsv}", sep='\\t')
    seq_dict = SeqIO.to_dict(SeqIO.parse("${meta.cleaned_fasta}", "fasta"))

    # Group by protein and genofeature
    for protein in df['protein'].unique():
        protein_df = df[df['protein'] == protein]
        
        for genofeature in protein_df['genofeature'].unique():
            # Create nested directory structure
            dir_path = os.path.join("files_for_embeddings", protein, genofeature)
            os.makedirs(dir_path, exist_ok=True)
            
            # Filter records for this genofeature
            genofeature_df = protein_df[protein_df['genofeature'] == genofeature]
            
            # Save filtered TSV
            genofeature_df.to_csv(f"{dir_path}/filtered.tsv", sep='\\t', index=False)
            
            # Extract matching sequences
            matching_seqs = []
            for orf_id in genofeature_df['orf_id']:
                if orf_id in seq_dict:
                    matching_seqs.append(seq_dict[orf_id])
            
            # Write matching sequences to FASTA
            if matching_seqs:
                SeqIO.write(matching_seqs, f"{dir_path}/{protein}_{genofeature}.fasta", "fasta")
    """
}

process COMBINE_DATASETS {

    conda "/home/nolanv/.conda/envs/phidra"
    publishDir "${params.outdir}",
        mode: 'copy',
        pattern: "*.tsv"

    input:
        path(tsv_files)

    output:
        path "combined_datasets.tsv", emit: combined

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os

    dfs = []
    input_files = "${tsv_files}".split()

    for tsv_file in input_files:
        if os.path.exists(tsv_file):
            print(f"Reading file: {tsv_file}")
            try:
                df = pd.read_csv(tsv_file, sep='\\t')
                if not df.empty:
                    dfs.append(df)
                    print(f"Added {len(df)} rows from {tsv_file}")
                else:
                    print(f"Warning: Empty DataFrame from {tsv_file}")
            except Exception as e:
                print(f"Error reading {tsv_file}: {str(e)}")
        else:
            print(f"Warning: File not found: {tsv_file}")

    if dfs:
        combined = pd.concat(dfs, ignore_index=True)
        print(f"Combined DataFrame has {len(combined)} rows")
        combined.to_csv("combined_datasets.tsv", sep='\\t', index=False)
    else:
        print("No data to combine, creating empty output")
        pd.DataFrame(columns=['contig_id', 'orf_id', 'identified', 'genofeature', 'protein', 'dataset']).to_csv("combined_datasets.tsv", sep='\\t', index=False)
    """
}

workflow {
    // Create dataset channels
    Channel
        .fromList(params.datasets.entrySet())
        .map { entry -> 
            def meta = [
                id: entry.key,      // Dataset identifier
                path: entry.value   // Dataset file path
            ]
            tuple(meta, meta.path)
        }
        .set { ch_datasets }

    // Process each dataset through the pipeline
    ch_datasets | PROCESS_DATASET
}


workflow PROCESS_DATASET {
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
            .branch { meta, fasta, init_search, rec_search -> 
                annotation: meta.protein.toLowerCase() in ['polb', 'helicase']
                pasv: meta.protein.toLowerCase() in ['pola', 'rnr']
                phidra_only: true
            }


        ch_pasv = ch_branched.pasv
            .map { meta, fasta, init_search, rec_search ->
                tuple(
                    meta,
                    fasta,
                    meta.pasv_align_refs // align_refs is already in meta from earlier
                )
            }

        ANNOTATE_HITS(ch_branched.annotation)
        PASV(ch_pasv)

        PASV_POST(PASV.out.results)
        PHIDRA_ONLY_SUMMARY(ch_branched.phidra_only)

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

        COMBINE_PHIDRA_TSV(ch_phidra_results).view()

        ANALYZE_AND_PLOT(COMBINE_PHIDRA_TSV.out.combined)


        ch_pasv_map = PASV_POST.out
            .map { meta, processed_file, stats_file, plot_file ->
                tuple(meta, processed_file)
            }



        ch_combined_phidra = COMBINE_PHIDRA_TSV.out.combined ?: Channel.empty()
        ch_pasv_map = ch_pasv_map ?: Channel.empty()

        ch_files_to_standardize = ch_combined_phidra
            .mix(ch_pasv_map)
            .map { meta, file -> 
                tuple(meta.id, meta, file)
            }
            .groupTuple(by: 0)
            .map { id, metas, files ->
                tuple(metas[0], files.flatten())  // Flatten the files array
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
}
