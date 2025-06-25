nextflow.enable.dsl=2

process SPLIT_BY_GENOFEATURE {
    publishDir "${params.outdir}/files_for_embeddings",
        mode: 'copy'

    input:
        tuple path(filtered_tsv), path(input_fasta)

    output:
        path "*", emit: fastas_for_embeddings

    script:
    def selected_list = params.selected_genofeatures.collect { "\'${it}\'" }.join(', ')
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import pandas as pd
    import os

    df = pd.read_csv("${filtered_tsv}", sep='\\t')
    seq_dict = SeqIO.to_dict(SeqIO.parse("${input_fasta}", "fasta"))
    

    selected = [gf.lower() for gf in [${selected_list}]]

    # Group by Protein and Genofeature
    for protein in df['Protein'].unique():
        protein_df = df[df['Protein'] == protein]
        
        for genofeature in protein_df['Genofeature'].unique():
            if genofeature.lower() in selected:
                # Create nested directory structure
                dir_path = os.path.join(protein, genofeature)
                os.makedirs(dir_path, exist_ok=True)
                
                # Filter records for this genofeature
                genofeature_df = protein_df[protein_df['Genofeature'] == genofeature]
                
                # Save filtered TSV
                genofeature_df.to_csv(f"{dir_path}/filtered.tsv", sep='\\t', index=False)
                
                # Extract matching sequences
                matching_seqs = []
                for orf_id in genofeature_df['ORF_ID']:
                    if orf_id in seq_dict:
                        matching_seqs.append(seq_dict[orf_id])
                
                # Write matching sequences to FASTA
                if matching_seqs:
                    SeqIO.write(matching_seqs, f"{dir_path}/{protein}_{genofeature}.fasta", "fasta")
                
                print(f"Processed {protein}/{genofeature}: {len(matching_seqs)} sequences")
    """
}

process EMBEDDINGS {
    publishDir "${params.outdir}/embeddings",
        mode: 'copy'
        
    label 'gpu'  
    conda "/home/nolanv/.conda/envs/esm-umap"
    
    memory { 50.GB }
    
    input:
        path split_dir
        
    output:
        path "*/*/batch_*/**.pt", emit: embeddings
        
    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    extract_script="/mnt/VEIL/tools/embeddings/models/extract.py"
    model_path="/mnt/VEIL/tools/embeddings/models/esm2_t36_3B_UR50D.pt"

    for protein_dir in */; do
        if [ -d "\$protein_dir" ]; then
            echo "Processing protein directory: \$protein_dir"
            protein=\$(basename "\$protein_dir")
            
            for genofeature_dir in "\${protein_dir}"*/; do
                if [ -d "\$genofeature_dir" ]; then
                    echo "Processing genofeature directory: \$genofeature_dir"
                    genofeature=\$(basename "\$genofeature_dir")
                    
                    for fasta in "\${genofeature_dir}"*.fasta; do
                        if [ -f "\$fasta" ]; then
                            # Create final output directory structure
                            output_dir="\${protein}/\${genofeature}"
                            mkdir -p "\$output_dir"
                            
                            echo "Processing FASTA: \$fasta -> \$output_dir"
                            
                            # ESM embeddings 
                            /home/nolanv/.conda/envs/esm-umap/bin/python \\
                                "\$extract_script" \\
                                "\$model_path" \\
                                "\$fasta" \\
                                "\$output_dir" \\
                                --repr_layers 36 \\
                                --include mean || exit 1
                        fi
                    done
                fi
            done
        fi
    done
    """
}

workflow {
    // Create input channels
    ch_filtered_tsv = Channel.fromPath(params.second_run)
    ch_input_fasta = Channel.fromPath(params.input_fasta)
    
    // Combine the input files into a tuple
    ch_inputs = ch_filtered_tsv.combine(ch_input_fasta)

    ch_split_genofeatures = SPLIT_BY_GENOFEATURE(ch_inputs)

    EMBEDDINGS(SPLIT_BY_GENOFEATURE.out.fastas_for_embeddings)
}