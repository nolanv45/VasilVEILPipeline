process SPLIT_BY_GENOFEATURE {
    publishDir "${params.outdir}/split_by_genofeature",
        mode: 'copy'

    input:
        tuple path(filtered_tsv), path(input_fasta)

    output:
        path "*/sequences.fasta", emit: fastas
        path "*/filtered.tsv", emit: tsvs

    script:
    def selected_list = params.selected_genofeatures.collect { "\'${it}\'" }.join(', ')
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import pandas as pd
    import os

    # Read the filtered TSV and input FASTA
    df = pd.read_csv("${filtered_tsv}", sep='\\t')
    seq_dict = SeqIO.to_dict(SeqIO.parse("${input_fasta}", "fasta"))
    
    # Convert selected genofeatures to lowercase for case-insensitive comparison
    selected = [gf.lower() for gf in [${selected_list}]]
    
    # Process each unique Genofeature that matches our selection
    for genofeature in df['Genofeature'].unique():
        if genofeature.lower() in selected:
            # Create directory for this genofeature
            os.makedirs(genofeature, exist_ok=True)
            
            # Filter records for this genofeature
            genofeature_df = df[df['Genofeature'] == genofeature]
            
            # Save the filtered TSV
            genofeature_df.to_csv(f"{genofeature}/filtered.tsv", sep='\\t', index=False)
            
            # Extract matching sequences
            matching_seqs = []
            for orf_id in genofeature_df['ORF_ID']:
                if orf_id in seq_dict:
                    matching_seqs.append(seq_dict[orf_id])
            
            # Write matching sequences to FASTA
            if matching_seqs:
                SeqIO.write(matching_seqs, f"{genofeature}/sequences.fasta", "fasta")
            
            print(f"Processed {genofeature}: {len(matching_seqs)} sequences")
    """
}

workflow {
    // Create input channels
    ch_filtered_tsv = Channel.fromPath(params.second_run)
    ch_input_fasta = Channel.fromPath(params.input_fasta)
    
    // Combine the input files into a tuple
    ch_inputs = ch_filtered_tsv.combine(ch_input_fasta)

    SPLIT_BY_GENOFEATURE(ch_inputs)
}