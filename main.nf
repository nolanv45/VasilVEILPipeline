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

// process PASV {
//     conda "/home/nolanv/.conda/envs/phidra"

//     input:
//     tuple val(protein), path(fasta)

//     output:
//     path "pasv_output/", emit: "pasv_output"

//     script:
//     """
//     # Use current working directory
//     WORK_DIR=\$(pwd)

//     # Create directories
//     mkdir -p \$WORK_DIR/pasv_output/{input,output,pasv}

//     # Define mapping
//     if [ "${protein}" == "PolA" ]; then
//         mapped_name="POL_putative"
//     elif [ "${protein}" == "RNR" ]; then
//         mapped_name="RNR_putative"
//     else
//         echo "Error: Unknown protein ${protein}"
//         exit 1
//     fi

//     # Copy input file with mapped name
//     cp ${fasta} "\$WORK_DIR/pasv_output/input/\${mapped_name}.fasta"

//     # Run PASV
//     /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/pasv.sh \$WORK_DIR/pasv_output
//     """
// }


// process PASV {
//     conda "/home/nolanv/.conda/envs/phidra"
    
//     input:
//     tuple val(protein), path(fasta)

//     output:
//     tuple val(protein), path("pasv_output")

//     script:
//     """
//     #!/usr/bin/env bash
//     set -e  # Exit on error

//     # Setup directories
//     mkdir -p pasv_output/{input,output,pasv}

//     # Convert protein name to lowercase for comparison
//     protein_lower=\$(echo "${protein}" | tr '[:upper:]' '[:lower:]')

//     # Map protein names to PASV expected names and signatures
//     if [ "\${protein_lower}" == "pola" ]; then
//         mapped_name="POL_putative"
//     elif [ "\${protein_lower}" == "rnr" ]; then
//         mapped_name="RNR_putative"
//     else
//         echo "Error: Unknown protein ${protein}"
//         exit 1
//     fi

//     # Copy and prepare input files
//     cp ${fasta} "pasv_output/input/\${mapped_name}.fasta"
//     cp ${align_refs} "pasv_output/input/align_refs.fa"

//     #################################################################################
//     ### SCRIPT SETTINGS

//     ## SOFTWARE REQUIREMENTS

//     ## FILE LOCATIONS

//     # MARKER GENES
//     PROTEIN="${protein}"
//     #"RNR POL RecA UvsX RecB RecD PcrA UvrD Dda UvsW UvrB RecG SNF2 Gp4 Gp41 DnaB"


//     # reference directory for PASV alignment references
//     REF_DIR=/mnt/VEIL/references/pasv/pasv_profile
//     # output
//     OUTDIR="output"

//     ## OTHER PARAMETERS
//     ALIGNER=mafft
//     #clustalo, mafft

//     #PASV FILE LOCATIONS
//     PASV_DIR="pasv" 

//     #RNR
//     RNR_ALIGN=rnr__classIa_classII__best_practices.fa
//     RNR_ROI_START=437
//     RNR_ROI_END=625
//     RNR_CAT_SITES="437,439,441,462,438"
//     RNR_PASV_VER="\$2 ~/N/ && \$3 ~/C/ && \$4 ~/E/ && \$5 ~/C/ && \$10 ~/Both/"
//     RNR_NO_ROI="\$2 ~/N/ && \$3 ~/C/ && \$4 ~/E/ && \$5 ~/C/ "

//     #POL
//     POL_ALIGN=pola_16_k12_ref.fa
//     POL_ROI_START=521
//     POL_ROI_END=923
//     POL_CAT_SITES="668,705,758,762"
//     POL_PASV_VER="\$2 ~/R/ && \$3 ~/D/ && \$4 ~/K/ && \$9 ~/Both/"
//     POL_NO_ROI="\$2 ~/R/ && \$3 ~/D/ && \$4 ~/K/"

//     #################################################################################
//     ### SCRIPT BODY

//     ## COMMANDS TO RUN

 

//     # Make sure pasv is available
//     if [ ! -f pasv_output/pasv/pasv ]
//     then
//     cd pasv_output/pasv
//     wget https://github.com/mooreryan/pasv/releases/download/2.0.2/pasv-2.0.2-alpine-static.zip
//     unzip pasv_output/pasv/pasv-2.0.2-alpine-static.zip
//     chmod 755 pasv_output/pasv/pasv
//     ./pasv --help
//     cd pasv_output/output
//     else
//     fi

//     for p in \${PROTEIN}
//     do
//     INPUT_FILE=\${p}_INPUT
//     ALIGN=\${p}_ALIGN
//     ROI_START=\${p}_ROI_START
//     ROI_END=\${p}_ROI_END
//     CAT_SITES=\${p}_CAT_SITES
//     PRED_SITE=\${p}_PRED_SITE
//     PASV_VER=\${p}_PASV_VER
//     CTG_ID_COLUMN=\${p}_CTG_ID_COLUMN
//     NO_ROI=\${p}_NO_ROI

//     echo "Generate PASV signatures: \${p}"
//     pasv_output/pasv/pasv msa \
//     --outdir=pasv_output/output \
//     --force \
//     --roi-start=\${!ROI_START} \
//     --roi-end=\${!ROI_END} \
//     --jobs=20 \
//     --aligner=\${ALIGNER} \

    
//     pasv_output/input/\${!mapped_name}.fasta \
//     ${REF_DIR}/${!ALIGN} \
//     ${!CAT_SITES}

//     echo "Rename output file: \${p}"
//     mv pasv_output/output/\${!mapped_name}.pasv_signatures.tsv \
//     pasv_output/output/\${p}_putative.pasv_signatures.tsv
//     echo "renamed output file"

//    python /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/pasv_post.py \
//     pasv_output/output/\${p}_putative.pasv_signatures.tsv \
//     pasv_output/output/\${p}_pasv_boxplots.png \
//     \${p}



//     done

//     echo "=====JOB FINISH====="
//     """
// }

process PASV {
    conda "/home/nolanv/.conda/envs/phidra"
    cpus 4

    input:
    tuple val(protein), path(fasta), path(align_refs)

    output:
    tuple val(protein), path("pasv_output/output/${protein}_putative.pasv_signatures.tsv")

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
import seaborn as sns
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
          path(initial_search, stageAs: 'initial_search.tsv'), 
          path(recursive_search, stageAs: 'recursive_search.tsv')

    output:
    path "annotated_hits_${protein}.tsv"

    script:
    """
    #!/usr/bin/env python3
import sys
import os

protein = "${protein}"
file1_path = "initial_search.tsv"
file2_path = "recursive_search.tsv"
output_file = "annotated_hits_${protein}.tsv"

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
            def proteinConfig = params.proteins.find { it.protein.toLowerCase() == protein.toLowerCase() }
            def align_refs = file(proteinConfig.pasv_align_refs)
            tuple(protein, fasta, align_refs)
        }
        .view { "[DEBUG] To PASV: $it" }
        | PASV
    PASV.out
        .map { protein, pasv_signatures ->
            tuple(protein, pasv_signatures)
        }
        .view { "[DEBUG] PASV output: $it" }
        | PASV_post
}