 #################################################################################
    ### SCRIPT SETTINGS

    ## SOFTWARE REQUIREMENTS

    ## FILE LOCATIONS
    BASEDIR="$1"

    # MARKER GENES
    PROTEIN="POL"
    #"RNR POL RecA UvsX RecB RecD PcrA UvrD Dda UvsW UvrB RecG SNF2 Gp4 Gp41 DnaB"

    # input directory and file for query sequences
    INPUT_DIR="input"
    INPUT_EXT=fasta
    #INPUT_FILE=putative
    RNR_INPUT=RNR_putative
    POL_INPUT=POL_putative
    RecA_INPUT=HEL_putative
    UvsX_INPUT=HEL_putative
    RecB_INPUT=HEL_putative
    RecD_INPUT=HEL_putative
    PcrA_INPUT=HEL_putative
    UvrD_INPUT=HEL_putative
    Dda_INPUT=HEL_putative
    UvsW_INPUT=HEL_putative
    UvrB_INPUT=HEL_putative
    RecG_INPUT=HEL_putative
    SNF2_INPUT=HEL_putative
    Gp4_INPUT=HEL_putative
    Gp41_INPUT=HEL_putative
    DnaB_INPUT=HEL_putative

    # reference directory for PASV alignment references
    REF_DIR=/mnt/VEIL/references/pasv/pasv_profile
    # output
    OUTDIR="output"

    ## OTHER PARAMETERS
    ALIGNER=mafft
    #clustalo, mafft

    #PASV FILE LOCATIONS
    PASV_DIR="pasv" 

    #RNR
    RNR_ALIGN=rnr__classIa_classII__best_practices.fa
    RNR_ROI_START=437
    RNR_ROI_END=625
    RNR_CAT_SITES="437,439,441,462,438"
    RNR_PASV_VER="\$2 ~/N/ && \$3 ~/C/ && \$4 ~/E/ && \$5 ~/C/ && \$10 ~/Both/"
    RNR_NO_ROI="\$2 ~/N/ && \$3 ~/C/ && \$4 ~/E/ && \$5 ~/C/ "

    #POL
    POL_ALIGN=pola_16_k12_ref.fa
    POL_ROI_START=521
    POL_ROI_END=923
    POL_CAT_SITES="668,705,758,762"
    POL_PASV_VER="\$2 ~/R/ && \$3 ~/D/ && \$4 ~/K/ && \$9 ~/Both/"
    POL_NO_ROI="\$2 ~/R/ && \$3 ~/D/ && \$4 ~/K/"

    #################################################################################
    ### SCRIPT BODY

    ## JOB LOG HEADER - adds information about the start time, node name, and other useful information to the slurm output file
    perl -E 'say"="x120'; echo "JOB STARTED: `date`"; echo "NODE: `hostname`"; echo "SCRIPT EXECUTED: ${0}"; perl -E 'say"="x120'; echo ""

    ## COMMANDS TO RUN

    # outdir
    mkdir -p ${BASEDIR}/${OUTDIR}

    # Make sure pasv is available
    if [ ! -f ${BASEDIR}/${PASV_DIR}/pasv ]
    then
    cd ${BASEDIR}/${PASV_DIR}
    wget https://github.com/mooreryan/pasv/releases/download/2.0.2/pasv-2.0.2-alpine-static.zip
    unzip ${BASEDIR}/${PASV_DIR}/pasv-2.0.2-alpine-static.zip
    chmod 755 ${BASEDIR}/${PASV_DIR}/pasv
    ./pasv --help
    cd ${BASEDIR}/${OUTDIR}
    fi

    for p in ${PROTEIN}
    do
    INPUT_FILE=${p}_INPUT
    ALIGN=${p}_ALIGN
    ROI_START=${p}_ROI_START
    ROI_END=${p}_ROI_END
    CAT_SITES=${p}_CAT_SITES
    PRED_SITE=${p}_PRED_SITE
    PASV_VER=${p}_PASV_VER
    CTG_ID_COLUMN=${p}_CTG_ID_COLUMN
    NO_ROI=${p}_NO_ROI

    echo "Generate PASV signatures: ${p}"
    ${BASEDIR}/${PASV_DIR}/pasv msa \
    --outdir=${BASEDIR}/${OUTDIR} \
    --force \
    --roi-start=${!ROI_START} \
    --roi-end=${!ROI_END} \
    --jobs=20 \
    --aligner=${ALIGNER} \
    ${BASEDIR}/${INPUT_DIR}/${!INPUT_FILE}.${INPUT_EXT} \
    ${REF_DIR}/${!ALIGN} \
    ${!CAT_SITES}

    echo "Rename output file: ${p}"
    mv ${BASEDIR}/${OUTDIR}/${!INPUT_FILE}.pasv_signatures.tsv \
    ${BASEDIR}/${OUTDIR}/${p}_putative.pasv_signatures.tsv
    echo "renamed output file"

   python /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/pasv_post.py \
    ${BASEDIR}/${OUTDIR}/${p}_putative.pasv_signatures.tsv \
    ${BASEDIR}/${OUTDIR}/${p}_pasv_boxplots.png \
    ${p}



    done

    echo "=====JOB FINISH====="
