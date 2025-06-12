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

    #Add if additional proteins
    # RecA
    RecA_ALIGN=RecA_refs_2023-07-24.fasta
    RecA_ROI_START=10
    RecA_ROI_END=223
    RecA_CAT_SITES="69,70,75,195,215,216"
    RecA_PASV_VER="\$8 ~/ESTQAL/ && \$11 ~/Both/"
    RecA_NO_ROI="\$8 ~/ESTQAL/"

    # UvsX
    UvsX_ALIGN=UvsX_refs_2023-06-27.fasta
    UvsX_ROI_START=23
    UvsX_ROI_END=196
    UvsX_CAT_SITES="23,59,61,62,63,65,66,149,150,195"
    UvsX_PASV_VER="\$12 ~/SAPSKFKASH/ && \$15 ~/Both/"
    UvsX_NO_ROI="\$12 ~/SAPSKFKASH/"

    # RecB
    RecB_ALIGN=RecB_refs_2023-06-27.fasta
    RecB_ROI_START=23
    RecB_ROI_END=1082
    RecB_CAT_SITES="24,1020,1067,1082"
    RecB_PASV_VER="\$6 ~/SEDK/ && \$9 ~/Both/"
    RecB_NO_ROI="\$6 ~/SEDK/"

    # RecD
    RecD_ALIGN=RecD_refs_2023-06-27.fasta
    RecD_ROI_START=154
    RecD_ROI_END=495
    RecD_CAT_SITES="154,241,272,487,495"
    RecD_PASV_VER="\$7 ~/QHSNN/ && \$10 ~/Both/"
    RecD_NO_ROI="\$7 ~/QHSNN/"

    # PcrA
    PcrA_ALIGN=PcrA_refs_2023-07-20.fasta
    PcrA_ROI_START=12
    PcrA_ROI_END=570
    PcrA_CAT_SITES="223,224,226,227,254,286,493,497,504"
    PcrA_PASV_VER="\$11 ~/DEQDQYYLE/ && \$14 ~/Both/"
    PcrA_NO_ROI="\$11 ~/DEQDQYYLE/"

    # UvrD
    UvrD_ALIGN=UvrD_refs_2023-07-20.fasta
    UvrD_ROI_START=21
    UvrD_ROI_END=652
    UvrD_CAT_SITES="238,239,241,242,276,486,611,652"
    UvrD_PASV_VER="\$10 ~/DEQDQIVS/ && \$13 ~/Both/"
    UvrD_NO_ROI="\$10 ~/DEQDQIVS/"

    # Dda
    Dda_ALIGN=Dda_refs_2023-07-20.fasta
    Dda_ROI_START=9
    Dda_ROI_END=436
    Dda_CAT_SITES="115,116,118,121,146,252,402,432"
    Dda_PASV_VER="\$10 ~/DESDQISG/ && \$13 ~/Both/"
    Dda_NO_ROI="\$10 ~/DESDQISG/"

    # UvsW
    UvsW_ALIGN=UvsW_refs_2023-07-20.fasta
    UvsW_ROI_START=114
    UvsW_ROI_END=470
    UvsW_CAT_SITES="114,211,234,409,410"
    UvsW_PASV_VER="\$7 ~/PTCST/ && \$10 ~/Both/"
    UvsW_NO_ROI="\$7 ~/PTCST/"

    # UvrB
    UvrB_ALIGN=UvrB_refs_2023-06-27.fasta
    UvrB_ROI_START=26
    UvrB_ROI_END=415
    UvrB_CAT_SITES="38,39,41,42,141,143"
    UvrB_PASV_VER="\$8 ~/LGTGSS/ && \$11 ~/Both/"
    UvrB_NO_ROI="\$8 ~/LGTGSS/"

    # RecG
    RecG_ALIGN=RecG_refs_2023-06-30.fasta
    RecG_ROI_START=283
    RecG_ROI_END=448
    RecG_CAT_SITES="295,297,298,299,300,332,426"
    RecG_PASV_VER="\$9 ~/QDVGSQM/ && \$12 ~/Both/"
    RecG_NO_ROI="\$9 ~/QDVGSQM/"

    # SNF2
    SNF2_ALIGN=SNF2_refs_2023-10-17.fasta
    SNF2_ROI_START=81
    SNF2_ROI_END=546
    SNF2_CAT_SITES="81,109,110,197,198,538,542"
    SNF2_PASV_VER="\$9 ~/QGKDEQR/ && \$12 ~/Both/"
    SNF2_NO_ROI="\$9 ~/QGKDEQR/"

    # Gp4
    Gp4_ALIGN=Gp4_refs_2023-06-27.fasta
    Gp4_ROI_START=151
    Gp4_ROI_END=548
    Gp4_CAT_SITES="158,316,378,437"
    Gp4_PASV_VER="\$6 ~/GMFD/ && \$9 ~/Both/"
    Gp4_NO_ROI="\$6 ~/GMFD/"

    # Gp41
    Gp41_ALIGN=Gp41_refs_2023-06-27.fasta
    Gp41_ROI_START=62
    Gp41_ROI_END=445
    Gp41_CAT_SITES="114,220,234,299"
    Gp41_PASV_VER="\$6 ~/QNAL/ && \$9 ~/Both/"
    Gp41_NO_ROI="\$6 ~/QNAL/"

    # DnaB
    DnaB_ALIGN=DnaB_refs_2023-06-27.fasta
    DnaB_ROI_START=200
    DnaB_ROI_END=467
    DnaB_CAT_SITES="232,263,316,390"
    DnaB_PASV_VER="\$6 ~/RMSE/ && \$9 ~/Both/"
    DnaB_NO_ROI="\$6 ~/RMSE/"

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

   #  echo "Filter putative proteins to those with cat sites and roi: ${p}"
   #  awk 'NR==1 {print }; \
   #  '"${!PASV_VER}"' {print }' \
   #  ${BASEDIR}/${OUTDIR}/${p}_putative.pasv_signatures.tsv \
   #  > ${BASEDIR}/${OUTDIR}/${p}_validated.pasv_signatures.tsv


   python /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/pasv_post.py \
    ${BASEDIR}/${OUTDIR}/${p}_putative.pasv_signatures.tsv \
    ${BASEDIR}/${OUTDIR}/${p}_pasv_boxplots.png \
    ${p}


   #  cut -d$'\t' -f 1 ${BASEDIR}/${OUTDIR}/${p}_validated.pasv_signatures.tsv | \
   #  sort -u > ${BASEDIR}/${OUTDIR}/${p}_validated_orfids.txt

   #  seqkit grep -f ${BASEDIR}/${OUTDIR}/${p}_validated_orfids.txt \
   #  ${BASEDIR}/${INPUT_DIR}/${!INPUT_FILE}.${INPUT_EXT} \
   #  -o ${BASEDIR}/${OUTDIR}/${p}_validated.${INPUT_EXT}

   #  echo "Filter putative proteins to those with cat sites: ${p}"
   #  awk 'NR==1 {print }; \
   #  '"${!NO_ROI}"' {print }' \
   #  ${BASEDIR}/${OUTDIR}/${p}_putative.pasv_signatures.tsv \
   #  > ${BASEDIR}/${OUTDIR}/${p}_cat_sites_regardless_roi.pasv_signatures.tsv

   #  if [ ${p}="RNR" ]
   #  then
   #  echo "Count number of proteins with 438 amino acid"
   #  awk '{count[$6]++} \
   #  END {for (word in count) print word, count[word]}' \
   #  ${BASEDIR}/${OUTDIR}/${p}_validated.pasv_signatures.tsv \
   #  | tee -a ${BASEDIR}/${OUTDIR}/${p}_stats.tsv
   #  fi

   #  if [ ${p}="POL" ]
   #  then
   #  echo "Count number of proteins with 762 amino acid"
   #  awk '{count[$5]++} \
   #  END {for (word in count) print word, count[word]}' \
   #  ${BASEDIR}/${OUTDIR}/${p}_validated.pasv_signatures.tsv \
   #  | tee -a ${BASEDIR}/${OUTDIR}/${p}_stats.tsv
   #  fi

   #  echo "Count putative and valid sequences: ${p}"
   #  seqkit stats -T -a ${BASEDIR}/${INPUT_DIR}/${!INPUT_FILE}.${INPUT_EXT} | \
   #  tee -a ${BASEDIR}/${OUTDIR}/${p}_stats.tsv
   #  seqkit stats -T -a ${BASEDIR}/${OUTDIR}/${p}_validated.${INPUT_EXT} | \
   #  tee -a ${BASEDIR}/${OUTDIR}/${p}_stats.tsv

   #  echo "done: ${p}"

    done

    echo "=====JOB FINISH====="
