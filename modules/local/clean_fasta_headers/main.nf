process CLEAN_FASTA_HEADERS {
    tag "${meta.id}"
    label "process_single"
    conda "${moduleDir}/environment.yml"

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("${meta.id}_cleaned.fasta"), emit: fasta
        path("versions.yml"), emit: versions
    script:
    """
    awk 'BEGIN{RS=">"; ORS=""} NR>1{n=index(\$0, "\\n"); header=substr(\$0,1,n-1); sub(/[[:space:]].*/, "", header); gsub(/\\./, "-", header); gsub(/[^A-Za-z0-9_-]/, "", header); if (header == "") header = "seq_" NR; seq=substr(\$0,n+1); gsub(/\\n/, "", seq); gsub(/\\*/, "", seq); print ">"header"\\n"seq"\\n"}' ${fasta} > "${meta.id}_cleaned.fasta"

    python3 - <<'PY'
import subprocess

awk_version = subprocess.check_output("awk --version | head -n1 | sed 's/.*Awk //; s/,.*//'", shell=True, text=True).strip()

with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    awk: {awk_version}\\n')
PY
    """
}