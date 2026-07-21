process EMBEDDING_PLAN {
    label "process_single"
    conda "${moduleDir}/environment.yml"
    input:
    path(combined_tsv)

    output:
    path("**/*.fasta"), emit: planned_fastas, optional: true
    path("versions.yml"), emit: versions

    script:
    def datasetsJson = groovy.json.JsonOutput.toJson(params.datasets)
    """
    #!/usr/bin/env bash
    set -euo pipefail

    python3 - <<'PY'
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

combined_tsv = "${combined_tsv}"
datasets = json.loads('''${datasetsJson}''')

groups = defaultdict(set)

def normalize_id(value: str) -> str:
    token = (value or '').strip().split()[0]
    token = token.replace('.', '-')
    token = re.sub(r'[^A-Za-z0-9_-]', '', token)
    return token

with open(combined_tsv, 'r', encoding='utf-8', newline='') as combined_handle:
    reader = csv.DictReader(combined_handle, delimiter='\\t')
    fieldnames = set(reader.fieldnames or [])
    required_cols = {'dataset', 'protein', 'genofeature', 'orf_id'}
    missing_cols = required_cols - fieldnames
    if missing_cols:
        raise KeyError(f"Missing required columns in combined TSV: {sorted(missing_cols)}")

    for row in reader:
        dataset = (row.get('dataset') or '').strip()
        protein = (row.get('protein') or '').strip()
        genofeature = (row.get('genofeature') or '').strip()
        orf_id = normalize_id(row.get('orf_id') or '')
        if dataset and protein and genofeature and orf_id:
            groups[(dataset, protein, genofeature)].add(orf_id)

for (dataset, protein, genofeature), expected_ids in groups.items():

    fasta_path = Path(datasets.get(dataset, ''))
    if not str(fasta_path):
        raise KeyError(f"No fasta configured for dataset {dataset}")
    if not fasta_path.exists():
        raise FileNotFoundError(f"Configured fasta does not exist for {dataset}: {fasta_path}")

    output_root = Path("${params.outdir}") / "embeddings" / dataset / protein / genofeature

    existing_ids = {
        normalize_id(pt.stem)
        for pt in output_root.rglob('*.pt')
    } if output_root.exists() else set()

    # Fallback guard for slight naming differences between TSV IDs and pt stems.
    if len(existing_ids) >= len(expected_ids):
        continue

    missing_ids = sorted(expected_ids.difference(existing_ids))

    if not missing_ids:
        continue

    missing_set = set(missing_ids)
    out_fasta = Path(dataset) / protein / genofeature / f"{dataset}_{protein}_{genofeature}.fasta"
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    with fasta_path.open('r', encoding='utf-8') as fasta_handle, out_fasta.open('w', encoding='utf-8') as out_handle:
        header = None
        sequence_lines = []
        written = {'count': 0}

        def flush_record(record_header, record_lines):
            if record_header is None:
                return
            record_id = normalize_id(record_header)
            if record_id in missing_set:
                out_handle.write(f">{record_header}\\n")
                for seq_line in record_lines:
                    out_handle.write(seq_line)
                    if not seq_line.endswith('\\n'):
                        out_handle.write('\\n')
                written['count'] += 1

        for raw_line in fasta_handle:
            if raw_line.startswith('>'):
                flush_record(header, sequence_lines)
                header = raw_line[1:].strip()
                sequence_lines = []
            else:
                sequence_lines.append(raw_line)

        flush_record(header, sequence_lines)

    # Do not emit empty FASTA files; that would trigger unnecessary EMBEDDINGS tasks.
    if out_fasta.stat().st_size == 0 or written['count'] == 0:
        out_fasta.unlink(missing_ok=True)
PY

    python3 - <<'PY'
import platform

process_name = "${task.process}"
python_version = platform.python_version()

try:
    import pandas as pd
    pandas_version = pd.__version__
except Exception:
    pandas_version = 'not_installed'

with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"{process_name}":\\n')
    handle.write(f'    python: {python_version}\\n')
    handle.write(f'    pandas: {pandas_version}\\n')
PY
    """
}