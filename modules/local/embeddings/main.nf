process EMBEDDINGS {
    publishDir "${params.outdir}",
        mode: 'copy'
        
    label 'process_gpu' 
    conda "/mnt/biostore-all/Polson/users/nolanv/pipeline_project/VasilVEILPipeline/containers/umap/umap.yml" 
    // container "containers/umap/umap.sif"
    
    input:
        tuple val(dataset), val(protein), val(genofeature), path(fasta)
        path model_folder
        
    output:
        path "embeddings/${dataset}/${protein}/${genofeature}", emit: embeddings_dirs
        path "versions.yml", emit: versions
        

    script:
    """
    set -euo pipefail

    out_dir="embeddings/${dataset}/${protein}/${genofeature}"
    mkdir -p "\$out_dir"

    python3 ${model_folder}/extract.py \
        ${model_folder}/esm2_t36_3B_UR50D.pt \
        "${fasta}" \
        "\$out_dir" \
        --repr_layers 36 \
        --include mean

    for batch_dir in "\$out_dir"/batch_*; do
        [ -d "\$batch_dir" ] || continue
        mv "\$batch_dir"/*.pt "\$out_dir"/
        rmdir "\$batch_dir" || true
    done

    python3 - <<'PY'
import platform
from importlib import import_module

process_name = "${task.process}"
versions = [
    f'"{process_name}":',
    f'    python: {platform.python_version()}'
]

for package_name in ['numpy', 'pandas', 'torch', 'umap']:
    try:
        module = import_module(package_name)
        version = getattr(module, '__version__', None)
        if version:
            versions.append(f'    {package_name}: {version}')
    except Exception:
        pass

with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write('\\n'.join(versions) + '\\n')
PY

    """
}