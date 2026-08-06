process GENERATE_COORDINATES {
    tag "nn${nn}_md${md}"
    publishDir "${params.outdir}",
        mode: 'copy'
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pytorch_numpy_pandas_pruned:02ddcbc0a6f7a925' :
        'community.wave.seqera.io/library/python_pytorch_numpy_pandas_pruned:80951f99909b8c30' }"

    input:
        tuple path(embeddings_dirs), val(excluded_genofeatures), val(nn), val(md), val(publish_subdir)
        
    output:
        path "${publish_subdir}", emit: coordinates_files
        path "versions.yml", emit: versions

    script:
"""
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMBA_DEFAULT_NUM_THREADS=1
export NUMBA_THREADING_LAYER=workqueue
export PYTHONHASHSEED=42

python3 <<'PY'

import os

os.environ["NUMBA_CACHE_DIR"] = os.path.join(
    os.getcwd(),
    ".numba_cache"
)

os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)

import random
import hashlib
from pathlib import Path

import torch
import numpy as np
import pandas as pd
import umap



# -------------------------
# Seeds
# -------------------------

random.seed(42)
np.random.seed(42)
torch.manual_seed(42)


# -------------------------
# Load embeddings
# -------------------------

def load_embedding(file_path):
    try:
        embedding_data = torch.load(
            file_path,
            map_location='cpu'
        )

        embedding = (
            embedding_data["mean_representations"][36]
            .numpy()
            .astype(np.float32)
        )

        embedding_id = embedding_data.get("label")

        return embedding, embedding_id

    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None, None


def find_pt_files(base_dir):

    pt_files = []

    for root, dirs, files in os.walk(base_dir):

        dirs.sort()
        files.sort()

        for file in files:

            if file.endswith(".pt"):
                pt_files.append(
                    os.path.join(root, file)
                )

    return pt_files



embeddings = []
embedding_ids = []


base_dirs = [
${embeddings_dirs.collect { "'${it}'" }.join(',')}
]


pt_files = []

for base_dir in base_dirs:

    if os.path.isdir(base_dir):

        pt_files.extend(
            find_pt_files(base_dir)
        )


print(f"Found {len(pt_files)} .pt files")



# -------------------------
# Read embeddings
# -------------------------

excluded_raw = "${excluded_genofeatures}".strip()

excluded_subdirs = []

if excluded_raw and excluded_raw.lower() not in ['null','none']:

    excluded_subdirs = [
        x.strip()
        for x in excluded_raw.split(',')
    ]


for file_path in pt_files:

    genofeature = Path(file_path).parent.parent.name


    if genofeature.lower() not in [
        x.lower()
        for x in excluded_subdirs
    ]:

        embedding, embedding_id = load_embedding(file_path)

        if embedding is not None:

            embeddings.append(embedding)
            embedding_ids.append(embedding_id)



print(
    f"Processing {len(embeddings)} embeddings"
)



# -------------------------
# Sort deterministically
# -------------------------

sorted_pairs = sorted(
    zip(
        embedding_ids,
        embeddings
    ),
    key=lambda x: x[0]
)


embedding_ids, embeddings = zip(
    *sorted_pairs
)


embedding_ids = list(embedding_ids)

embeddings = np.stack(
    embeddings
).astype(np.float32)



print("Embedding ID hash:",
      hashlib.sha256(
          "\\n".join(embedding_ids).encode()
      ).hexdigest())


print("Embedding matrix hash:",
      hashlib.sha256(
          embeddings.tobytes()
      ).hexdigest())



# -------------------------
# UMAP
# -------------------------

nn = int(${nn})
md = float(${md})


print(
    f"Generating UMAP nn={nn}, md={md}"
)


reducer = umap.UMAP(

    n_components=2,

    n_neighbors=nn,

    min_dist=md,

    metric="cosine",

    random_state=42,

    transform_seed=42,

    n_epochs=200,

    low_memory=False,

    force_approximation_algorithm=False,

    angular_rp_forest=False

)



reducer.fit(embeddings)


coordinates = reducer.embedding_

print(
    "Coordinate hash:",
    hashlib.sha256(
        coordinates.tobytes()
    ).hexdigest()
)

# -------------------------
# Save coordinates
# -------------------------

md_int = int(md * 10)

coord_df = pd.DataFrame({
    "embedding_id": embedding_ids,
    "x": coordinates[:, 0],
    "y": coordinates[:, 1]
})

# Build connection table
groups = {}

for i, eid in enumerate(embedding_ids):
    contig_id = "_".join(eid.split("_")[:-3])

    if contig_id not in groups:
        groups[contig_id] = []

    groups[contig_id].append(i)

connections = []

for contig_id, indices in groups.items():
    if len(indices) > 1:
        for i in range(len(indices) - 1):
            for j in range(i + 1, len(indices)):
                connections.append([
                    embedding_ids[indices[i]],
                    embedding_ids[indices[j]]
                ])

connections_df = pd.DataFrame(
    connections,
    columns=["id1", "id2"]
)

# Create output directory
os.makedirs("${publish_subdir}", exist_ok=True)

coord_df.to_csv(
    "${publish_subdir}/coordinates_nn{}_md{}.tsv".format(
        nn,
        md_int
    ),
    sep="\\t",
    index=False
)

connections_df.to_csv(
    "${publish_subdir}/connections.tsv",
    sep="\\t",
    index=False
)

import platform

with open("versions.yml", "w") as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f"    python: {platform.python_version()}\\n")
    handle.write(f"    numpy: {np.__version__}\\n")
    handle.write(f"    pandas: {pd.__version__}\\n")
    handle.write(f"    torch: {torch.__version__}\\n")
    handle.write(f"    umap: {umap.__version__}\\n")
PY
"""
}