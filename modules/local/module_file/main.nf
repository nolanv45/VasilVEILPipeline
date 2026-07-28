process MODULE_FILE {
    label "process_single"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas:2a5ca2e7dd4ced9c' :
        'community.wave.seqera.io/library/python_pandas:eb68d9296e3f036a' }"
        
    publishDir "${params.outdir}/05_final_analysis/", 
        mode: 'copy',
        pattern: "*.tsv"

    input:
    path(tsv_file)
    path(metadata_file)

    output:
    path("rep_module_df.tsv"), emit: module_file
    path("versions.yml"), emit: versions

    script:
    def excluded_list = params.excluded_genofeatures.collect { feature -> "\'${feature}\'" }.join(', ')
    """
#!/usr/bin/env python3
# will take tsv files that output from each dataset per orf basis.
# will filter genofeatures based on metadata file only, so user input required
import pandas as pd
import glob
import os

df = pd.read_csv("${tsv_file}", sep='\t')
metadata = pd.read_csv("${metadata_file}", sep='\t')
excluded_features = [${excluded_list}]
genofeatures = []
for feature in metadata['genofeature']:
    if feature not in excluded_features:
        genofeatures.append(feature)

df = df[df['genofeature'].isin(genofeatures)]
contig_df = pd.DataFrame({'contig_id': df['contig_id'].unique()})

for genofeature in genofeatures:
    feature_data = df[df['genofeature'] == genofeature]
    if len(feature_data) > 0:
        feature_orfs = feature_data.set_index('contig_id')['orf_id']
        contig_df[genofeature] = contig_df['contig_id'].map(
            feature_orfs.groupby('contig_id').first()
        )
    else:
        contig_df[genofeature] = None
genofeature_columns = [col for col in contig_df.columns if col != 'contig_id']
contig_df['module'] = contig_df[genofeature_columns].apply(
    lambda row: '_'.join([col for col, val in zip(genofeature_columns, row) if pd.notna(val)]), 
    axis=1
)
contig_df['dataset'] = contig_df['contig_id'].map(
    df.groupby('contig_id')['dataset'].first()
)
cols = ['dataset', 'contig_id'] + genofeature_columns + ['module']
contig_df = contig_df[cols]

contig_df.to_csv("rep_module_df.tsv", sep='\t', index=False)

# stats file
module_counts = contig_df['module'].value_counts().reset_index()
module_counts.columns = ['module', 'count']
module_counts.to_csv("module_stats.tsv", sep='\t', index=False)

import platform
with open('versions.yml', 'w', encoding='utf-8') as handle:
    handle.write(f'"${task.process}":\\n')
    handle.write(f'    python: {platform.python_version()}\\n')
    handle.write(f'    pandas: {pd.__version__}\\n')
    """
}