# VEIL Nextflow Pipeline Tutorial
This guide walks you through installing dependencies, downloading containers, and running the VEIL Nextflow pipeline on HPC or local systems.

# Initial setup

## 1. Install Nextflow
Install Nextflow according to latest installation guide. This can be done through the self-install or through a new conda environment. [https://www.nextflow.io/docs/latest/install.html](https://www.nextflow.io/docs/latest/install.html) 

## 2. Clone the Repository
Change directory to desired download location.

git clone https://github.com/nolanv45/VasilVEILPipeline.git

## 3. Download Pipeline Dependencies

For VEIL members with access to Biomix, simply copy dependencies into your cloned pipeline directory with the commands below.

cd VasilVEILPipeline

mkdir -p tools

cp -r /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/tools/* tools/

mkdir -p containers

cp -r /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/containers/* containers/

If not a VEIL member, perform steps below:

chmod +x download_dependencies.sh

The dependency folder contains large files, running the comand below could take over an hour, it is suggested to use a job scheduler.

./download_dependencies.sh

This script downloads:
Singularity container files (*.sif)
ESM Viral Fine-tuned Model (*.pt)
All files will automatically be placed in necessary directories.


# Configuration setup

## Output directory
outdir = "/path/to/output/directory"

## Dataset input selection
Edit this block to this example format:

Can enter as many datasets as needed.

params {
  datasets = [
    EX1: "/absolute/path/to/EX1",
    EX2: "/absolute/path/to/EX2
  ]
}

## Phidra and Pfam database paths
Copy absolute path of PHIDRA directory and pfamDB (found inside PHIDRA directory)

Example:
phidra_dir = "/path/to/phidra"
pfamDB     = "/path/to/phidra/pfam_database"

## Last method of genofeature identification
By default, proteins are sorted into three different tracks of final genofeature identification.
PASV proteins are assigned based on active sites

pasv_proteins = ["RNR", "PolA"] \
tophit_proteins = ["helicase"] \
domain_proteins = ["PolB"] 

## PASV signature configuration
VEIL standard PASV settings for both PolA and RNR are adjustable if needed. Expected signatures (sigs) will always be in output if found in input data. If any signatures make up more than 1% of the total protein data, but aren't in expected_sigs, will be included in the output. This 1% threshed can be modified with pasv_threshold, in a decimal format.

  pasv_threshold = 0.01

  pasv_settings = [
    PolA: [
      mapped_name: 'PolA_putative',
      roi_start: 521,
      roi_end: 923,
      cat_sites: "668,705,758,762",
      expected_sigs: ['RDKY','RDKF','RDKL']
    ],
    RNR: [
      mapped_name: 'RNR_putative',
      roi_start: 437,
      roi_end: 625,
      cat_sites: "437,439,441,462,438",
      expected_sigs: ['NCECL','NCECV']
    ]
  ]

## Domain genofeature sorted proteins
Proteins that are assigned genofeatures by matching domains can be put with their Pfam IDs.

pfam_annotation_map = [
  "rPolB" : "PF00136",
  "pPolB" : "PF03175"
]

## Protein configuration
Each protein block controls PHIDRA and PASV behavior. Enter each protein as a block separated by a comma in this format. If the protein will not use PASV, can omit pasv_align_refs.

proteins = [
    [ protein: "RNR", 
    pfamDomain: "/path/to/pfamDomain.tsv", 
    subjectDB: "/path/to/subjectDB.fa",
    pasv_align_refs: "/path/to/pasv_refs.fa"
    ],     
    [ protein: "PolB", 
    pfamDomain: "/path/to/pfamDomain.tsv", 
    subjectDB: "/path/to/subjectDB.fa",
    ]
]

## UMAP and HDBSCAN parameters
By default, several UMAP and HDBSCAN parameters are used for UMAP and HDBSCAN tiling. Additional values can be added. These include minimum distance (md), nearest neighbors (nn), and minimum cluster (mc) values. The tiled output can help decide which final md and nn values to use.


md = [0.3, 0.5, 0.7]
nn = [75, 100]
mc = [10, 30]

## Final analysis parameters
Leave the below parameter alone and keep it false.

final_analysis = false

If there are any datasets that should be removed for the final analysis and plotting, enter them here:

exclude_datasets = ['ENA', 'BATS']

If there are any genofeatures that should be removed for the final analysis and plotting, enter them here:

excluded_genofeatures = ['pPolB', 'PolA']

The last two parameters decide which md and nn values to use for the genofeature centric UMAP plot analysis.

final_md = 0.7
final_nn = 100

# Running the pipeline

## First run
When the pipeline is run in a clean designated output directory, it will stop after generating all relevant UMAP and HDBSCAN plots.

With the configuration setup complete, ensure you are in the pipeline directory and run:

nextflow run main.nf

This command can be implemented into any job scheduler script as well.

## Second run
When the pipeline is run for a second time with existing embeddings in the output directory, it will skip the protein/genofeature identification steps and immediately generate UMAP and HDBSCAN plots, since embeddings are the most resource intense step.

## Final analysis
The genofeature centric plots can be generated by adding a flag to the run comand.

nextflow run main.nf --final_analysis

The pipeline will then regenerate coordinates if necessary, such as excluding datasets and/or genofeatures, and produce genofeature centric plots.