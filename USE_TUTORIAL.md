# VEIL Nextflow Pipeline Tutorial

This guide walks you through installing dependencies, downloading containers, and running the VEIL Nextflow pipeline on HPC or local systems.

---

# Initial Setup

## 1. Install Nextflow

Install Nextflow according to latest installation guide. Install Nextflow onto your personal HPC profile if using one. This can be done through the self-install or through a new conda environment: [https://www.nextflow.io/docs/latest/install.html](https://www.nextflow.io/docs/latest/install.html)

## 2. Clone the Repository

Change directory to desired download location:

```bash
git clone https://github.com/nolanv45/VasilVEILPipeline.git
```

## 3. Download Pipeline Dependencies

### For VEIL members with access to Biomix

Simply copy dependencies into your cloned pipeline directory with the commands below:

```bash
cd VasilVEILPipeline
mkdir -p tools
cp -r /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/tools/* tools/
mkdir -p containers
cp -r /mnt/VEIL/users/nolanv/pipeline_project/VasilVEILPipeline/containers/* containers/
```

Then, perform these commands.

```bash
chmod +x pfam_download.sh
./pfam_download.sh
```

### For non-VEIL members

Perform steps below:

```bash
chmod +x download_dependencies.sh
./download_dependencies.sh
```

```bash
chmod +x pfam_download.sh
./pfam_download.sh
```

**Note:** The dependency folder contains large files. Running the command below could take over an hour, it is suggested to use a job scheduler.

This script downloads:
- Singularity container files (`*.sif`)
- ESM Viral Fine-tuned Model (`*.pt`)

All files will automatically be placed in necessary directories.

---

# Configuration Setup

All of the needed input data and parameters discussed below are located in the `nextflow.config` file inside the pipeline directory.

## Output Directory

The output directory can be any empty directory of choice:

```groovy
outdir = "/path/to/output/directory"
```

## Dataset Input Selection

Edit this block to this example format. Can enter as many datasets as needed. Naming convention is up to the user:

```groovy
params {
  datasets = [
    Dataset1: "/absolute/path/to/Dataset1",
    Dataset2: "/absolute/path/to/Dataset2"
  ]
}
```

Additionally, the underscore-delimited fields can be extraced for orf start/stop positions on their relative contigs. The number put into orf_coord_fields needs to be the start orf position. Write a negative index; -1 being last field, -2 second to last, etc.

```groovy
  orf_coord_fields = [BATS: -2, GOV: -3]
```

## PHIDRA and Pfam Database Paths

PHIDRA conducts homology search and Pfam domain-based validation to identify proteins of interest. The path of PHIDRA directory and pfamDB (found inside PHIDRA directory) is by default included in the pipeline directory.

Example:

```groovy
phidra_dir = "/path/to/PHIDRA"
pfamDB     = "/path/to/PHIDRA/pfam_database"
```

**Note:** The variable names `phidra_dir` and `pfamDB` must remain the same as they are referenced in the pipeline code.

## Genofeature Identification Methods

A **genofeature** (genomic feature) is the level of genomic information relevant to enzyme biochemistry and/or phenotype prediction. Genofeatures vary by protein - sometimes it's the protein itself, sometimes it is a Pfam domain, and sometimes it is even the amino acid at a specific active site.

By default, proteins are sorted into three different methods for genofeature identification:

- **PASV proteins** are assigned genofeatures based on active site signatures
- **Tophit proteins** are assigned genofeatures based on top BLAST hits
- **Domain proteins** are assigned genofeatures based on Pfam domain matches

All of the proteins go through PHIDRA validation, but the genofeature is assigned either through PASV, BLAST tophit, or its domain.

```groovy
pasv_proteins   = ["RNR", "PolA"]
tophit_proteins = ["helicase"]
domain_proteins = ["PolB"]
```

**Note:** Additional proteins can be added to the parameter above, and the naming convention must stay consistent throughout the rest of the configuration file.

## PASV Signature Configuration

PASV assigns genofeatures to proteins based on active site signatures. VEIL standard PASV settings for both PolA and RNR are adjustable if needed.

Expected signatures (sigs) will always be in output if found in input data. If any signatures make up more than 1% of the total protein data, but aren't in expected_sigs, will be included in the output. This 1% threshold can be modified with `pasv_threshold`, in a decimal format. This is by default included in the configuration file.

```groovy
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
```


## Domain Genofeature Sorted Proteins

Some proteins are assigned a genofeature based on Pfam domain identification. For these domain-sorted proteins, specify the mapping between genofeature names and their corresponding Pfam domain IDs. Use a comma-separated list format with genofeature name followed by Pfam domain ID.

```groovy
pfam_annotation_map = [
  "rPolB" : "PF00136",
  "pPolB" : "PF03175"
]
```

## Protein Configuration

Each protein block controls PHIDRA and PASV behavior. Enter each protein as a block separated by a comma in this format.

### Required parameters for each protein:

- `protein`: The name of the protein to analyze (keep consistent with the genofeature identification methods section, it is case sensitive)
- `pfamDomain`: Path to a TSV file containing InterPro Domain Architectures for use in PHIDRA
- `subjectDB`: Path to a FASTA file of known proteins used for homology search for use in PHIDRA

### PASV protein additional parameter:

- `pasv_align_refs`: Path to a FASTA file of reference sequences for PASV alignment and active site analysis. Only required for proteins using PASV genofeature assignment. (VEIL only uses PASV for RNR and PolA)

```groovy
proteins = [
  [ 
    protein: "RNR", 
    pfamDomain: "/path/to/pfamDomain.tsv", 
    subjectDB: "/path/to/subjectDB.fa",
    pasv_align_refs: "/path/to/pasv_refs.fa"
  ],     
  [ 
    protein: "PolB", 
    pfamDomain: "/path/to/pfamDomain.tsv", 
    subjectDB: "/path/to/subjectDB.fa"
  ]
]
```

## UMAP and HDBSCAN Parameters

By default, several UMAP and HDBSCAN parameters are used for UMAP and HDBSCAN tiling. Additional values can be added. These include:

- **md** - minimum distance
- **nn** - nearest neighbors  
- **mc** - minimum cluster

The tiled output can help decide which final md and nn values to use.

```groovy
md = [0.3, 0.5, 0.7]
nn = [75, 100]
mc = [10, 30]
```

## Final Analysis Parameters

Leave the below parameter alone and keep it false:

```groovy
final_analysis = false
```

If there are any datasets that should be removed for the final analysis and plotting, enter them here:

```groovy
exclude_datasets = ['ENA', 'BATS']
```

If there are any genofeatures that should be removed for the final analysis and plotting, enter them here:

```groovy
excluded_genofeatures = ['pPolB', 'PolA']
```

The last two parameters decide which md and nn values to use for the genofeature centric UMAP plot analysis:

```groovy
final_md = 0.7
final_nn = 100
```

---

# Running the Pipeline

## First Run

When the pipeline is run in a clean designated output directory, it will stop after generating all relevant UMAP and HDBSCAN plots.

With the configuration setup complete, ensure you are in the pipeline directory and run:

```bash
nextflow run main.nf
```

This command can be implemented into any job scheduler script as well.

## Second Run

When the pipeline is run for a second time with existing embeddings in the output directory, it will skip the protein/genofeature identification steps and immediately generate UMAP and HDBSCAN plots, since embeddings are the most resource intense step. This becomes time-saving when adding additional minimum distance or nearest neighbor values for the UMAP plot grid.

## Final Analysis

The genofeature centric plots can be generated by adding a flag to the run command:

```bash
nextflow run main.nf --final_analysis
```

The pipeline will then regenerate coordinates if necessary, if the config file excludes any datasets and/or genofeatures, and produce genofeature centric plots.