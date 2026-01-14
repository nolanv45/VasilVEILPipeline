#  PHIDRA ![Version](https://img.shields.io/badge/version-2.0.0_stable-green)
**P**rotein **H**omology **I**dentification via **D**omain-**R**elated **A**rchitecture

A simple way to search and validate identified Pfam domains of interest against a curated InterProScan Domain Architecture (IDA) file to check whether or not your proteins match a domain composition found in the InterPro Database.

## Description

A Python-based package of scripts that performs an initial homology search using `MMseqs2` to identify targeted proteins of interest in small or large datasets. Top hits are searched against the Pfam database using `pfam_scan`, and verified domains are checked and compared against a custom input InterProScan Domain Architecture (IDA) file relative to your target protein of interest. 
A recursive search is then performed using the full-length proteins with validated IDAs as the subject database and the original input as the query, with initial matches filtered out. This process captures potentially more distant proteins that may have been missed in the initial homology search but are functionally relevant.

➡️  Official documentation is hosted on [PHIDRA](https://zschreib.github.io/tool-instructions/PHIDRA/).

## Getting started

### Dependencies

* `phidra` can be run either on a machine with a properly set up Python environment or by creating a custom [Conda](https://docs.conda.io/en/latest/) environment if you do not wish to change your current setup (recommended).
* [MMseqs2](https://github.com/soedinglab/MMseqs2) required for initial and recursive homology search. 
* [HMMER](http://hmmer.org/documentation.html) required for Pfam database creation and protein domain identification through `pfam_scan`.
* [Python](https://www.python.org/downloads/) >= 3.8 

### Installation

* Conda setup (recommended) [Help](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
```
#Create environment with >= Python 3.8
conda create --name phidra python=3.8

#Activate environment
conda activate phidra

#Install hmmer/mmseqs2 (Required)
conda install -c conda-forge -c bioconda mmseqs2
conda install bioconda::hmmer

#Clone tool into working directory
git clone https://github.com/zschreib/phidra

cd phidra

#Grabs required python packages
pip install -r requirements.txt

```
* Python setup >= Python 3.8, requires MMseqs2 and HMMER to be already set up.
```
git clone https://github.com/zschreib/phidra
cd phidra
pip install -r requirements.txt
```

* Tool and version check. If you receive any errors here, do not proceed until fixed.
```
mmseqs -h
hmmscan -h
python phidra_run.py -v
```
### Setting up Pfam Database (Optional)
* If you already have a custom or existing Pfam-HMM formatted database you can skip this step. 
* You can optionally include the **Pfam-B** database to perform a less restrictive search, which can help identify more remote or novel domain relationships that may be missed by the curated Pfam-A profiles.
```
mkdir pfam_database
cd pfam_database

wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

gunzip -c Pfam-A.hmm.dat.gz > Pfam-A.hmm.dat
gunzip -c Pfam-A.hmm.gz > Pfam-A.hmm
rm Pfam-A.hmm.gz Pfam-A.hmm.dat.gz

hmmpress Pfam-A.hmm
```


## Creating an InterProScan Domain Architecture (IDA) file

### 1. Identify essential domains
* For DNA polymerase A, I have identified PF00476 (DNA_pol_A) as the core domain.
* Search domain on InterPro: `https://www.ebi.ac.uk/interpro/entry/pfam/PF00476/domain_architecture/` and grab the full IDA profile.

### 2. Download domain architecture data in TSV format for all or selected domain combinations

File should include (by default):
* IDA ID (unique hash)
* Domain combinations
* Protein counts within IDA
* Representative sequence
* Representative length
* Domain positions/coordinates

### 3. Sample IDA file in TSV format:
```
IDA ID	IDA Text	Unique Proteins	Representative Accession	Representative Length	Representative Domains
82911c7e8cf0ed5121595d5944b2a2f9a2c4f49e	PF02739:IPR020046-PF01367:IPR020045-PF01612:IPR002562-PF00476:IPR001098	22766	P00582	928	PF02739{5_3_exonuc_N}:IPR020046{5-3_exonucl_a-hlix_arch_N}[9-170],PF01367{5_3_exonuc}:IPR020045{DNA_polI_H3TH}[171-272],PF01612{DNA_pol_A_exo1}:IPR002562{3'-5'_exonuclease_dom}[330-516],PF00476{DNA_pol_A}:IPR001098{DNA-dir_DNA_pol_A_palm_dom}[551-925]
d4033548d92ed83469d80faa39cea59a849060a8	PF00476:IPR001098	18435	P00581	704	PF00476{DNA_pol_A}:IPR001098{DNA-dir_DNA_pol_A_palm_dom}[333-701]
16983c3a921f1409678537d9977eca4ed176e7c2	PF02739:IPR020046-PF01367:IPR020045-PF22619:IPR054690-PF00476:IPR001098	10607	Q04957	877	PF02739{5_3_exonuc_N}:IPR020046{5-3_exonucl_a-hlix_arch_N}[4-170],PF01367{5_3_exonuc}:IPR020045{DNA_polI_H3TH}[172-266],PF22619{DNA_polI_exo1}:IPR054690{DNA_polI_exonuclease}[318-457],PF00476{DNA_pol_A}:IPR001098{DNA-dir_DNA_pol_A_palm_dom}[498-875]
ce54c2b3dc9b223aa1f4992353c0972accb51c74	PF02739:IPR020046-PF01367:IPR020045-PF00476:IPR001098	6038	O84500	866	PF02739{5_3_exonuc_N}:IPR020046{5-3_exonucl_a-hlix_arch_N}[3-166],PF01367{5_3_exonuc}:IPR020045{DNA_polI_H3TH}[167-259],PF00476{DNA_pol_A}:IPR001098{DNA-dir_DNA_pol_A_palm_dom}[495-865]
f54dc0def957e931720f1460942192debd6092ca	PF01612:IPR002562-PF00476:IPR001098	4576	Q05254	595	PF01612{DNA_pol_A_exo1}:IPR002562{3'-5'_exonuclease_dom}[19-210],PF00476{DNA_pol_A}:IPR001098{DNA-dir_DNA_pol_A_palm_dom}[243-587]
```
### 4. Analyze architectures

Core domain (PF00476) commonly associates with:
* PF02739 (5' exonuclease N-terminal)
* PF01367 (5' exonuclease)
* PF01612 (DNA polymerase A exonuclease)

### 5. Use for validation

* Compare unknown or known sequence hits to the subject database against known IDA profiles.
* Explore domain organization and composition.

Domain architectures help validate protein predictions by ensuring essential functional domains are present, and they also provide structured, functionally meaningful features that can be leveraged as high-quality inputs for machine-learning models to improve training and downstream classification.

## Running the tool
```
usage: phidra_run.py [-h] [-v] -i INPUT_FASTA -db SUBJECT_DB -pfam PFAM_HMM_DB -ida IDA_FILE -f FUNCTION -o OUTPUT_DIR
                     [-t THREADS] [-e EVALUE]

Identifies homologous proteins and associated Pfam domains from input protein sequences, while comparing against
InterPro Domain Architectures to analyze domain-level similarities and functional relationships.

Help:
  -h, --help            Show this help message and exit
  -v, --version         Show program version and exit

Required arguments:
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        Query FASTA for mmseqs search (default: None)
  -db SUBJECT_DB, --subject_db SUBJECT_DB
                        Subject FASTA for mmseqs createdb (default: None)
  -pfam PFAM_HMM_DB, --pfam_hmm_db PFAM_HMM_DB
                        Pfam HMM format database path (default: None)
  -ida IDA_FILE, --ida_file IDA_FILE
                        IDA TSV file (default: None)
  -f FUNCTION, --function FUNCTION
                        User label for this run (default: None)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Base output directory (default: None)

Optional arguments:
  -t THREADS, --threads THREADS
                        Threads for tools supporting -cpu/--threads (default: 1)
  -e EVALUE, --evalue EVALUE
                        E-value threshold for mmseqs easy-search (e.g., 1E-3, 1e-5) (default: 1E-3)
```

Example run using the provided examples directory:

```
python phidra_run.py -i examples/query/polA_test.fa -db examples/subject/pola_16_k12_ref.fa -pfam [pfam_DB_location] -ida examples/IDA/polA_IDA_list.tsv -f examples/output -t 5 -e 1E-3
```

## Results output summary
```
output/
├── final_results/
│   ├── pfam_coverage_report.tsv                # Summary of Pfam domain coverage across all search hits
│   ├── summary.tsv                             # Summary table combining counts for each iteration
│   ├── unvalidated_ida_pfams/               # Pfam domains found but not validated by iterative domain architecture (IDA)
│   │   ├── domains.fa                          # FASTA of individual Pfam domain hits
│   │   ├── full_proteins.fa                    # Full-length proteins containing those domains
│   │   └── pfam_unvalidated_merged_report.tsv  # Merged table of unvalidated domain results
│   └── validated_ida_pfams/                 # Pfam domains validated by IDA recursion
│       ├── domains.fa                          # FASTA of validated Pfam domain hits
│       ├── full_proteins.fa                    # Full-length proteins containing validated domains
│       └── pfam_validated_merged_report.tsv    # Merged table of validated domain results
│
├── mmseqs/
│   ├── initial/                             # First-pass MMseqs2 search output
│   │   ├── bits.tsv                             # Top hit table by best bitscore
│   │   ├── hits.fa                              # FASTA of sequences with best e-value
│   │   ├── hits.tsv                             # Top hit table by best e-value
│   │   └── res.m8                               # MMseqs2 table output (m8 format) for all significant initial hits
│   └── recursive/                           # Secondary MMseqs2 search using hits as queries
│       └── res.m8                               # MMseqs2 table output (m8 format) for all significant recursive hits
│                                                # No recursive hits so tophit/FASTA not created
└── pfam/
    ├── initial/                             # Initial Pfam HMMER domain search results
    │   ├── pfam_coverage_report.tsv             # Coverage summary of initial Pfam search
    │   ├── unvalidated_ida_report/          # Domains not validated by user IDA but hit Pfam domain
    │   │   ├── domains.fa                       # FASTA of unvalidated Pfam domain hits
    │   │   ├── full_proteins.fa                 # Full-length proteins for unvalidated hits
    │   │   └── pfam_unvalidated_report.tsv      # Tabular report of unvalidated domain results
    │   └── validated_ida_report/            # Domains validated by IDA
    │       ├── domains.fa                       # FASTA of validated Pfam domain hits
    │       ├── full_proteins.fa                 # Full-length proteins for validated hits
    │       └── pfam_validated_report.tsv        # Tabular report of validated domain results
    └── recursive/                               # Results from recursive Pfam search
        (empty)                                  # No recursive hits so pfam files not created
```

* The named output directory contains the complete output of the `phidra` pipeline.
* Includes homology search results, IDA-validated and unvalidated Pfam domain calls with both sequence-level and domain-level FASTA files, and comprehensive summary tables. 
* A detailed description of each file is provided above to help navigate and interpret the results.

## Looking for Version 1?

- Browse the stable **[v1.x branch](https://github.com/zschreib/phidra/tree/v1.x)**  
- Or check the **[v1.0.0 release](https://github.com/zschreib/phidra/releases/tag/v1.0.0)**  

## Citation
If you found this tool useful, please cite: 

Schreiber, Z. D. PHIDRA: *Protein Homology Identification via Domain-Related Architecture*. GitHub. Retrieved Sept 27, 2025, from [https://github.com/zschreib/phidra/](https://github.com/zschreib/phidra/).

## Authors

Contributor’s name and contact info

zschreib.dev@gmail.com

## License

This project is licensed under the GNU General Public License v3.0 see the LICENSE file for more details.

## Acknowledgments

Inspiration, code snippets, etc.
* [MMseqs2](https://github.com/soedinglab/MMseqs2)
* [InterPro](https://www.ebi.ac.uk/interpro/)
* [pfam_scan](https://github.com/aziele/pfam_scan)
