# VEIL Environmental Metagenomic Phage Replication Module Pipeline
A Nextflow pipeline for identifying and analyzing phage replication proteins and replication modules in environmental viral metagenomes.

The Viral Ecology and Informatics Lab (VEIL) at the University of Delaware investigates viral communities by identifying viral populations through the replication proteins encoded in their genomes. These replication proteins are critical to the biology of phages (viruses that infect microbes), influencing phenotypic characteristics of infection, such as replication speed, burst size, and infection strategy (virulent vs. temperate), which ultimately impact microbial host communities and nutrient cycling. 

This repository hosts an automated pipeline built with Nextflow which identifies and characterizes viral replication proteins and modules, specifically targeting DNA Polymerase A (PolA), DNA Polymerase B (PolB), Ribonucleotide Reductase (RNR), and helicases using user-established reference resources. The pipeline provides annotated outputs and visualizations.

## Features

- Automated detection of phage replication proteins (PolA, PolB, RNR, and helicases)
- Annotation of genofeature metadata through PHIDRA and PASV
- Replication module detection (co-occurring proteins on contigs)
- Protein embeddings via an ESM2 Protein Language Model (PLM) fine-tuned on diverse viral proteomes (Sawhney et al., 2025), which captures sequence, structural, and functional information specific to viral proteins
- Outputs include annotated tsvs, summary statistics, and UMAP projection visualizations
- Built with Nextflow for portability and HPC/cloud compatibility

## Acknowledgements

- This material is supported by the National Science Foundation under Grant No. 1736030.
- Use of the BioMix cluster was supported by the Delaware INBRE program, with a grant from the National Institute of General Medical Sciences - NIGMS (P20 GM103446) from the National Institutes of Health.
- This work is supported by the Undergraduate Research Program, as part of the Summer Scholars program at the University of Delaware.

## References

- Nasko, D. J., Chopyk, J., Sakowski, E. G., Ferrell, B. D., Polson, S. W., & Wommack, K. E. (2018). Family A DNA Polymerase Phylogeny Uncovers Diversity and Replication Gene Organization in the Virioplankton. Frontiers in microbiology, 9, 3053. https://doi.org/10.3389/fmicb.2018.03053
- Sawhney, R., Ferrell, B., Dejean, T., Schreiber, Z., Harrigan, W., Polson, S. W., Wommack, K. E., & Belcaid, M. (2025). Fine-Tuning Protein Language Models Unlocks the Potential of Underrepresented Viral Proteomes. bioRxiv : the preprint server for biology, 2025.04.17.649224. https://doi.org/10.1101/2025.04.17.649224
- https://github.com/zschreib/phidra
- Moore, R. M., Harrison, A. O., Nasko, D. J., Chopyk, J., Cebeci, M., Ferrell, B. D., Polson, S. W., & Wommack, K. E. (2021). PASV: Automatic protein partitioning and validation using conserved residues. bioRxiv. https://doi.org/10.1101/2021.01.20.427478
- https://zschreib.github.io/tool-instructions/
