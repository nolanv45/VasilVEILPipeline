#!/bin/bash
#SBATCH --job-name=ENA_pipeline
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-%j-%x-%N.out
#SBATCH --mem=80G 
#SBATCH --error=nextflow.err
#SBATCH --mail-user=nolanv@udel.edu

source /etc/profile.d/conda.sh
conda activate phidra
nextflow run main.nf
