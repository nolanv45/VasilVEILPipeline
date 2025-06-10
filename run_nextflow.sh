#!/bin/bash
#SBATCH --job-name=phidra
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-%j-%x-%N.out
#SBATCH --mem=50G
#SBATCH --time=20:00:00 
#SBATCH --error=nextflow.err

source /etc/profile.d/conda.sh
conda activate phidra
nextflow run main.nf
