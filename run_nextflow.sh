#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=32000
#SBATCH --time=04:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

source /etc/profile.d/conda.sh
conda activate phidra
nextflow run main.nf
