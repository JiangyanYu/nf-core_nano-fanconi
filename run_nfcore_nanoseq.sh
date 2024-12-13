#!/bin/bash

#SBATCH --job-name=nanoseq
# request short-gpu for the task
#SBATCH --partition=short-gpu
#SBATCH --cpus-per-task=4

# Load required modules
module load Nextflow
module load Docker  # Or Docker if you use Docker

# Run the nf-core/nanoseq pipeline
nextflow run nf-core/nanoseq -profile test,docker
