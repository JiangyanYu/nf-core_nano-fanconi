nf-core based workflow to analyse nanopore data for fanconi diagnosis.
Basic structure is from following repos (2024-12-10):
1. https://github.com/nf-core/nanoseq
2. https://github.com/dhslab/nf-core-wgsnano

to test the pipeline /root/.local/bin/nextflow run jiangyanyu/nf-core_nano-fanconi -profile test,docker --outdir ./

submit the job to julia

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
