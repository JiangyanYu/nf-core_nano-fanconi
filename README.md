nf-core based workflow to analyse nanopore data for fanconi diagnosis.
Basic structure is from following repos (2024-12-10):
1. https://github.com/nf-core/nanoseq
2. https://github.com/dhslab/nf-core-wgsnano

to test the pipeline /root/.local/bin/nextflow run jiangyanyu/nf-core_nano-fanconi -profile test,docker --outdir ./

submit the job to julia, by run_nfcore_nanoseq.sh

Reminder:
1) Deepvariant model_type is set as WGS
