# Copilot Instructions for nf-core_nano-fanconi

## Project Overview
- Nextflow DSL2 pipeline for nanopore long-read sequencing, focused on Fanconi anemia diagnosis.
- Modular workflow: local modules in `modules/local/`, nf-core modules in `modules/nf-core/`.
- Main workflow logic in `workflows/nano-fanconi.nf`.

## Architecture & Data Flow
- Input: `assets/samplesheet.csv` (absolute paths), reference genome in `profile.config`.
- Resource and tool selection: `nextflow.config` (true/false for modules, resource limits).
- Steps: input validation, basecalling, alignment, variant calling, per-chromosome phasing, merging, haplotagging, annotation, QC, reporting.
- Output: per-sample folders, MultiQC, pipeline info, intermediate files.

## Developer Workflows
- Run: `nextflow run ./nf-core_nano-fanconi/ -profile <profile>,docker --outdir <output>`
- Debug: Use Nextflow's reporting (execution_report, timeline, trace, DAG).
- Docker required for reproducibility.

## Project-Specific Patterns
- Channels set/passed via `.set {}` and `.out`. Ensure channel names match between steps.
- Conditional logic for input formats and module selection is common.
- Resource params (`max_cpus`, `max_memory`, `max_time`) in `nextflow.config`.
- Output grouped by sample and step; MultiQC and PycoQC for summary.

## Integration Points & Dependencies
- Tools: Nextflow, Docker, DeepVariant, AnnotSV, Dorado, PycoQC, MultiQC, Sawfish, Whatshap, Samtools, Bcftools, Tabix.
- Annotation DBs may need manual install (see README).
- Validate module inputs/outputs for file types/structure.

## Examples & Conventions
- Add analysis: create module in `modules/local/`, include in `nano-fanconi.nf`, wire channels.
- Change resources: update `nextflow.config`, `conf/base.config`.
- New input formats: add logic in `nano-fanconi.nf`, update `docs/usage.md`.

## Key Files
- `workflows/nano-fanconi.nf`: Workflow logic/channel wiring.
- `modules/local/`: Custom modules.
- `assets/samplesheet.csv`: Input samples.
- `profile.config`, `nextflow.config`, `conf/base.config`: Config/resource settings.
- `docs/usage.md`, `docs/output.md`: Documentation.
- `README.md`: Project overview/setup.

---
For questions about module integration, resource settings, or output structure, refer to workflow/config files. If unclear, ask for examples from maintainers.
