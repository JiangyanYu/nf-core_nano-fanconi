nf-core based workflow to analyse nanopore data for fanconi diagnosis.
Basic structure is from following repos (2024-12-10):
1. https://github.com/nf-core/nanoseq
2. https://github.com/dhslab/nf-core-wgsnano
2.1 since the firewall is restricted to docker, I have simply copied the dorado image from dhslab (https://github.com/dhslab/dhslab-docker-images/pkgs/container/docker-dorado) to my docker hub, version 241016.
2.2 same for https://github.com/dhslab/dhslab-docker-images/pkgs/container/docker-whatshap, version 240302.

Reminder:
1) Deepvariant model_type is set as WGS
2) the current annotsv docker image does not contain annotationsDir, thus need to be installed first (https://github.com/lgmgeo/AnnotSV/blob/master/bin/INSTALL_annotations.sh). Then manually change the directory in the test.config file.

Input files:
1. samplesheet.csv
1.1 id is used to merge samples from different flow cells
1.2 sample shall be different from id, otherwise merge_bam_sample will through error

