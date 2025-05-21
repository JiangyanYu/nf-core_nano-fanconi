nf-core based workflow to analyze nanopore long-read sequencing data for fanconi diagnosis.
Basic steps to use the workflow:
1) Download and unzip package from github
2) Adjusting files accordingly
   --
   --
   --
4) Perform the analysis by runing: 

Basic structure is from following repos (2024-12-10):
1. https://github.com/nf-core/nanoseq
2. https://github.com/dhslab/nf-core-wgsnano
2.1 since the firewall is restricted to docker, I have simply copied the dorado image from dhslab (https://github.com/dhslab/dhslab-docker-images/pkgs/container/docker-dorado) to my docker hub, version 241016.
2.2 same for https://github.com/dhslab/dhslab-docker-images/pkgs/container/docker-whatshap, version 240302.

Reminder:
1) Deepvariant model_type is set as WGS
2) the current annotsv docker image does not contain annotationsDir, thus need to be installed first (https://github.com/lgmgeo/AnnotSV/blob/master/bin/INSTALL_annotations.sh). Then manually change the directory in the test.config file.

Needed input files:
1. samplesheet.csv: raw sequencing files
1.1 id is used to merge samples from different flow cells
1.2 sample shall be different from id, otherwise merge_bam_sample will through error
2. genmexbfx.config: reference files, annotation database, dorado setting
3. 

