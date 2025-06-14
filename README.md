nano-fanconi is an nf-core based workflow to analyze nanopore long-read sequencing data for routine fanconi diagnosis.
![Description](https://github.com/JiangyanYu/nf-core_nano-fanconi/blob/main/docs/workflow_complete_graph.png)
**Basic steps to use the workflow**:
1. Install nextflow according to its manual (https://www.nextflow.io/docs/latest/install.html)
2. Download and unzip (or git clone) nano-fanconi package from github (https://github.com/JiangyanYu/nf-core_nano-fanconi.git)
```
git clone -b main https://github.com/JiangyanYu/nf-core_nano-fanconi
```
4. Adjusting file path accordingly:
   1) Sequencing data path in **/nf-core_nano-fanconi/assets/samplesheet.csv**. The file directory is the absolute path in your file system. id is used to merge samples from different flow cells. Therefore, sample shall be different from id, otherwise merge_bam_sample will throw errors.
   2) Reference genome path in **/nf-core_nano-fanconi/profile.config**. Besides the reference path, annotsvAnnotations database directory as well as the dorado model details can be specified in this file.
   3) Select packages to be used in **/nf-core_nano-fanconi/nextflow.config**. Annotate (true or false) the programs one likes to use. In the default setting, annotation part is marked as false. If needed, installation of annotation database is needed (see below).
   4) Resource specification in **/nf-core_nano-fanconi/conf/base.config**.
5. Run the analysis by following command:
```
   nextflow run ./nf-core_nano-fanconi/ \ # The path to the nano-fanconi package
      -profile nanofanconi,docker \ # Corresponding to the setting in /nf-core_nano-fanconi/profile.config. Beaware that no space between nanofanconi,docker 
      --outdir ./output # Specify the directory for output results.
```

**Notes**
1) the current annotsv docker image does not contain annotationsDir, thus need to be installed first (https://github.com/lgmgeo/AnnotSV/blob/master/bin/INSTALL_annotations.sh). Then manually change the directory in the profile.config file.
2) Deepvariant model_type is set as WGS

**Cite us (to-be-updated):**
Nano-Fanconi: A Nextflow framework for automated analysis of Nanopore based long-read sequencing data for Fanconi anemia diagnosis. 

   
<sub>Basic structure is from following repos (2024-12-10):</sub> \
<sub>1. https://github.com/nf-core/nanoseq</sub> \
<sub>2. https://github.com/dhslab/nf-core-wgsnano</sub> \
<sub>2.1 since the firewall is restricted to docker, I have simply copied the dorado image from dhslab (https://github.com/dhslab/dhslab-docker-images/pkgs/container/docker-dorado) to my docker hub, version 241016.</sub> \
<sub>2.2 same for https://github.com/dhslab/dhslab-docker-images/pkgs/container/docker-whatshap, version 240302.</sub>


