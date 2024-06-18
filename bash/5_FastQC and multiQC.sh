#!/bin/bash
#####Quality control - FastQC and multiQC

#####FastQC (can't do in screen)

cd /project/trunc
fastqc -t 6 *.fastq --outdir project/FASTQC



#####multiQC

cd /project/
multiqc .




