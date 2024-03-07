#!/bin/bash

#load modules
module load SAMtools/1.8-foss-2016b

#set variables to reference files
ref=<path>/human_g1k_v37_decoy.fasta
pos=<path>/dbsnp_138.b37.pos
snp=<path>/dbsnp_138.b37.snp

#activate conda for some(?) reason
source activate 

#create read pileup with sam tools and call pseudohaploid genotypes with pileupcaller into eigenstrat fileset
samtools mpileup -R -B -q 1 \
	--positions ${pos} \
	--fasta-ref ${ref} \
	--bam-list newBamList.txt \
	| pileupCaller --sampleNameFile newSamplelist.txt \
	--snpFile ${snp} \
	--randomHaploid \
	--samplePopName BENCHMARK \
	--eigenstratOut <name_prefix>.pseudohap
