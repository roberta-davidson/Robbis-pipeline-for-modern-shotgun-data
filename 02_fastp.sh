#!/bin/bash

#fastp script for one sample 

fastp --verbose \
--thread 16 \
-g -x -c -h PUN76.fastp.html -j PUN76_fastp.json \
--in1 ./PUN76_S1_L006_R1_001.fastq.gz \
--in2 ./PUN76_S1_L006_R2_001.fastq.gz \
--out1 ./PUN76_R1.fastp.fastq.gz \
--out2 ./PUN76_R2.fastp.fastq.gz

#flags explained by chatgpt:
#    -g or --cut_right: This option specifies the number of bases to trim from the 3' (right) end of each read. It is used to remove low-quality bases from the end of reads.
#    -x or --cut_tail: This option specifies the minimum length threshold for reads after trimming. Any reads shorter than this threshold after trimming will be discarded.
#    -c or --correction: This option enables base correction. Fastp will perform error correction on the input reads, attempting to fix base errors.
#    -h or --html: This option specifies the output HTML report file name. Fastp generates an HTML report summarizing the trimming and filtering statistics, which can be useful for visualizing the quality control results.