#!/bun/bash

#load modules
module load arch/haswell 
module load bcl2fastq2/2.19.1

#use BCL2fastq to convert bcl files from sequencer to fastq
bcl2fastq \
--runfolder-dir <path_to_directory_containing_bcl files> \
--output-dir <path_to_output_directory_for_fastq_files \
--sample-sheet <path_to>/SampleSheet.csv \
--ignore-missing-positions \
--ignore-missing-controls \
--ignore-missing-filter \
--ignore-missing-bcls \
--barcode-mismatches=0,1,2