# Modern-shotgun-pipeline

## Make Sample Sheet
Sample sheet format. (index should be reverse complement of actual index molecule) library prep check documentation to make sure this is right, sometimes they provide index already converted to reverse compliment, sometimes not. Can run both to verify you used the right one. 
```
[Header],,,,
IEMFileVersion,4,,,
Investigator Name,Robbi,,,
Experiment Name,IncaModern_shotgun,,,
Date,26/01/2021,,,
Workflow,GenerateFASTQ,,,
,,,,
[Reads],,,,
151,,,,
151,,,,

[Settings]
CreateFastqForIndexReads,1,,,
,,,,
[Data],,,,
Sample_ID,DNA_Extract_Name,I7_Index_ID,index
PUN76,Ap3,UDP0289,GAACAATTCC
PUN67,Ap10,UDP0290,TGTGGTCCGG
PUN68,Ap18,UDP0291,CTTCTAAGTC
K24,INCP-2,UDP0292,AATATTGCCA
K24a,INCP-6,UDP0293,TCGTGCATTC
K29,INCP-10,UDP0294,AAGATACACG
```

## Demultiplexing with bcl2fastq

```
module load arch/haswell #must load arch/haswell first to access bcl2fastq and other modules
module load bcl2fastq2/2.19.1

bcl2fastq \
--runfolder-dir <path_to_directory_containing_bcl files> \
--output-dir <path_to_output_directory_for_fastq_files \
--sample-sheet <path_to>/SampleSheet.csv \
--ignore-missing-positions \
--ignore-missing-controls \
--ignore-missing-filter \
--ignore-missing-bcls \
--barcode-mismatches=0,1,2
```
Expected output files: 
```
<>_I1_001.fastq.gz 	#index used, irrelevant file
<>_R1_001.fastq.gz 	#forward read  
<>_R2_001.fastq.gz 	#reverse read
```

## fastp
This script is set up for just one of my samples, because at the time I did them one by one, which I wouldn't recommend.
```
cd /hpcfs/users/a1717363/IncaModern/Lane_6/PUN76/

fastp --verbose \
--thread 16 \
-g -x -c -h PUN76.fastp.html -j PUN76_fastp.json \ 	#g, x, c, h used because novaseq data
--in1 ./PUN76_S1_L006_R1_001.fastq.gz --in2 ./PUN76_S1_L006_R2_001.fastq.gz \ 
--out1 ./PUN76_R1.fastp.fastq.gz \
--out2 ./PUN76_R2.fastp.fastq.gz
```
	-g trim poly G commonly observed in Nextseq/NovaSeq sequencing \
	-x trim_poly X \
	-c correction in overlapped region default 30 \

Expected output files: 
```
<>.fastp.html  		#individual report file
<>.json   		#report file input for multiqc which gives interactive report file 
<>_R1.fastp.fastq.gz   	#forward read  
<>_R2.fastp.fastq.gz 	#reverse read
```
## multiqc
Requires multiqc installed with conda \
Outputs cool html file that compares all your samples in the run directory \
Can set up to compare different seq lanes if necessary \
Put all .json files in one directory and run:
```
multiqc .
```

## Mapping (bwa-mem)
```
module load BWA/0.7.17-foss-2016b
module load SAMtools/1.9-foss-2016b

ref=/hpcfs/groups/acad_users/Refs/Homo_sapiens/GATK/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
fastpdir1=/hpcfs/users/a1717363/IncaModern/01-fastp/Lane_6/
mapped1=/hpcfs/users/a1717363/IncaModern/02-mapped/Lane6/
samples=(K24 K24a K29 PUN67 PUN68 PUN76)

i=1
for sample in ${samples[@]}
do
    echo =========================================
    echo Sample $i/${#samples[@]} $sample
    echo =========================================
    bwa mem $ref \
    $fastpdir1/*${sample}*R1*gz \
    $fastpdir1/*${sample}*R2*gz \
    -t 54 -T 0 \
    -R "@RG\tID:${sample}_L6\tSM:$sample\tLIB:NexFlexRD\tPL:Illumina" | samtools view -OBAM > $mapped1/${sample}.bam
    ((i=i+1))
done
```
Index files required to run bwa:
```
<ref>.dict 	#dictionary
<ref>.fa 	#actual fasta file, run 'bwa index ref.fa' for the rest of the files
<ref>.fa.amb 	#text file, to record appearance of N (or other non-ATGC) in the ref fasta.
<ref>.fa.ann 	#text file, to record ref sequences, name, length, etc.
<ref>.fa.bwt 	#the Burrows-Wheeler transformed sequence.
<ref>.fa.fai 	#the fasta index, tab delimited file where columns are: 
		#1 - The contig name, 
		#2 - number of bases in the contig, 
		#3 - byte index of the file where the contig sequence begins. 
		#4 - bases per line in the FASTA file, 
		#5 - bytes per line in the FASTA file
<ref>.fa.pac 	#binary, packaged sequence (four base pairs encode one byte).
<ref>.fa.sa 	#binary, suffix array index
```

## Sort bam files
sort mapped reads by genome coordinates
```
cd /hpcfs/users/a1717363/IncaModern

mapped1=./02-mapped/Lane6
mapped2=./02-mapped/Lane7
sortdir=./03-sort
samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6

#must install biobambam on conda for bamsort to work

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    bamsort threads=24 <$mapped1/${individ}.bam> $sortdir/Lane6/${individ}sorted.bam &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    bamsort threads=24 <$mapped2/${individ}.bam> $sortdir/Lane7/${individ}sorted.bam &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done
```

## Merge .bam files
Merge data from different lanes. \
Will ouput one .bam per sample.
```
module load SAMtools/1.9-foss-2016b

cd /hpcfs/users/a1717363/IncaModern/

sort1=./03-sort/Lane6
sort2=./03-sort/Lane7
merge=./04-mergebam

samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    samtools merge -c -@ 5 $merge/${individ}_merged.bam $sort1/${individ}sorted.bam $sort2/${individ}sorted.bam &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done
```

## Mark Duplicates
Identifies and marks duplicated reads in bam file
```
cd /hpcfs/users/a1717363/IncaModern/

merge=./04-mergebam
markdup=./05-markdup
samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    # run up to $parallel fastp processes in parallel
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    bammarkduplicates2 index=1 markthread=8 I=$merge/${individ}_merged.bam O=$markdup/${individ}_markdup.bam &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done
```
Expected output files:
```
<>.bam
<>.bam.bai
```

## BQSR (Base Quality Score Recalibration)
Uses GATK (Genome Analysis ToolKit) **must use GATK/3.5-Java-1.8.0_71 \
Detects systematic errors made by the sequencing machine in assigning quality scores of each base call and recalibrates the base call.
```
ml arch/haswell
ml arch/arch/haswell
ml modulefiles/arch/haswell
ml GATK/3.5-Java-1.8.0_71

cd /hpcfs/users/a1717363/IncaModern/

ref=/hpcfs/groups/acad_users/Refs/Homo_sapiens/GATK/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
knownsites=/hpcfs/users/a1717363/mapping_resources/00-All.vcf.gz
indata=./05-markdup
outdata=./06-bqsr
recal=./06-bqsr/recal
tmp=./06-bsqr/tmp
samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    java -Xmx16G -Djava.io.tmpdir=$tmp/ -jar $EBROOTGATK/GenomeAnalysisTK.jar -T BaseRecalibrator \
    -nct 8 -R $ref -I $indata/${individ}_markdup.bam --knownSites $knownsites -o $recal/${individ}.table &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done

samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6
for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    java -Xmx16G -Djava.io.tmpdir=$tmp/ -jar $EBROOTGATK/GenomeAnalysisTK.jar -T PrintReads \
    -nct 16 -R $ref -I $indata/${individ}_markdup.bam -BQSR $recal/${individ}.table -o $outdata/${individ}_bqsr.bam &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done

module load SAMtools/1.9-foss-2016b
samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    samtools index $outdata/${individ}_bqsr.bam &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done
```
Expected output files
```
.bam
.bai
.bam.bai
```

## QC Mapping
Requires qualimap installed with conda \
Outputs a .html file per sample to download and view, providing quality control metrics on the mapping
```
ml arch/haswell
ml arch/arch/haswell
ml modulefiles/arch/haswell
ml GATK/3.5-Java-1.8.0_71

cd /hpcfs/users/a1717363/IncaModern/

indata=./06-bqsr
outdata=.06-bqsr/qualimap

samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    java -jar $EBROOTGATK/GenomeAnalysisTK.jar qualimap bamqc -bam $indata/${individ}_bqsr.bam --java-mem-size=16G -outdir $outdata/${individ} -outformat HTML &
    pids[$j]=$!
    done
for pid in ${pids[@]}; do wait $pid; done

done
```

## gVCF
genomiic Variant Call Format \
similar to VCF file but stores information for variant and non-variant sites
```
module load Java/1.8.0_121

cd /hpcfs/users/a1717363/IncaModern/

ref=/hpcfs/groups/acad_users/Refs/Homo_sapiens/GATK/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
knownsites=/hpcfs/users/a1717363/mapping_resources/00-All.vcf.gz
indata=./06-bqsr
outdata=./07-gvcf
tmp=./07-gvcf/tmp

samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    java -Xmx32g -Djava.io.tmpdir=$tmp/ -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller \
    --genotyping_mode DISCOVERY -ERC GVCF --pcr_indel_model NONE \
    -R $ref -I $indata/${individ}_bqsr.bam -o $outdata/${individ}.g.vcf.gz \
    -rf BadCigar \
    -A DepthPerAlleleBySample \
    -A DepthPerSampleHC \
    -A InbreedingCoeff \
    -A StrandBiasBySample \
    -A Coverage \
    -A FisherStrand \
    -A MappingQualityRankSumTest \
    -A MappingQualityZero \
    -A QualByDepth \
    -A RMSMappingQuality \
    -A ReadPosRankSumTest &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done
```
Expected output files:
```
<>.g.vcf.gz
<>.g.vcf.gz.tbi
```
## convert gVCF to VCF
```
module load Java/1.8.0_121

cd /hpcfs/users/a1717363/IncaModern/

ref=/hpcfs/groups/acad_users/Refs/Homo_sapiens/GATK/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
indata=./07-gvcf
outdata=./08-gVCF-VCF
tmp=./08-gVCF-VCF/tmp

samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    java -Xmx32g -Djava.io.tmpdir=$tmp -jar $EBROOTGATK/GenomeAnalysisTK.jar -T GenotypeGVCFs \
    -nt 16 -R $ref -V $indata/${individ}.g.vcf.gz \
    --includeNonVariantSites \		#includes all sites whether there is a variant in the sample or not. Bigger files, more time. Can be useful if filtering to a certain snp set.
    -o $outdata/${individ}.vcf.gz &
    pids[$j]=$!
    done
for pid in ${pids[@]}; do wait $pid; done

done
```
## Filter VCF
Filter VCF by sites in the 1240K dataset
Requires input file of positions to filter by, tab delimited with chromosome in the first column and position in bp in the second. No header.
eg. format:
```
#CHR	#POS
1	75266
X	83890
Y	73952
```
```
module load vcftools/0.1.12a-GCC-5.3.0-binutils-2.25 

cd /hpcfs/users/a1717363/IncaModern

vcfdir=./08-VCF
filterdir=./09-filter_VCF_1240K
sites=./1240K-sites
samples=(Ka24 K24 K29 PUN67 PUN76 PUN68)
parallel=6

i=1
for sample in ${samples[@]}
do
vcftools --gzvcf ${vcfdir}/${sample}.vcf.gz\
	--positions ${sites}/1240K_sites.txt \	#coordinates of SNP sites
	--recode \				#essential to write output file, otherwise no output file is written
	--out ${filterdir}/${sample}_filtered_1240K
	pids[$j]=$!
((i=i+1))
done
```
## Merge VCFs
Check actual vcf file format by `htsfile <file>` \
Should output: `VCF version 4.2 BGZF-compressed variant calling data` and be accompanied by a *.vcf.gz.csi index file. \
(htsfile on this outputs: `CSI version 1 compressed index data`) \
If this is not the case run:
```
mv file.vcf.gz plain.vcf
bcftools view -Oz -o compressed.vcf.gz plain.vcf
htsfile compressed.vcf.gz
bcftools index compressed.vcf.gz
```
To merge files:
```
bcftools merge -m all -l input_list.txt -Oz -o IncaModern.vcf.gz
```
Where `input_list.txt` is a text file with one vcf per line.

