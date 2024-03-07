# Modern-shotgun-pipeline

## Make Sample Sheet

Sample sheet format. (index should be reverse complement of actual index molecule) library prep check documentation to make sure this is right, sometimes they provide index already converted to reverse compliment, sometimes not. Can run both to verify you used the right one.

```bash
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

must load arch/haswell first to access bcl2fastq and other modules

```bash
module load arch/haswell 
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

```bash
<>_I1_001.fastq.gz  #index used, irrelevant file
<>_R1_001.fastq.gz  #forward read  
<>_R2_001.fastq.gz  #reverse read
```

## fastp

This script is set up for just one of my samples, because at the time I did them one by one, which I wouldn't recommend.

```bash
cd /hpcfs/users/a1717363/IncaModern/Lane_6/PUN76/

fastp --verbose \
--thread 16 \
-g -x -c -h PUN76.fastp.html -j PUN76_fastp.json \
--in1 ./PUN76_S1_L006_R1_001.fastq.gz --in2 ./PUN76_S1_L006_R2_001.fastq.gz \ 
--out1 ./PUN76_R1.fastp.fastq.gz \
--out2 ./PUN76_R2.fastp.fastq.gz
```

 -g trim poly G commonly observed in Nextseq/NovaSeq sequencing \
 -x trim_poly X \
 -c correction in overlapped region default 30 \

Expected output files:

```bash
<>.fastp.html    #individual report file
<>.json     #report file input for multiqc which gives interactive report file 
<>_R1.fastp.fastq.gz    #forward read  
<>_R2.fastp.fastq.gz  #reverse read
```

## multiqc

Requires multiqc installed with conda \
Outputs cool html file that compares all your samples in the run directory \
Can set up to compare different seq lanes if necessary \
Put all .json files in one directory and run:

```bash
multiqc .
```

## Mapping (bwa-mem)

```bash
module load BWA/0.7.17-foss-2016b
module load SAMtools/1.9-foss-2016b

ref=/<path>/GRCh38_full_analysis_set_plus_decoy_hla.fa
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

```bash
<ref>.dict  #dictionary
<ref>.fa  #actual fasta file, run 'bwa index ref.fa' for the rest of the files
<ref>.fa.amb  #text file, to record appearance of N (or other non-ATGC) in the ref fasta.
<ref>.fa.ann  #text file, to record ref sequences, name, length, etc.
<ref>.fa.bwt  #the Burrows-Wheeler transformed sequence.
<ref>.fa.fai  #the fasta index, tab delimited file where columns are: 
  #1 - The contig name, 
  #2 - number of bases in the contig, 
  #3 - byte index of the file where the contig sequence begins. 
  #4 - bases per line in the FASTA file, 
  #5 - bytes per line in the FASTA file
<ref>.fa.pac  #binary, packaged sequence (four base pairs encode one byte).
<ref>.fa.sa  #binary, suffix array index
```

## Sort bam files

sort mapped reads by genome coordinates

```bash
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

```bash
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

Identifies and marks duplicated reads in bam file. \
Could simply remove duplicates instead (see below).

```bash
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

```bash
<>.bam
<>.bam.bai
```

## Remove duplicates

Remove duplicated reads in bam file. \
The first time I processed data I skipped this step and had duplicate variants to weed out later. \

```bash
#conda-dependency: biobambam 2.0.87

cd /hpcfs/users/a1717363/IncaModern

indir=./07-indelRealign
outdir=./08-removedup
samples=(Ka24 K24 K29 PUN68 PUN76 PUN67)
parallel=6

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    bammarkduplicates2 index=1 \
 rmdup=1 \
 markthread=8 \
 I=$indir/${individ}_indelReal.bam \
 O=$outdir/${individ}_indelReal.removedup.bam \
 D=$outdir/${individ}.removed &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done
```

## Realign around INDELs

Realign bam file around indels. The script first writes a file of interval regions and then realigns. \
I also accidentally skipped this step on my first run through of the pipeline.

```bash
module load parallel/20191022
module load picard/2.23.3
idat=/hpcfs/users/a1717363/IncaModern/06-bqsr
odir=/hpcfs/users/a1717363/IncaModern/07-indelRealign
ref=/hpcfs/users/a1717363/mapping_resources/GRCh37/human_g1k_v37_decoy.fasta
gatk3_jar=$EBROOTGATK/GenomeAnalysisTK.jar

sample=(PUN67 Ka24 K24 K29 PUN68 PUN76)
parallel=6

for (( i=0 ; i<${#sample[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
 java -Xmx8G -jar ${gatk3_jar} \
   -T RealignerTargetCreator \
   -R ${ref} \
   --num_threads 8 \
   --mismatchFraction 0.30 \
   --maxIntervalSize 650 \
   --allow_potentially_misencoded_quality_scores \
   -I ${idat}/${sample}_bqsr.bam \
   -o ${odir}/${sample}_indelReal.intervals
        pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done

for (( i=0 ; i<${#sample[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
 java -Xmx8G -jar ${gatk3_jar} \
   -T IndelRealigner \
   -R ${ref} \
   -model USE_READS \
   -compress 0 \
   --filter_bases_not_stored \
   --allow_potentially_misencoded_quality_scores \
   -I ${idat}/${sample}_bqsr.bam \
   -targetIntervals ${odir}/${sample}_indelReal.intervals \
   -o ${odir}/${sample}_indelReal.bam

 cp -v \
   ${odir}/${sample}_indelReal.bai \
  ${odir}/${sample}_indelRealign.bam.bai &
  pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done
```

## BQSR (Base Quality Score Recalibration)

Uses GATK (Genome Analysis ToolKit) **must use GATK/3.5-Java-1.8.0_71 \
Detects systematic errors made by the sequencing machine in assigning quality scores of each base call and recalibrates the base call.

```bash
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

```bash
.bam
.bai
.bam.bai
```

## QC Mapping

Requires qualimap installed with conda \
Outputs a .html file per sample to download and view, providing quality control metrics on the mapping

```bash
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

```bash
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

```bash
<>.g.vcf.gz
<>.g.vcf.gz.tbi
```

## convert gVCF to VCF

```bash
module load Java/1.8.0_121

cd /hpcfs/users/a1717363/IncaModern/

ref=/<path>/GRCh38_full_analysis_set_plus_decoy_hla.fa
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

```bash
#CHR #POS
1 75266
X 83890
Y 73952
```

```bash
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
 --positions ${sites}/1240K_sites.txt \ #coordinates of SNP sites
 --recode \    #essential to write output file, otherwise no output file is written
 --out ${filterdir}/${sample}_filtered_1240K
 pids[$j]=$!
((i=i+1))
done
```

## Merge VCFs

Check actual vcf file format by `htsfile <file>`
Should output: `VCF version 4.2 BGZF-compressed variant calling data` and be accompanied by a *.vcf.gz.csi index file. \
(htsfile on this outputs: `CSI version 1 compressed index data`) \
If this is not the case run:

```bash
mv file.vcf.gz plain.vcf
bcftools view -Oz -o compressed.vcf.gz plain.vcf
htsfile compressed.vcf.gz
bcftools index compressed.vcf.gz
```

To merge files:

```bash
bcftools merge -m all -l input_list.txt -Oz -o <output_name>.vcf.gz
```

Where `input_list.txt` is a text file with one vcf per line.

## Pseudohaploid Variant Call

Need to input: \
A reference `fasta`, same file as used for mapping.\
`.pos` and `.snp` file of the SNPs you want to call. \
The `.snp` file has the format:

```bash
rs376007522   1   0   10109   A   T
rs368469931   1   0   10139   A   T
rs371194064   1   0   10150   C   T
```

The `.pos` file should have one variant pe line with the Chromosome listed in the first column and position in the second. \
You can write this file from the `.snp` file with the awk command: `awk '{print $2 "\t" $4}' *.snp > *.pos`

`newBamList.txt` is a list of all `.bam` files you wish to call pseudohaploid variants from. One file per line with the absolute path specified. \
`newSamplelist.txt` is a list of sample names, in the same order as the bam file names. Usually the same wile but with the path and file extensions removed.

Now run the script:

```bash
module load SAMtools/1.8-foss-2016b

ref=<path>/human_g1k_v37_decoy.fasta
pos=<path>/dbsnp_138.b37.pos
snp=<path>/dbsnp_138.b37.snp

source activate 
samtools mpileup -R -B -q 1 \
 --positions ${pos} \
 --fasta-ref ${ref} \
 --bam-list newBamList.txt \
 | pileupCaller --sampleNameFile newSamplelist.txt \
 --snpFile ${snp} \
 --randomHaploid \
 --samplePopName BENCHMARK \
 --eigenstratOut <name_prefix>.pseudohap
```

Expect EIGENSTRAT output filest (`.snp`, `.geno`, `.ind` files).
