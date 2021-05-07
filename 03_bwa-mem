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
