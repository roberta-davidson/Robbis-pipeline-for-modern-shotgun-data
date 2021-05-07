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
