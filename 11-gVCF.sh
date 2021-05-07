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
