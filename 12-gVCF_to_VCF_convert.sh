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
