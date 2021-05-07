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
