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
