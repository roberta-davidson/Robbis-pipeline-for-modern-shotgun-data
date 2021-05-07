cd /hpcfs/users/a1717363/IncaModern/Lane_6/PUN76/

fastp --verbose \
--thread 16 \
-g -x -c -h PUN76.fastp.html -j PUN76_fastp.json \
--in1 ./PUN76_S1_L006_R1_001.fastq.gz \
--in2 ./PUN76_S1_L006_R2_001.fastq.gz \ 
--out1 ./PUN76_R1.fastp.fastq.gz \
--out2 ./PUN76_R2.fastp.fastq.gz
