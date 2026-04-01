#!/bin/bash

for file in $(ls *.fa | sed -e 's/.contigs.fa//' | sort -u)
do
	gunzip -k ~/work_space01/cleandata/gzipfile/${file}_clean_R1.fq.gz -c > ${file}_clean_1.fastq
	gunzip -k ~/work_space01/cleandata/gzipfile/${file}_clean_R2.fq.gz -c > ${file}_clean_2.fastq
	time metawrap binning -o ./initial_bining/${file} -t 64 -a ${file}.contigs.fa --metabat2 --maxbin2 --concoct ${file}_clean_1.fastq ${file}_clean_2.fastq
	sleep 10
	rm ${file}_clean_*.fastq
	echo "############# FASTQ File Removed #################"
done
