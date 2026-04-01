#!/bin/bash

for file in $(ls *.fq.gz | sed -e 's/_clean_R1.fq.gz//' -e 's/_clean_R2.fq.gz//' | sort -u)
do 
	time megahit -t 48 -m 0.9 --k-min 29 --min-contig-len 1000 -1 ${file}_clean_R1.fq.gz -2 ${file}_clean_R2.fq.gz -o /home/hxl/work01/my_contigs/${file} --out-prefix ${file}
done
