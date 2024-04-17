#!/bin/bash 

cd ~/Fudge016/TrappingCh2/bamfiles
module load samtools/1.10

#this should not generate extraneous files, as whenever <something>.sorted is generated, 
#it will overwrite the previous instance of <something>.sorted
for f in *.bam; 
	do  samtools sort -o ${f/bam/sorted} $f #https://www.biostars.org/p/9537028/
done
