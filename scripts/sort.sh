#!/usr/bin/bash 
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --tmp=10g

cd ~/rhizobium_trapping/data/bamfiles
module load samtools/1.10

#this should not generate extraneous files, as whenever <something>.sorted is generated, 
#it will overwrite the previous instance of <something>.sorted
#however, this script doesn't fare well with leftover .tmp files from 
#unfinished sort processes. be sure to remove if this is the case.
for f in *.bam; 
	do  samtools sort -o ${f/bam/sorted} $f #https://www.biostars.org/p/9537028/
done
