#!/usr/bin/bash 
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --tmp=10g


cd ~/Fudge016/TrappingCh2/rhizobium_trapping/data/bamfiles
#this is for the HPC we're working with- using an old version of samtools (1.10) gave a 
# 'libcrypto.so.10' error 
module load samtools/1.16.1-gcc-8.2.0-egljrr3

#this should not generate extraneous files, as whenever <something>.sorted is generated, 
#it will overwrite the previous instance of <something>.sorted
#however, this script doesn't fare well with leftover .tmp files from 
#unfinished sort processes. be sure to remove if this is the case.
for f in *.bam; 
	do  samtools sort -o ${f/bam/sorted} $f #https://www.biostars.org/p/9537028/
done
