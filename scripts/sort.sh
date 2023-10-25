#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fudge016@umn.edu
cd ~/Fudge016/TrappingCh2/bamfiles
module load samtools/1.10

#note that unlike the command we give to make the LIST of bamfiles, 
#this will not generate extraneous files, as whenever <something>.sorted is generated, 
#it will overwrite the previous instance of <something>.sorted
for f in *.bam; 
	do  samtools sort -o ${f/bam/sorted} $f #https://www.biostars.org/p/9537028/
done
