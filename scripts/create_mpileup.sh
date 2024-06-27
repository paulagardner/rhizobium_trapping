#!/bin/sh
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shibainu372@gmail.com
module load samtools/1.16.1-gcc-8.2.0-egljrr3
#the long file in paula_dir is the reference genome, against which samtools will compare the sorted .bam files. i
samtools mpileup -b ~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt \
-f ~/Fudge016/TrappingCh2/rhizobium_trapping/data/refseq_reference_coded_with_N/ncbi_dataset/data/GCF_000009265.1/GCF_000009265.1_ASM926v1_genomic.fna \
-o ~/Fudge016/TrappingCh2/rhizobium_trapping/data/documentation-test.mpileup #alter output file names as necessary 
