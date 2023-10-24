#!/bin/sh
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shibainu372@gmail.com
module load samtools/1.10
samtools mpileup -b ~/Fudge016/TrappingCh2/workflow/data/sorted_bam_files.txt -f ~/Fudge016/TrappingCh2/paula_dir/refseq_reference_coded_with_N/ncbi_dataset/data/GCF_000009265.1/GCF_000009265.1_ASM926v1_genomic.fna -o ~/Fudge016/TrappingCh2/workflow/data/all_bams.mpileup
