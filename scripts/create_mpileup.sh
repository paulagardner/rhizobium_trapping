#!/bin/sh 

#the long file in paula_dir is the reference genome, against which samtools will compare the sorted .bam files. i
samtools mpileup -b /mnt/c/Users/Paula/Desktop/TrappingCh2/rhizobium_trapping/data/2025sorted_bamfiles.txt \
-f  /mnt/c/Users/Paula/Desktop/TrappingCh2/rhizobium_trapping/data/refseq_reference_coded_with_N/ncbi_dataset/data/GCF_000009265.1/GCF_000009265.1_ASM926v1_genomic.fna \
-o /home/paula/rhizobium/rhizobium_trapping/data/2025revisions.mpileup #alter output file names as necessary 