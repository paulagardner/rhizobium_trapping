#!/bin/sh

module load samtools/1.10

#the long file in paula_dir is the reference genome, against which samtools will compare the sorted .bam files. i
samtools mpileup -b ~/rhizobium_trapping/data/sorted_bam_files.txt \
-f ~/rhizobium_trapping/data/refseq_reference_coded_with_N/ncbi_dataset/data/GCF_000009265.1/GCF_000009265.1_ASM926v1_genomic.fna \
-o ~/rhizobium_trapping/data/documentation-test.mpileup #alter output file names as necessary 
