#!/bin/sh

#two ways of approaching it. The faster way:
#popoolation2_path = ~/popoolation2_1201
java -ea -Xmx7g -jar ~/software/popoolation2_1201/mpileup2sync.jar \
	--input /mnt/c/Users/Paula/Desktop/TrappingCh2/rhizobium_trapping/data/documentation-test.mpileup \
	--output /mnt/c/Users/Paula/Desktop/TrappingCh2/rhizobium_trapping/data/old_documentation-test_revisions.sync \
	--fastq-type sanger --min-qual 20 --threads 8
	
#perl ~/software/popoolation2_1201/snp-frequency-diff.pl \
#	--fastq-type sanger \
#	--min-qual 20 \
#	--input /mnt/c/Users/Paula/Desktop/TrappingCh2/rhizobium_trapping/data/documentation-test.mpileup \
#	--output /mnt/c/Users/Paula/Desktop/TrappingCh2/rhizobium_trapping/data/old_documentation-test_revisions.sync


