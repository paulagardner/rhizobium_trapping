#!/bin/sh

#two ways of approaching it. The faster way:
#popoolation2_path = ~/popoolation2_1201
#java -ea -Xmx7g -jar ~/popoolation2_1201/mpileup2sync.jar \
#	--input ~/rhizobium_trapping/data/documentation-test.mpileup \
#	--output ~/rhizobium_trapping/data/sync_all_bams.sync \
#	--fastq-type sanger --min-qual 20 --threads 8

perl ~/popoolation2_1201/mpileup2sync.pl \
	--fastq-type sanger \
	--min-qual 20 \
	--input ~/rhizobium_trapping/data/documentation-test.mpileup \
	--output ~/rhizobium_trapping/data/sync_all_bams.sync


