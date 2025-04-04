#!/bin/sh

module load OpenJDK/jdk-20.0.2

#two ways of approaching it. The faster way:
#popoolation2_path = ~/popoolation2_1201
#java -ea -Xmx7g -jar /gpfs/software/ada/build/popoolation2-code/mpileup2sync.jar \
#	--input /gpfs/bio/xrq24scu/rhizobium_trapping/data/documentation-test.mpileup \
#	--output /gpfs/bio/xrq24scu/rhizobium_trapping/data/documentation-test_revisions.sync \
#	--fastq-type sanger --min-qual 20 --threads 8

perl /gpfs/bio/xrq24scu/tools/popoolation2_1201/mpileup2sync.pl \
	--fastq-type sanger \
	--min-qual 20 \
	--input /gpfs/bio/xrq24scu/rhizobium_trapping/data/2025revisions.mpileup \
	--output /gpfs/bio/xrq24scu/rhizobium_trapping/data/2025revisions.sync
