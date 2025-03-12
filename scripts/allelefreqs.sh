#!/bin/bash

#running into some issues: https://sourceforge.net/p/popoolation2/tickets/32/?page=1

perl ~/software/popoolation2_1201/snp-frequency-diff.pl \
	--input ~/rhizobium_trapping/data/sync_all_bams.sync \
	--output-prefix ~/rhizobium_trapping/data/samtools_1_10_test \
	--max-coverage 200
