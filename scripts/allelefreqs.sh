#!/bin/bash

#running into some issues: https://sourceforge.net/p/popoolation2/tickets/32/?page=1

perl /gpfs/bio/xrq24scu/tools/popoolation2_1201/snp-frequency-diff.pl \
	--input /gpfs/bio/xrq24scu/rhizobium_trapping/data/2025revisions.sync \
	--output-prefix /gpfs/bio/xrq24scu/rhizobium_trapping/data/2025revisions \
	--max-coverage 200
