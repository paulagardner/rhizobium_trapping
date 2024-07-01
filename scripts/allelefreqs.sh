#!/bin/bash

#SBATCH --job-name=allele_freq_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=4g

#running into some issues: https://sourceforge.net/p/popoolation2/tickets/32/?page=1

perl ~/Fudge016/Software/PoPoolation/popoolation2_1201/snp-frequency-diff.pl \
	--input ~/Fudge016/TrappingCh2/rhizobium_trapping/data/sync_all_bams.sync \
	--output-prefix ~/Fudge016/TrappingCh2/rhizobium_trapping/data/snp_freq_diffs_test_mincov2 \
	--min-count 1 \
	--min-coverage 2 \
	--max-coverage 200
