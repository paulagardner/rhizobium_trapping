#!/bin/sh
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shibainu372@gmail.com

#popoolation2_path = ~/Fudge016/Software/PoPoolation/popoolation2_1201/
java -ea -Xmx7g -jar ~/Fudge016/Software/PoPoolation/popoolation2_1201/mpileup2sync.jar --input ~/Fudge016/TrappingCh2/rhizobium_trapping/data/all_bams.mpileup --output ~/Fudge016/TrappingCh2/rhizobium_trapping/data/sync_all_bams.sync --fastq-type sanger --min-qual 20 --threads 8
