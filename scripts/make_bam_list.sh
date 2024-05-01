#!/bin/bash
#SBATCH --time = 2:00:00
#SBATCH --ntasks = 1
#SBATCH --mem = 10g
#SBATCH --tmp = 10g
#SBATCH --mail-type = ALL
#SBATCH --mail-user = shibainu372@gmail.com

#cd ~/Fudge016/TrappingCh2/bamfiles 
if ! test -f $~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt; then
	echo "file does not exist, making file"	
	touch ~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt
	echo "file made"
fi

if [ -s ~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt ]; then #this says: if the contents of the file in question are empty, then do the for() loop
	echo "file is not empty! It may already contain what you need. it has" $( wc -l ~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt | cut -f1 -d' ') "lines";
	
else 
	echo "file is empty, adding entries to list";

	for file in ~/Fudge016/TrappingCh2/bamfiles/*sorted; do
			echo "$file" >> ~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt;
	done #done is needed whenever there's a for() loop. shellcheck.net is SUPER helpful here 
	echo "list populated with " $( wc -l ~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt | cut -f1 -d' ') "lines"
fi #fi is needed to complete a if() statement

#note about echo: >> will append lines, > will overwrite the file. This code is (hopefully) set up to not run if the file is not empty. 
