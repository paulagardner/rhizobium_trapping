#!/bin/bash

<<<<<<< Updated upstream
cd ~/Fudge016/TrappingCh2/rhizobium_trapping/data/bamfiles #change shells to where the sorted bam files are

if [ ! -f $~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt ]; then
        echo "file does not exist, making file"
        touch ~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt
     	echo "file made"
else
	echo "file exists"
fi

if [ ! -s $~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt ]; then #this says: if the contents of the file in question are empty, then do the for() loop
	for file in ~/Fudge016/TrappingCh2/rhizobium_trapping/data/bamfiles/*sorted; do
		echo "$file" >> ~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt; #this path is where the .txt is going to live. alter as necessary. 
	done #done is needed whenever there's a for() loop. shellcheck.net is SUPER helpful here 
	echo "list populated with" $( wc -l ~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt | cut -f1 -d' ') "lines"
else 
	echo "file is not empty! It may already contain what you need. it has" $( wc -l ~/Fudge016/TrappingCh2/rhizobium_trapping/data/sorted_bam_files.txt | cut -f1 -d' ') "lines, implying that many BAM files."
=======
cd ~/rhizobium_trapping/data/bamfiles || exit  #change shells to where the sorted bam files are

if [ ! -f $~/rhizobium_trapping/data/sorted_bam_files.txt ]; then
        echo "file does not exist, making file"
        touch ~/rhizobium_trapping/data/sorted_bam_files.txt
        echo "file made"
fi

if [ ! -s $~/rhizobium_trapping/data/sorted_bam_files.txt ]; then #this says: if the contents of the file in question are empty, then do the for() loop
        for file in ~/rhizobium_trapping/data/bamfiles/*sorted; do
                echo "$(pwd)/$file" >> ~/rhizobium_trapping/data/sorted_bam_files.txt; #this path is where the .txt is going to live. alter as necessary. 
        done #done is needed whenever there's a for() loop. shellcheck.net is SUPER helpful here 
        echo "list populated with " $( wc -l ~/rhizobium_trapping/data/sorted_bam_files.txt | cut -f1 -d' ') "lines"
else
        echo "file is not empty! It may already contain what you need. it has" $( wc -l ~/rhizobium_trapping/data/sorted_bam_files.txt | cut -f1 -d' ') "lines, implying that many BAM files."
>>>>>>> Stashed changes

fi #fi is needed to complete a if() statement

#note about echo: >> will append lines, > will overwrite the file. This script should not be able to overwrite a file 
