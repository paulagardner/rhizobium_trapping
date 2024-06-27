# START HERE 

This is a repo that contains scripts and notes that enable, wherever feasible, easily-to-modify generation of downstream files, plots, and statistical analyses. 

To do this, we used:

### Files:
* #### a directory containing our .BAM files.
ours were created by tools at kbase.us. This required aligning sequence data from each treatment to a chosen reference genome (for us, Rhizobium leguminosarum ==[SPECIFIC REFERENCE INFO HERE]==). In our case, each .BAM file represented one of all the combinations of:
    * vetch ecotype
    * warm/cool treatment
    * soil origin
    * [forgetting one more detail
* #### a reference genome

### Software: 
* Rstudio version ==FILL VERSION HERE==
* PoPoolation2 version ==FILL VERSION HERE==
* samtools version 1.16
* ==ETC I'M FORGETTING==
At time of publishing, the following should create a conda environment with the software versions one might use to perform these analyses:
```
```
________________






## File info
[meta.sh](https://github.com/paulagardner/rhizobium_trapping/blob/main/meta.sh) contains a sequential list of all bash commands that might be run to generate or re-generate the data. Currently, it is not configured for use in slurm. 

[sort.sh](scripts/sort.sh): recursively sort all the .bams in one directory and generate sorted bam files (.bam.sorted) from them. The script puts these files in that directory. Sorting is a step necessary for further samtools operations. Otherwise, mpileup may not work. 

[make_bam_list.sh](scripts/make_bam_list.sh): generate a .txt file that contains all the paths of the sorted .bams . This is necessary for creating an mpileup file.

[create_mpileup.sh](scripts/create_mpileup.sh): generate a samtools mpileup file, which collates the reads from all .bam files into one file


Why do we need an mpileup? (I FORGET, EDIT THIS) 





[To-do list](https://github.com/paulagardner/rhizobium_trapping/blob/main/Todo.md)

[Methods notes](https://github.com/paulagardner/rhizobium_trapping/blob/main/methods.md)
