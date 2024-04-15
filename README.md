# START HERE 

This is a repo that contains scripts and notes that enable, wherever feasible, easily-to-modify generation of downstream files, plots, and statistical analyses. 

To do this, we used:

### Files:
- #### a directory containing our .BAM files.
  -ours were created by tools at kbase.us. This required aligning sequence data from each treatment to a chosen reference genome (for us, Rhizobium leguminosarum ==[SPECIFIC REFERENCE INFO HERE]==). In our case, each .BAM file represented one of all the combinations of:
    - vetch ecotype
    - warm/cool treatment
    - soil origin
    - [forgetting one more detail
- #### a reference genome

### Software: 
- Rstudio version ==FILL VERSION HERE==
- PoPoolation2 version ==FILL VERSION HERE==
- ==ETC I'M FORGETTING==
At time of publishing, the following should create a conda environment with the software versions one might use to perform these analyses:
```
```
________________






## File info
[meta.sh](https://github.com/paulagardner/rhizobium_trapping/blob/main/meta.sh) contains a sequential list of all bash commands that might be run to generate or re-generate the data. Currently, it is not configured for use in slurm. 

[To-do list](https://github.com/paulagardner/rhizobium_trapping/blob/main/Todo.md)

[Methods notes](https://github.com/paulagardner/rhizobium_trapping/blob/main/methods.md)
