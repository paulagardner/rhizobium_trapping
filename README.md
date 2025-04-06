This is a repo that contains scripts and notes that enable generation of downstream files, plots, and statistical analyses. 

To do this, we used:

### Files:

#### a directory containing .BAM files
created in kbase.us. This required aligning sequence data from each treatment to a chosen reference genome (for us, Rhizobium leguminosarum bv. viciae 3841). In our case, each .BAM file represented one of each of the permutations of:

 - host plant (hairy vetch) ecotype
 - warm/cool treatment
 - soil source

#### Rhizobium leguminosarium bv. viciae 3841 alignment
   - available through [this](https://www.ncbi.nlm.nih.gov/nuccore/NC_008380.1) NCBI link (other websites will have different encoding of the same data)  


### Software: 
* R version 4.4.1
* R packages: vegan, ggplot2, ggpubr, dplyr, tidyverse, tidyr, patchwork, wesanderson
* PoPoolation2 version [1201](https://sourceforge.net/projects/popoolation2/files/popoolation2_1201.zip/download)
* samtools version 1.10


***


## File info, in order of how you run them
[sort.sh](scripts/sort.sh): recursively sort all the .bams in one directory and generate sorted bam files (.bam.sorted) from them. The script puts these files in that directory. Sorting is a step necessary for further samtools operations. Otherwise, mpileup may not work. 
- input:    
   - a directory containing .bam files 
- output:   
   - sorted .bams. We call our sorted bams .bam.sorted, one could change the file extension as desired. 

[make_bam_list.sh](scripts/make_bam_list.sh): generate a .txt file that contains all the paths of the sorted .bams. This is necessary for creating an mpileup file.
- input:
   - directory containing sorted .bam files (the way sort.sh is written, they went into the same directory as the .bams)
- output:
   - .txt file

[create_mpileup.sh](scripts/create_mpileup.sh): generate a samtools mpileup file, which collates the reads from all .bam files into one file.
- input:   
   - directory of sorted .bam files
   - .txt file containing baths to those bam files
- output:   
   - .mpileup file


[sync.sh](scripts/sync.sh): create a sync file from the mpileup, so that popoolation2 can use it.
- input:   
   - .mpileup file
- output:   
   - .sync file

[alellefreqs.sh](scripts/allelefreqs.sh): calculate allele frequency differences. 
- input:
   - .sync file 
- output: 
   - rc, pwc, and params file

[SNP_RDA_RF.R](scripts/SNP_RDA_RF.R): take allele frequency differences, and perform RDA and permanova analyses in R. 

- input: 
   - rc file
- output:
   - plots, statistics


[mpileup_to_vcf.sh](scripts/mpileup_to_vcf.sh): take .mpileup file and convert to vcf for ease of filtering

- input: 
   - .mpileup file
- output:
   - .vcf (with many default filters)
   
[missingness_filter.sh](scripts/missingness_filter.sh): remove variants that have more than 30% missingness 
- input: 
   - .vcf file
- output:
   - filtered .vcf






