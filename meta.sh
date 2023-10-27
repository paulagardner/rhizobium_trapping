
#unzip your reference files: 
#unzip GCF_000009265.1.zip

cd Fudge016/TrappingCh2/rhizobium_trapping/shell_scripts
#sort t	e .bam files (samtools mpileup might throw an error if this step is skipped:
sbatch sort.sh

#create a list of these sorted .txt files
bash make_bam_list.sh

#create an mpileup file from a list of .bam files of the populations:
sbatch create_mpileup.sh

#create a sync file 
bash sync.sh

cd Fudge016/TrappingCh2/rhizobium_trapping
#plot things in R:
module load R
#instructions for installing any needed R packages: msi.umn.edu/sw/r

#take sync file and calculate paired FST values from it (currently stored in a .rds file in /data) 
Rscript create_paired_fsts.R


#FILES:
#AUMerit_Rosemount_Cold.FASTQ.zip -- kbase file rebecca uploaded, a sample of the R1 and R2 reads combined into one fasta
#e49a2507-e60c-4b46-bfb7-8e81cbf0076a.inter.fastq -- the unzipped version of the above file
#rleg3841.zip : ncbi reference, with a different encoding than the AUMerit files
