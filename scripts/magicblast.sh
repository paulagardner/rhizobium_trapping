#### Set requirements ####
salloc -N 1 -n 20 --mem-per-cpu=2048 -t 12:00:00

#### Activate environment ####
# conda create -n magicblast
# conda activate magicblast 
# conda install -c bioconda magicblast

conda activate magicblast 

REF=/storage/home/jus1990/burghardt/Jeremy/AVT_Fa20_Sp21_Nods/working/USDA1106/GCF_000346065.1_ASM34606v1_genomic.fna
total_files=`ls /storage/home/jus1990/burghardt/Jeremy/AVT_Fa20_Sp21_Nods/trimmed_fastqs/*.atria.fastq.gz | wc -l`
arr=( $(ls /storage/home/jus1990/burghardt/Jeremy/AVT_Fa20_Sp21_Nods/trimmed_fastqs/*.atria.fastq.gz) )


# create blast database from REF for mapping to
makeblastdb -in $REF -out USDA1106_reference -parse_seqids -dbtype nucl

#Map and name sort fastq files by name for fixmate
for ((i=0; i<$total_files; i+=2))
{
sample_name=`echo ${arr[$i]} | awk -F "/" -v OFS='/' '{print $9}' | awk -F "." -v OFS='.' '{print $1}' | awk -F "_" -v OFS='_' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5}'`
echo "[mapping running for] $sample_name"
printf "\n"
magicblast -query ${arr[$i]} -query_mate ${arr[$i+1]} -db USDA1106_reference -num_threads 20 -infmt fastq | samtools sort -n | samtools view -bS > $sample_name.nsort.bam 
}

# Errors for these bams 
rm A16_V_Sp21_4_S86.nsort.bam
rm A20_V_Fa20_3_S213.nsort.bam


# fixmate (fills in mate coordinates and insert size fields)
for ((i=0; i<$total_files; i+=2))
{
sample_name=`echo ${arr[$i]} | awk -F "/" -v OFS='/' '{print $9}' | awk -F "." -v OFS='.' '{print $1}' | awk -F "_" -v OFS='_' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5}'`
echo "[mapping running for] $sample_name"
printf "\n"
samtools fixmate -m $sample_name.nsort.bam $sample_name.fixmate.bam 
}

# coordinate sort | markdup (mark duplicate alignments in a coordinate sorted file)

for ((i=0; i<$total_files; i+=2))
{
sample_name=`echo ${arr[$i]} | awk -F "/" -v OFS='/' '{print $9}' | awk -F "." -v OFS='.' '{print $1}' | awk -F "_" -v OFS='_' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5}'`
echo "[mapping running for] $sample_name"
printf "\n"
samtools sort $sample_name.fixmate.bam -o $sample_name.fixsort.bam
} 

for ((i=0; i<$total_files; i+=2))
{
sample_name=`echo ${arr[$i]} | awk -F "/" -v OFS='/' '{print $9}' | awk -F "." -v OFS='.' '{print $1}' | awk -F "_" -v OFS='_' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5}'`
echo "[mapping running for] $sample_name"
printf "\n"
samtools markdup $sample_name.fixsort.bam $sample_name.markdup.bam
}

#parallel index bams 
ls *markdup.bam | parallel samtools index '{}'

# flagstat and check coverage	
for ((i=0; i<$total_files; i+=2))
{
sample_name=`echo ${arr[$i]} | awk -F "/" -v OFS='/' '{print $9}' | awk -F "." -v OFS='.' '{print $1}' | awk -F "_" -v OFS='_' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5}'`
echo "[mapping running for] $sample_name"
printf "\n"
samtools flagstat $sample_name.markdup.bam 
samtools coverage $sample_name.markdup.bam 
}

ls *markdup.bam > bam.list

#create vcf file per chrom (resource intensive)
bcftools mpileup --threads 20 -Ou -f $REF --max-depth 500 -b bam.list --annotate FORMAT/AD -r NC_020528.1 | bcftools call -mv -Oz -o sino_mel_NC_020528.1.vcf.gz
bcftools mpileup --threads 20 -Ou -f $REF --max-depth 500 -b bam.list --annotate FORMAT/AD -r NC_020560.1 | bcftools call -mv -Oz -o sino_mel_NC_020560.1.vcf.gz
bcftools mpileup --threads 20 -Ou -f $REF --max-depth 500 -b bam.list --annotate FORMAT/AD -r NC_020527.1 | bcftools call -mv -Oz -o sino_mel_NC_020527.1.vcf.gz

ls sino_mel_NC* > sino_mel_list

bcftools concat -f sino_mel_list -o sino_mel_concat.vcf.gz


### IGNORE BELOW THIS LINE, Unless it's useful ###



#sbatch -N 1 -n 20 --mem=500GB -t 48:00:00 -p open pileup.sh 

#rename chroms
samtools idxstats A16_O_Fa20_1_S186.markdup.bam # get chr names
# NC_020528.1	
# NC_020527.1	
# NC_020560.1

sed -i -e 's/NC_020528.1/chr/g' sino_mel_concat.vcf
sed -i -e 's/NC_020527.1/psyma/g' sino_mel_concat.vcf
sed -i -e 's/NC_020560.1/psymb/g' sino_mel_concat.vcf

#zip and index
bgzip -c sino_mel_concat.vcf > sino_mel_concat.vcf.gz
tabix sino_mel_concat.vcf.gz

#filter by allele freq and chrom
bcftools view -r chr -q 0.01:minor sino_mel_concat.vcf.gz -o sino_mel_001_chr.vcf.gz 
bcftools view -r psyma -q 0.01:minor sino_mel_concat.vcf.gz -o sino_mel_001_psyma.vcf.gz 
bcftools view -r psymb -q 0.01:minor sino_mel_concat.vcf.gz -o sino_mel_001_psymb.vcf.gz 

bcftools view -r chr -q 0.05:minor sino_mel_concat.vcf.gz -o sino_mel_005_chr.vcf.gz 
bcftools view -r psyma -q 0.05:minor sino_mel_concat.vcf.gz -o sino_mel_005_psyma.vcf.gz 
bcftools view -r psymb -q 0.05:minor sino_mel_concat.vcf.gz -o sino_mel_005_psymb.vcf.gz 


#samples_to_remove = "L" samples
vcftools --remove-indv -f samples_to_remove --vcf sino_mel_concat.vcf.gz --recode --out sino_mel_filtered

bgzip -c sino_mel_filtered.recode.vcf > sino_mel_filtered.recode.vcf.gz
tabix sino_mel_filtered.recode.vcf.gz

bcftools view -r chr -q 0.01:minor sino_mel_filtered.recode.vcf.gz -o sino_mel_001_chr.vcf.gz 
bcftools view -r psyma -q 0.01:minor sino_mel_filtered.recode.vcf.gz -o sino_mel_001_psyma.vcf.gz 
bcftools view -r psymb -q 0.01:minor sino_mel_filtered.recode.vcf.gz -o sino_mel_001_psymb.vcf.gz 

bcftools view -r chr -q 0.05:minor sino_mel_filtered.recode.vcf.gz -o sino_mel_005_chr.vcf.gz 
bcftools view -r psyma -q 0.05:minor sino_mel_filtered.recode.vcf.gz -o sino_mel_005_psyma.vcf.gz 
bcftools view -r psymb -q 0.05:minor sino_mel_filtered.recode.vcf.gz -o sino_mel_005_psymb.vcf.gz 
bcftools view -q 0.05:minor sino_mel_filtered.recode.vcf.gz -o sino_mel_filtered_005.vcf.gz 

