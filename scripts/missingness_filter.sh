#!/bin/bash

# Compress the VCF (overwrite original or rename as desired)
bgzip -f /home/paula/rhizobium/rhizobium_trapping/data/documentation-test.vcf

# Index the compressed VCF
tabix -p vcf /home/paula/rhizobium/rhizobium_trapping/data/documentation-test.vcf.gz

# Filter by missingness
bcftools +fill-tags /home/paula/rhizobium/rhizobium_trapping/data/documentation-test.vcf.gz \
  -Oz -o /home/paula/rhizobium/rhizobium_trapping/data/documentation-test.filled.vcf.gz -- -t F_MISSING

tabix -p vcf /home/paula/rhizobium/rhizobium_trapping/data/documentation-test.filled.vcf.gz

bcftools view -i 'F_MISSING<=0.3' \
  /home/paula/rhizobium/rhizobium_trapping/data/documentation-test.filled.vcf.gz \
