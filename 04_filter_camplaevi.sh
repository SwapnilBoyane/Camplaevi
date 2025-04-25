#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-31


conda activate vcftools

# define main working directory
workdir=/lustre/scratch/sboyane/camplaevi

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

# filter data for ADMIXTURE, PCA, EEMS, IBD. 
#pca20kbp
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep keeplist.txt  --max-missing 1.0 --min-alleles 2 --max-alleles 2 --max-maf 0.49 --mac 2 --thin 20000 --remove-indels --recode --recode-INFO-all --out ${workdir}/05_filtered_vcf/${region_array}

# run vcftools with SNP and invariant site output, 20% max missing data, no indels
#for phylogeny and observed heterozygosity remove outgroup
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep keeplist.txt --max-missing 0.8 --max-alleles 2  --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/08_OH/${region_array}

# run bcftools to simplify the vcftools output for observed heterozygosity analysis
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/08_OH/${region_array}.recode.vcf > ${workdir}/08_OH/${region_array}.simple.vcf


# compress files with bgzip and index with tabix
~/anaconda3/bin/bgzip ${workdir}/05_filtered_vcf/${region_array}.recode.vcf
~/anaconda3/bin/tabix ${workdir}/05_filtered_vcf/${region_array}.recode.vcf.gz
