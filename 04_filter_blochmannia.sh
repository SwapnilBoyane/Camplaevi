#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G


source activate vcftools

# define main working directory\
workdir=/lustre/scratch/sboyane/camplaevi/01_blochmannia


#run vcftools with SNP output spaced 20kbp\
#for PCA, EEMS, IBD \
#Filter 1 (in 04_vcf folder)\

vcftools --vcf ${workdir}/04_vcf/blochmannia_camplaevi.vcf  --max-missing 1.0 --keep keeplist_pca.txt  --min-alleles 2 --max-alleles 2 --mac 2  --max-maf 0.49 \
 --recode --recode-INFO-all --out ${workdir}/05_filtered_vcf/bloch_camplaevi_20kbp_final

bgzip ${workdir}/05_filtered_vcf/bloch_camplaevi_20kbp_final.vcf
tabix ${workdir}/05_filtered_vcf/bloch_camplaevi_20kbp_final.vcf.gz

