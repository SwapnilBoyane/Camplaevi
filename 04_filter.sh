#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-31

conda activate bcftools
conda activate vcftools

# define main working directory
workdir=/lustre/scratch/sboyane/camplaevi

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

#run vcftools with SNP output spaced 20kbp
#for ADMIXTURE, PCA, EEMS, IBD 
vcftools --vcf ${workdir}/04_vcf/${region_array}.g.vcf --keep keeplist.txt --max-missing 1.0   --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 2 --thin 20000 --max-maf 0.49  --recode --recode-INFO-all --out ${workdir}/05_20kbp_vcf/${region_array}

# run vcftools with SNP output spaced 20kbp 
#for relatedness 
vcftools --vcf ${workdir}/04_vcf/${region_array}.g.vcf --keep keeplist.txt --max-missing 1.0   --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 1 --thin 20000 --max-maf 0.49  --recode --recode-INFO-all --out ${workdir}/04_related_vcf_20kbp/${region_array}

# run vcftools with SNP and invariant site output, 20% max missing data, no indels 
#for Phylogeny 
vcftools --vcf ${workdir}/04_vcf/${region_array}.g.vcf --max-missing 0.8 --minGQ 20 --minDP 6 --max-meanDP 50 --max-alleles 2  --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/phylo/${region_array}

# run vcftools with SNP and invariant site output, 20% max missing data, no indels
#for observed heterozygosity remove outgroup
vcftools --vcf ${workdir}/04_vcf/${region_array}.g.vcf --keep keeplist.txt --max-missing 0.8 --minGQ 20 --minDP 6 --max-meanDP 50 --max-alleles 2  --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/08_OH/${region_array}

# bgzip and tabix index files that will be subdivided into windows
# directory 1
~/anaconda3/bin/bgzip ${workdir}/phylo/${region_array}.recode.vcf
#tabix
~/anaconda3/bin/tabix -p vcf ${workdir}/phylo/${region_array}.recode.vcf.gz
