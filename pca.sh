#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=pca
#SBATCH --partition=nocona
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G


cd /lustre/scratch/sboyane/camplaevi/05_filtered_vcf/20kbp_pca

grep "#" scaffold0001.recode.vcf > pca_20kbp.all.vcf

for i in $( ls *recode.vcf ); do
grep -v "#" $i >> pca_20kbp.all.vcf;
done


#for PCA 
plink --vcf pca_20kbp.all.vcf --allow-extra-chr --double-id --set-missing-var-ids @:# --pca --out pca_20kbp 
