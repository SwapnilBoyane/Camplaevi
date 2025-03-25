#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

cd /lustre/scratch/sboyane/camplaevi/03_vcf

###total genotyped sites per individual

for file in *.gz; do
  count=$(zcat "$file" | wc -l)
  echo "$file: $count" >> total_sites_all.txt
done
