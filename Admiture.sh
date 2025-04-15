#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=admixture
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=16G

source activate admixture
cd /lustre/scratch/sboyane/camplaevi/05_filtered_vcf/admixture

# make one vcf
grep "#" scaffold0001.recode.vcf > admix_20kbp_total.vcf

#removing campomotus modoc from the vcf list
for i in $( ls *recode.vcf); do grep -v "#" $i >> admix_20kbp_total.vcf; done

# make chromosome map for the vcf
grep -v "#" admix_20kbp_total.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt

# run vcftools for the combined vcf
vcftools --vcf admix_20kbp_total.vcf  --plink --chrom-map chrom_map.txt --out total

# convert  with plink
plink --file total --recode12 --out total2 --allow-extra-chr

# run admixture
for K in 1 2 3 4 5; do admixture --cv total2.ped $K  | tee log_${K}.out; done

# check cv
grep -h CV log_*.out

