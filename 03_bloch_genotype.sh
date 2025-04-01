#!/bin/bash
#SBATCH --chdir=.
#SBATCH --job-name=genotype
#SBATCH --partition=nocona
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-30

# Prepend bcftools path to PATH
source activate vcftools

# define main working directory\
workdir=/lustre/scratch/sboyane/camplaevi/01_blochmannia

# base name of fastq files, intermediate files, and output files\
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames.txt | tail -n1 )

# define the reference genome\
refgenome=/lustre/work/sboyane/ref/bloch/C145_quercicola.fasta

# run bcftools to genotype
bcftools mpileup --skip-indels -C 0 -d 200 --min-MQ 10 --threads 4 -f ${refgenome} ${workdir}/01_bam_files/${basename_array}_final.bam | bcftools call --ploidy 1 -m --threads 4 -o ${workdir}/02_vcf/${basename_array}.vcf

# bgzip
~/anaconda3/bin/bgzip ${workdir}/02_vcf/${basename_array}.vcf

#tabix
~/anaconda3/bin/tabix ${workdir}/02_vcf/${basename_array}.vcf.gz

# filter individual vcf files
bcftools view -i 'MIN(DP)>5' ${workdir}/02_vcf/${basename_array}.vcf.gz > ${workdir}/03_vcf/${basename_array}.vcf
# bgzip
~/anaconda3/bin/bgzip ${workdir}/03_vcf/${basename_array}.vcf

#tabix
~/anaconda3/bin/tabix ${workdir}/03_vcf/${basename_array}.vcf.gz
