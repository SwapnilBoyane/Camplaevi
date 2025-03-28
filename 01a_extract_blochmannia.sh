#!/bin/bash
#SBATCH --chdir=.
#SBATCH --job-name=extract
#SBATCH --partition= nocona
#SBATCH --nodes=1 --ntasks=12
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-30

source activate samtools


# define main working directory
workdir=/lustre/scratch/sboyane/camplaevi

#define file for array job  #Reads the nth line (matching SLURM_ARRAY_TASK_ID) from basenames.txt
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames.txt | tail -n1 )

# define the reference genome
ref=/lustre/work/sboyane/ref/bloch/C145_quercicola.fasta 

# run bbsplit blochmannia  ##BBSplit isolates Blochmannia bacterial reads from cleaned fastq files
/lustre/work/jmanthey/bbmap/bbsplit.sh in1=${workdir}/01_cleaned/${basename_array}_R1.fastq.gz in2=${workdir}/01_cleaned/${basename_array}_R2.fastq.gz ref=${bloch} basename=${workdir}/01_blochmannia/${basename_array}_%.fastq.gz outu1=${workdir}/01_blochmannia/${basename_array}_R1.fastq.gz outu2=${workdir}/01_blochmannia/${basename_array}_R2.fastq.gz

# remove unnecessary bbsplit output files
rm ${workdir}/01_blochmannia/${basename_array}_R1.fastq.gz
rm ${workdir}/01_blochmannia/${basename_array}_R2.fastq.gz

# Extract the reads from the interleaved file and separate them into R1 and R2 files
seqtk seq -1 ${workdir}/01_blochmannia/${basename_array}_C145_quercicola.fastq.gz > ${workdir}/01_blochmannia/01_cleaned/${basename_array}_R1.fastq
seqtk seq -2 ${workdir}/01_blochmannia/${basename_array}_C145_quercicola.fastq.gz > ${workdir}/01_blochmannia/01_cleaned/${basename_array}_R2.fastq

# Compress the files to .gz
gzip ${workdir}/01_blochmannia/01_cleaned/${basename_array}_R1.fastq
gzip ${workdir}/01_blochmannia/01_cleaned/${basename_array}_R2.fastq



