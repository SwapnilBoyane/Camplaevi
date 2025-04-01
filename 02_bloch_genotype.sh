#!/bin/bash
#SBATCH --chdir=.
#SBATCH --job-name=align
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=8
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-30

source activate samtools

# define main working directory\
workdir=/lustre/scratch/sboyane/camplaevi/01_blochmannia

# base name of fastq files, intermediate files, and output files\
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames.txt | tail -n1 )

# define the reference genome\
refgenome=/lustre/work/sboyane/ref/bloch/C145_quercicola.fasta

# run bwa mem
/home/sboyane/anaconda3/bin/bwa mem -t 8 ${refgenome} \
${workdir}/01_cleaned/${basename_array}_R1.fastq.gz \
${workdir}/01_cleaned/${basename_array}_R2.fastq.gz > \
${workdir}/01_bam_files/${basename_array}.sam

# filter for mapped reads\

/home/sboyane/anaconda3/bin/samtools view -b -f 2 \
-o ${workdir}/01_bam_files/${basename_array}.bam \
${workdir}/01_bam_files/${basename_array}.sam

# remove sam
rm ${workdir}/01_bam_files/${basename_array}.sam

# clean up the bam file
~/anaconda3/bin/picard \
CleanSam I=${workdir}/01_bam_files/${basename_array}.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned.bam

# remove the raw bam
rm ${workdir}/01_bam_files/${basename_array}.bam

# sort the cleaned bam file
~/anaconda3/bin/picard \
SortSam I=${workdir}/01_bam_files/${basename_array}_cleaned.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam SORT_ORDER=coordinate

# remove the cleaned bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned.bam

# add read groups to sorted and cleaned bam file
~/anaconda3/bin/picard \
AddOrReplaceReadGroups I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam \
RGLB=1 RGPL=illumina RGPU=unit1 RGSM=${basename_array}

# remove cleaned and sorted bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam


# remove duplicates to sorted, cleaned, and read grouped bam file (creates final bam file)
~/anaconda3/bin/picard \
MarkDuplicates REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 \
M=${workdir}/01_bam_files/${basename_array}_markdups_metric_file.txt \
I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam \
O=${workdir}/01_bam_files/${basename_array}_final.bam

# remove sorted, cleaned, and read grouped bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam

# index the final bam file
/home/sboyane/anaconda3/bin/samtools index ${workdir}/01_bam_files/${basename_array}_final.bam
