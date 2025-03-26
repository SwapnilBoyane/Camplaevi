#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam_processing     # Names the job "bam_processing” #
#SBATCH --partition nocona      #Specifies the compute partition (queue) where the job will run#
#SBATCH --nodes=1 --ntasks=8
#SBATCH --time=48:00:00   # Allocates 48 hours of runtime#
#SBATCH --mem-per-cpu=8G # Allocates 8GB of RAM per CPU
#SBATCH --array=1-30 # Runs the job as an array from 1 to 30. This means 30 different jobs will be processed, each with a different dataset#

#use module to load tools or use direct file path
source activate samtools

# define main working directory #Sets the main directory where input files are stored
workdir=/lustre/scratch/<eraider>/camplaevi

#define file for array job  #Reads the nth line (matching SLURM_ARRAY_TASK_ID) from basenames.txt
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames.txt | tail -n1 )

# define the reference genome # define variable meaning that whenever we use "{refgenome}" it will call for the path define below
refgenome=/lustre/work/<eraider>/ref/camp_sp_genome_filtered.fasta

# define the location of the reference mitogenomes
mito=/home/jmanthey/denovo_genomes/formicinae_mitogenomes.fasta

# define the location of the reference blochmannia genomes
bloch=/lustre/work/<eraider>/ref/bloch/C019_modoc.fasta

# run bbduk ###Uses BBDuk (BBTools) to: Remove sequencing adapters, Trim low-quality bases (qtrim=rl trimq=10), Remove short reads (minlen=50), Filter contaminants using adapters.fa
/lustre/work/jmanthey/bbmap/bbduk.sh in1=${workdir}/00_fastq/${basename_array}_R1.fastq.gz in2=${workdir}/00_fastq/${basename_array}_R2.fastq.gz out1=${workdir}/01_cleaned/${basename_array}_R1.fastq.gz out2=${workdir}/01_cleaned/${basename_array}_R2.fastq.gz minlen=50 ftl=10 qtrim=rl trimq=10 ktrim=r k=25 mink=7 ref=/lustre/work/jmanthey/bbmap/resources/adapters.fa hdist=1 tbo tpe

# run bbsplit mitogenomes ##BBSplit separates mitochondrial reads from the cleaned data
/lustre/work/jmanthey/bbmap/bbsplit.sh in1=${workdir}/01_cleaned/${basename_array}_R1.fastq.gz in2=${workdir}/01_cleaned/${basename_array}_R2.fastq.gz ref=${mito} basename=${workdir}/01_mtDNA/${basename_array}_%.fastq.gz outu1=${workdir}/01_mtDNA/${basename_array}_R1.fastq.gz outu2=${workdir}/01_mtDNA/${basename_array}_R2.fastq.gz

# remove unnecessary bbsplit output files
rm ${workdir}/01_mtDNA/${basename_array}_R1.fastq.gz
rm ${workdir}/01_mtDNA/${basename_array}_R2.fastq.gz

# run bwa mem ##Uses BWA-MEM to align reads to the reference genome# in this step it uses the reference genome and align the raw reads and output will be saved in  01_bam_files directory
/home/<eraider>/anaconda3/bin/bwa mem -t 8 ${refgenome} ${workdir}/01_cleaned/${basename_array}_R1.fastq.gz ${workdir}/01_cleaned/${basename_array}_R2.fastq.gz > ${workdir}/01_bam_files/${basename_array}.sam

# convert sam to bam #This step Converts the SAM file (text format) to BAM (compressed binary format)
~/anaconda3/envs/samtools/bin/samtools view -b -S -o ${workdir}/01_bam_files/${basename_array}.bam ${workdir}/01_bam_files/${basename_array}.sam

# remove sam
rm ${workdir}/01_bam_files/${basename_array}.sam

# clean up the bam file # used Picard.jar for sorting and cleaning bam 
~/anaconda3/bin/picard CleanSam I=${workdir}/01_bam_files/${basename_array}.bam O=${workdir}/01_bam_files/${basename_array}_cleaned.bam

# remove the raw bam
rm ${workdir}/01_bam_files/${basename_array}.bam

# sort the cleaned bam file #Sorts reads by genome position
~/anaconda3/bin/picard SortSam I=${workdir}/01_bam_files/${basename_array}_cleaned.bam O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam SORT_ORDER=coordinate

# remove the cleaned bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned.bam

# add read groups to sorted and cleaned bam file
~/anaconda3/bin/picard AddOrReplaceReadGroups I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam RGLB=1 RGPL=illumina RGPU=unit1 RGSM=${basename_array}

# remove the raw bam
rm ${workdir}/01_bam_files/${basename_array}.bam

# sort the cleaned bam file #Sorts reads by genome position
~/anaconda3/bin/picard SortSam I=${workdir}/01_bam_files/${basename_array}_cleaned.bam O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam SORT_ORDER=coordinate

# remove the cleaned bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned.bam

# add read groups to sorted and cleaned bam file
~/anaconda3/bin/picard AddOrReplaceReadGroups I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam RGLB=1 RGPL=illumina RGPU=unit1 RGSM=${basename_array}

# remove cleaned and sorted bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam

# remove duplicates to sorted, cleaned, and read grouped bam file (creates final bam file) #Removes PCR duplicates (to avoid false positives in variant calling.
~/anaconda3/bin/picard MarkDuplicates REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 M=${workdir}/01_bam_files/${basename_array}_markdups_metric_file.txt I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam O=${workdir}/01_bam_files/${basename_array}_final.bam

# remove sorted, cleaned, and read grouped bam file #Creates a BAM index (.bai file) to enable faster access
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam

# index the final bam file
~/anaconda3/envs/samtools/bin/samtools index ${workdir}/01_bam_files/${basename_array}_final.bam
