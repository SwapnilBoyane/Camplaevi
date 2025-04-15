

camponotus_laevigatus phylogeography
# Camplaevi
1. set up directories: 
2. Run "01_align_camplaevi.sh" to align host fastq files.
3. Run "02_genotype.sh" to genotype host genomes.
4. Run "02a_genotype_site_per_indivi.sh" to conut total genotyped sites.
5. Run "03_merge.sh"to combines individual VCFs by scaffold
6. Run "04_filter_camplaevi.sh" to filter SNPs for downstream analysis

#Blochmannia pipeline 
1. Run "01_extract_blochmannia.sh" to extraxct blochmannia reads.
2. Run "01a_bloch_align.sh" to align blochmannia reads.
3. Run "02_bloch_genotype.sh" to genotype. Note: "--ploidy 1" was mentioned during genotyping.
4. Run "03_merge_vcf_bloch_camplaevi.sh" to merge.
5. Run "04_filter_blochmannia.sh" to filter SNPs for downstream analysis.
