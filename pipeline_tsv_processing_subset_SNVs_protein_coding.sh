#!/bin/bash

#SBATCH --job-name=tsv_processing
#SBATCH -D '.'
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --ntasks=10
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8


# Subset missense variants and sort --> ens estem quedant només amb SNVs
cd /slgpfs/scratch/cli79/cli79334/projects/other/1kGPhg38/results/20230804_200606_pipeline_bam_targetVariantCallingRed_coding_1kGPhg38/

#cat mutations_coding_1kGPhg38.tsv | grep -E 'protein_coding' | sort -k1,1V -k2,2n -k3,3n > tmp_subset_variants_protein_coding_1kGPhg38.tsv



# Add header
#cd /slgpfs/scratch/cli79/cli79334/projects/other/1kGPhg38/results/20230804_200606_pipeline_bam_targetVariantCallingRed_coding_1kGPhg38/

#cat <(head -1 mutations_coding_1kGPhg38.tsv) tmp_subset_variants_protein_coding_1kGPhg38.tsv > subset_variants_protein_coding_1kGPhg38.tsv


# Only snvs
#$4 --> REF
#$5 --> ALT
#cd /slgpfs/scratch/cli79/cli79334/projects/other/1kGPhg38/results/20230804_200606_pipeline_bam_targetVariantCallingRed_coding_1kGPhg38/

#cat subset_variants_protein_coding_1kGPhg38.tsv | awk '{if(length($4)==1 && length($5)==1) print}' > tmp_subset_snvs_protein_coding_1kGPhg38.tsv
cat <(head -1 mutations_coding_1kGPhg38.tsv) tmp_subset_snvs_protein_coding_1kGPhg38.tsv > subset_snvs_protein_coding_1kGPhg38.tsv
# snvs --> 564752
# indels --> 79    

# Remove intermediate files
cd /slgpfs/scratch/cli79/cli79334/projects/other/1kGPhg38/results/20230804_200606_pipeline_bam_targetVariantCallingRed_coding_1kGPhg38/

rm subset_variants_protein_coding_1kGPhg38.tsv
rm tmp_subset_variants_protein_coding_1kGPhg38.tsv
rm tmp_subset_snvs_protein_coding_1kGPhg38.tsv




#https://pcingola.github.io/SnpEff/se_inputoutput/#ann-field-vcf-output-files

# Count each category
cat subset_snvs_protein_coding_1kGPhg38.tsv | sed '1d' |  cut -f11 | sort | uniq -c > count_effects.txt