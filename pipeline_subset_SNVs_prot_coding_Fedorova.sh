#!/bin/bash

#SBATCH --job-name=tsv_processing
#SBATCH -D '.'
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --ntasks=4
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4


# Subset protein coding variants and sort 
cat mutations_Fedorova_test4.tsv | grep -E 'protein_coding' | sort -k1,1V -k2,2n -k3,3n > tmp_subset_variants_protein_coding_Fedorova_test4.tsv

# Add header
cat <(head -1 mutations_Fedorova_test4.tsv) tmp_subset_variants_protein_coding_Fedorova_test4.tsv > subset_variants_protein_coding_Fedorova_test4.tsv

# Select only SNVs
cat subset_variants_protein_coding_Fedorova_test4.tsv | awk '{if(length($4)==1 && length($5)==1) print}' > tmp_subset_snvs_protein_coding_Fedorova_test4.tsv
cat <(head -1 mutations_Fedorova_test4.tsv) tmp_subset_snvs_protein_coding_Fedorova_test4.tsv > subset_snvs_protein_coding_Fedorova_test4.tsv
   
# Remove intermediate files
rm subset_variants_protein_coding_Fedorova_test4.tsv
rm tmp_subset_variants_protein_coding_Fedorova_test4.tsv
rm tmp_subset_snvs_protein_coding_Fedorova_test4.tsv