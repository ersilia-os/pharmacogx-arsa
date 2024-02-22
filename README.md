# PharmacoGx ARSA
Look for abundant region specific alleles for pharmacogenomics

This repository is related to [pharmacogx-embeddings](https://github.com/ersilia-os/pharmacogx-embeddings). Here, a set of scripts and executable notebooks are offered to map Africa-abundant and Africa-specific variants across the 1000 Genomes Project, PharmGKB and others. The repository is related to Anna Montaner's MSc Thesis.

The most important output files are the following (found in `results/` folder):
* `abundance_prec_recall_curves.tsv` and `specificity_prec_recall_curves.tsv`: These files are used to find the ideal cutoffs for abundance (AFR 20%) and specificity (10x). 
* `subset_snvs_protein_coding_1kGPhg38_gene_level.tsv`: SNVs found in the 1000 Genomes Project, grouped by gene. For each gene, the number of Africa-abundant and Africa-specific mutations are counted. SNVs are also stratified by "intron", "missense", or "other".
* `pharmgkb_mutations_gene_level.tsv`: Mutations found in PharmGKB, annotated as described above.

These are the only files that are eventually used in the [pharmacogx-embeddings](https://github.com/ersilia-os/pharmacogx-embeddings) repository. The rest can be used for consultation and reproducibility.