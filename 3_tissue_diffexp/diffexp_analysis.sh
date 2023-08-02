#!/usr/bin/bash

########################################
# Script name: diffexp_analysis.sh
# Script purpose: Run differential expression R-scripts
# Author: Jessica A. Goodheart
# Last edited: 30 May 2023
# Usage: bash diffexp_analysis.sh
########################################

cld=[PATH TO]/tissue_diffexp/deseq2_clades
tid=[PATH TO]/tissue_diffexp/deseq2_tissues

######### Clade-specific analyses  ################
# Create All_other_genes.txt file for all non-clade-specific genes
cat "$cld"/berghia/Berghia_stephanieae_genes.txt "$cld"/aeolidoidea/Aeolidioidea_only_genes.txt "$cld"/nudibranchia/Nudibranchia_only_genes.txt "$cld"/gastropoda/Gastropoda_only_genes.txt "$cld"/mollusca/Mollusca_only_genes.txt | sort | uniq > clade_specific_genes.txt
comm -23 BsV1_geneids.txt clade_specific_genes.txt > "$cld"/all_others/All_other_genes.txt

# Run DESeq2 setup analysis
cd "$cld"
Rscript Bs_genome_DESeq2.R

# Loop through R-scripts
clades=(berghia aeolidoidea nudibranchia gastropoda mollusca all_others all_genes)

for i in "${clades[@]}"; do
    cd "$cld"/"$i"
    Rscript Bs_genome_DESeq2_$i.R 
done

######### Tissue expression analyses  #############
# Loop through R-scripts
tissues=(distceras foot oraltentacle rhinophore tail proximalceras brain)

for i in "${tissues[@]}"; do
    cd "$tid"/"$i"
    Rscript Bs_genome_DESeq2_$i.R 
done

############# Post-DEseq2 analyses  ###############
cd "$tid"
Rscript Bs_genome_DESeq2_tissues.R

cd "$cld"
Rscript Bs_genome_clade_specific_genes.R
