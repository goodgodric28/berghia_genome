This repository contains the scripts and commands for the pipeline used to construct and process the *Berghia stephanieae* genome and perform associated analyses (Goodheart et al., in prep.; link will be added when available). The commands are  included below, and the submission scripts in this repository include  information on the requested resources for each job.

# Genome Annotation (1_genome_annotation)

Commands that use these scripts are provided in Additional_File_12_Commands_Used_Annotation.docx, which can be found inside the relevant folder here on GitHub and in the supplementary materials of the manuscript.

# OrthoFinder and KinFin Analyses (2_OF_and_kinfin)

## OrthoFinder Analysis
Performed using OrthoFinder version 2.5.4 and MAFFT v7.453:
```
orthofinder -f [PATH TO]/orthofinder -s species_tree_jessica_June2022.tre -t 32 -a 32
```

## KinFin Analysis
Performed using KinFin v1.0:
```
bash prepare_kinfin.sh [PATH_TO]/orthofinder/OrthoFinder/Results_Jun23/ berghia_orthofinder_June2022 

sbatch submit_kinfin.sh

awk -F"\t" '$10 != "" { print $1, $10 }' Orthogroups_UnassignedGenes.tsv > Berghia_UnassignedGenes.tsv

bash post_kinfin.sh Berghia_stephanieae Aeolidoidea Nudibranchia Gastropoda Mollusca
```

# Differential Expression Analyses (3_tissue_diffexp)

## Mapping and Counts
Performed with STAR v2.7.9a and 
HTSeq v1.99.2:
```
STAR --runThreadN 16 --genomeDir berghia_18scaffolds_full --sjdbGTFfile 
augustus.RM_iso.anysupport.gtf --readFilesIn 
OT1_CRRA200005335-1a_HV572DSXX_L1_1.fq.gz 
OT1_CRRA200005335-1a_HV572DSXX_L1_2.fq.gz --outFileNamePrefix 
OT1_augustus_RM_iso_anysupport_starmappedtwopass --readFilesCommand zcat 
--outSAMtype BAM SortedByCoordinate --twopassMode Basic 
--sjdbGTFfeatureExon 'CDS'

htseq-count -s no -r pos -t CDS --add-chromosome-info ${f} 
../OT1_augustus_RM_iso_anysupport_starmappedtwopass.bam
```

## Differential Expression Analysis
Commands must be run for each sample. Performed with R v3.6.2 DESeq2 v1.26.0:
```
bash diffexp_analysis.sh
```

# Citation

Goodheart JA, Rio RA, Taraporevala NF, Fiorenza RA, Barnes SR, Morrill K, Jacob MAC, Whitesel C, Masterson P, Batzel GO, Johnston HT, Ramirez MD, Katz, PS, and Lyons DC. In Preparation. A new chromosome-level genome 
for the nudibranch gastropod *Berghia stephanieae* helps parse clade-specific gene expression in novel and conserved phenotypes.

