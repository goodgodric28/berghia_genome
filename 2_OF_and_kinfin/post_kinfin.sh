#!/usr/bin/bash

########################################
# Script name: post_kinfin.sh
# Script purpose: Obtain lists of clade-specific genes
# Author: Jessica A. Goodheart
# Last edited: 29 June 2022
# Usage: bash post_kinfin.sh [species] [superfamily] [superfamily2] [order] [subclass] [class] [phylum]
# Example: bash post_kinfin.sh Berghia_stephanieae Aeolidoidea Fionoidea Nudibranchia Heterobranchia Gastropoda Mollusca
########################################

# Set up clade names of interest
sp=$1
su=$2
#suF=$3
or=$3
sc=$4
cl=$5
ph=$6

############# Species level ################
# Get list of species-specific clusters
awk -F"\t" '($2=="present" && $3=="specific") { print $1 }' TAXON/TAXON."${sp}"_*.cluster_metrics.txt > "${sp}"_cluster_numbers.txt
awk '{ print $1 }' ../"${sp}"_UnassignedGenes.tsv >> "${sp}"_cluster_numbers.txt

# Get species-specific cluster data
awk 'NR==FNR {a[$1]++; next} $1 in a' "${sp}"_cluster_numbers.txt ../Orthogroups_edited.txt > "${sp}"_clusters.txt
tail -n +2 ../"${sp}"_UnassignedGenes.tsv >> "${sp}"_clusters.txt

# Get lists of species-specific proteins and genes
cat "${sp}"_clusters.txt | tr " " "\n" | grep -o "jg[0-9]*.t[0-9]*\|m1-g[0-9]*.t[0-9]*" | sort | uniq > "${sp}"_proteins.txt
cat "${sp}"_proteins.txt | sed "s/.t[0-9]*//g" | sort | uniq > "${sp}"_genes.txt

######### Superfamily level ################
# Get list of superfamily-specific clusters
awk -F"\t" '($2=="present" && $3=="specific") { print $1 }' superfamily/superfamily."${su}".cluster_metrics.txt  > "${su}"_cluster_numbers.txt
#awk -F"\t" '($2=="present" && $3=="specific") { print $1 }' superfamily/superfamily."${suF}".cluster_metrics.txt >> "${su}"_cluster_numbers.txt

# Get superfamily-specific cluster data
awk 'NR==FNR {a[$1]++; next} $1 in a' "${su}"_cluster_numbers.txt ../Orthogroups_edited.txt > "${su}"_clusters.txt

# Get lists of superfamily-specific proteins and genes
cat "${su}"_clusters.txt | tr " " "\n" | grep -o "jg[0-9]*.t[0-9]*\|m1-g[0-9]*.t[0-9]*" | sort | uniq > "${su}"_proteins.txt
comm -23 "${su}"_proteins.txt "${sp}"_proteins.txt > "${su}"_only_proteins.txt
cat "${su}"_only_proteins.txt | sed "s/.t[0-9]*//g" | sort | uniq > "${su}"_only_genes.txt

############## Order level #################
# Get list of order-specific clusters
awk -F"\t" '($2=="present" && $3=="specific") { print $1 }' order/order."${or}".cluster_metrics.txt > "${or}"_cluster_numbers.txt

# Get order-specific cluster data
awk 'NR==FNR {a[$1]++; next} $1 in a' "${or}"_cluster_numbers.txt ../Orthogroups_edited.txt > "${or}"_clusters.txt

# Get lists of order-specific proteins and genes
cat "${or}"_clusters.txt | tr " " "\n" | grep -o "jg[0-9]*.t[0-9]*\|m1-g[0-9]*.t[0-9]*" | sort | uniq > "${or}"_proteins.txt
comm -23 "${or}"_proteins.txt "${su}"_proteins.txt > tmp
comm -23 tmp "${sp}"_proteins.txt > "${or}"_only_proteins.txt
cat "${or}"_only_proteins.txt | sed "s/.t[0-9]*//g" | sort | uniq > "${or}"_only_genes.txt

############## Subclass level #################
# Get list of class-specific clusters 
awk -F"\t" '($2=="present" && $3=="specific") { print $1 }' subclass/subclass."${sc}".cluster_metrics.txt > "${sc}"_cluster_numbers.txt

# Get subclass-specific cluster data 
awk 'NR==FNR {a[$1]++; next} $1 in a' "${sc}"_cluster_numbers.txt ../Orthogroups_edited.txt > "${sc}"_clusters.txt

# Get lists of subclass-specific proteins and genes  
cat "${sc}"_clusters.txt | tr " " "\n" | grep -o "jg[0-9]*.t[0-9]*\|m1-g[0-9]*.t[0-9]*" | sort | uniq > "${sc}"_proteins.txt
comm -23 "${sc}"_proteins.txt "${or}"_proteins.txt > tmp
comm -23 tmp "${su}"_proteins.txt > tmp.1
comm -23 tmp.1 "${sp}"_proteins.txt > "${sc}"_only_proteins.txt
cat "${sc}"_only_proteins.txt | sed "s/.t[0-9]*//g" | sort | uniq > "${sc}"_only_genes.txt

############## Class level #################
# Get list of class-specific clusters
awk -F"\t" '($2=="present" && $3=="specific") { print $1 }' class/class."${cl}".cluster_metrics.txt > "${cl}"_cluster_numbers.txt

# Get class-specific cluster data
awk 'NR==FNR {a[$1]++; next} $1 in a' "${cl}"_cluster_numbers.txt ../Orthogroups_edited.txt > "${cl}"_clusters.txt

# Get lists of class-specific proteins and genes
cat "${cl}"_clusters.txt | tr " " "\n" | grep -o "jg[0-9]*.t[0-9]*\|m1-g[0-9]*.t[0-9]*" | sort | uniq > "${cl}"_proteins.txt
comm -23 "${cl}"_proteins.txt "${sc}"_proteins.txt > tmp
comm -23 tmp "${or}"_proteins.txt > tmp.1
comm -23 tmp.1 "${su}"_proteins.txt > tmp.2
comm -23 tmp.2 "${sp}"_proteins.txt > "${cl}"_only_proteins.txt
cat "${cl}"_only_proteins.txt | sed "s/.t[0-9]*//g" | sort | uniq > "${cl}"_only_genes.txt

############# Phylum level #################
# Get list of phylum-specific clusters
awk -F"\t" '($2=="present" && $3=="specific") { print $1 }' phylum/phylum."${ph}".cluster_metrics.txt > "${ph}"_cluster_numbers.txt

# Get phylum-specific cluster data
awk 'NR==FNR {a[$1]++; next} $1 in a' "${ph}"_cluster_numbers.txt ../Orthogroups_edited.txt > "${ph}"_clusters.txt

# Get lists of phylum-specific proteins and genes
cat "${ph}"_clusters.txt | tr " " "\n" | grep -o "jg[0-9]*.t[0-9]*\|m1-g[0-9]*.t[0-9]*" | sort | uniq > "${ph}"_proteins.txt
comm -23 "${ph}"_proteins.txt "${cl}"_proteins.txt > tmp
comm -23 tmp "${sc}"_proteins.txt > tmp.1
comm -23 tmp.1 "${or}"_proteins.txt > tmp.2
comm -23 tmp.2 "${su}"_proteins.txt > tmp.3
comm -23 tmp.3 "${sp}"_proteins.txt > "${ph}"_only_proteins.txt
cat "${ph}"_only_proteins.txt | sed "s/.t[0-9]*//g" | sort | uniq > "${ph}"_only_genes.txt

rm tmp*

########### Summary statistics  ############
# Only clade-specific proteins, not including nested proteins
spp=$(cat "${sp}"_proteins.txt | wc -l)
sup=$(cat "${su}"_only_proteins.txt | wc -l)
orp=$(cat "${or}"_only_proteins.txt | wc -l)
scp=$(cat "${sc}"_only_proteins.txt | wc -l)
clp=$(cat "${cl}"_only_proteins.txt | wc -l)
php=$(cat "${ph}"_only_proteins.txt | wc -l)

prots=($spp $sup $orp $scp $clp $php)

# Only clade-specific genes, not including nested genes
spg=$(cat "${sp}"_genes.txt | wc -l)
sug=$(cat "${su}"_only_genes.txt | wc -l)
org=$(cat "${or}"_only_genes.txt | wc -l)
scg=$(cat "${sc}"_only_genes.txt | wc -l)
clg=$(cat "${cl}"_only_genes.txt | wc -l)
phg=$(cat "${ph}"_only_genes.txt | wc -l)

genes=($spg $sug $org $scg $clg $phg)

# Create stats file
printf 'clade\tnum_proteins\tnum_genes\n' > kinfin_stats.txt
a=0
for i in "$@"; do
    printf "%s\t%s\t%s\n" "$i" "${prots[a]}" "${genes[a]}" >> kinfin_stats.txt
    ((a+=1))
done
