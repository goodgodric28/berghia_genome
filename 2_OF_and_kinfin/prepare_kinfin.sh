#!/usr/bin/bash

########################################
# Script name: prepare_kinfin.sh
# Script purpose: Prepare input files to kinfin inside kinfin directory
# Author: Jessica A. Goodheart
# Last edited: 28 June 2022
# Usage: bash prepare_kinfin.sh [path to orthofinder run] [desired output directory]
# Example bash prepare_kinfin.sh /expanse/lustre/projects/sio139/sluglife/crepidula_orthofinder/OrthoFinder/Results_Jun28 crepidula_orthofinder_June2022
########################################

# Set up environment
source /expanse/lustre/projects/sio139/shared/bashrc_files/kinfin.bashrc

# Set path variable
of=$1
out=$2

# Copy SequenceIDs and SpeciesIDs files
echo "Copying OrthoFinder files (SpeciesIDs.txt, SequenceIDs.txt, SpeciesTree_rooted.txt, Orthogroups.tsv)..."
cp $of/WorkingDirectory/SpeciesIDs.txt .
cp $of/WorkingDirectory/SequenceIDs.txt .

# Copy species tree
cp $of/Species_Tree/SpeciesTree_rooted.txt .

# Copy and edit OrthoFinder Orthogroups file for input to KinFin
cp $of/Orthogroups/Orthogroups.tsv .

echo "Preparing Orthogroups_edited.txt file..."
sed 's/\t\t/ /g' Orthogroups.tsv | sed 's/\t/ /g' | sed 's/,//g' | sed 's/^M$//g' | sed -r 's/(OG[0-9][0-9][0-9][0-9][0-9][0-9][0-9])\ +/\1\ /g' > Orthogroups_edited.txt

dos2unix Orthogroups_edited.txt

tail -n +2 Orthogroups_edited.txt > Orthogroups_edited.txt.temp
mv Orthogroups_edited.txt.temp Orthogroups_edited.txt

# Prepare config file
echo "Preparing config.txt file..."
num=($(awk -F': ' '{print $1}' SpeciesIDs.txt))
sp=($(awk '{print substr($2, 1, length($2)-4)}' SpeciesIDs.txt))

declare -a taxids
for i in "${sp[@]}";
do
    if [ $i == "Berghia_stephanieae_brakerRMplusISO" ]; then
	i="Berghia_stephanieae"
	taxid=$(python get_species_taxonomy.py $i)
        taxids+=( $taxid )
    elif [ $i == "Unidentia_sp" ]; then
	i="Unidentia_sp_JG-2016"
	taxid=$(python get_species_taxonomy.py $i)
        taxids+=( $taxid )
    else
	taxid=$(python get_species_taxonomy.py $i)
	taxids+=( $taxid )
    fi
done

printf '#IDX,TAXON,TAXID\n' > config.txt
for i in "${!taxids[@]}"; do
    printf "%s,%s,%s\n" "${num[i]}" "${sp[i]}" "${taxids[i]}" >> config.txt
done

# Prepare SLURM script 
echo "Preparing submit_kinfin.sh file..."
echo "#!/bin/sh                                                  

#SBATCH --job-name=kinfin
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=large-shared
#SBATCH --mem=200GB
#SBATCH --time=12:00:00           
#SBATCH --account=sio139
#SBATCH --output=%x.%j.out

date

kinfin -g Orthogroups_edited.txt -c config.txt -s SequenceIDs.txt -p SpeciesIDs.txt -f functional_annotation.txt -a ../pep_files -t SpeciesTree_rooted.txt -o $out --infer_singletons -r phylum,class,order,superfamily

date" > submit_kinfin.sh
echo "Finished."
