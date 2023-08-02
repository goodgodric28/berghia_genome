#!/bin/sh                                                  

#SBATCH --job-name=kinfin
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=large-shared
#SBATCH --mem=200GB
#SBATCH --time=12:00:00           
#SBATCH --account=sio139
#SBATCH --output=%x.%j.out

date

kinfin -g Orthogroups_edited.txt -c config.txt -s SequenceIDs.txt -p SpeciesIDs.txt -f functional_annotation.txt -a ../pep_files -t species_tree_jessica_June2022.tre -o berghia_kinfin_June2022 --infer_singletons -r phylum,class,order,superfamily

date
