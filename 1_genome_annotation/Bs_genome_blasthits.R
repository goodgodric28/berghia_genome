####################################
# Bs_genome_blasthits.R
# Written by: Jessica A. Goodheart
# Last Updated: 10 June 2022
# Purpose: To get UniProt info from genome blast hits for reference.
# Inputs used: UniProt blast hits to protein models from genome
####################################

####################################
# Initial setup
####################################

# Load packages
require(plyr)

# Set the working directory
directory <- "~/Dropbox (Personal)/Research/4Projects/berghia_genome/April_2021/annotations_analysis/genome_annot_nov2021/protein_annotations/filtered"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "Bs_protein"

####################################
# BLASTP hits
####################################
# Pull in blastp hits to UniProt
blastp <- read.table("blastp/blastp.uniprot-sprot.outfmt6", sep="\t")
blastp <- blastp[-c(3:10,12)]
colnames(blastp) <- c("prot_id","uniprot_id","e-value")
colnm <- c("sp", "uniprot_id", "uniprot_code")
blastp_edit <- cbind(blastp$prot_id, reshape2::colsplit(blastp$uniprot_id, "\\|", colnm), blastp$`e-value`)
blastp_edit <- blastp_edit[,-c(2,4)]
blastp_edit <- unique(blastp_edit)
colnames(blastp_edit) <- c("prot_id","db_id","evalue")

# Pull in blastp hits to RefSeq
blastp.ref <- read.table("blastp/blastp.refseq.outfmt6", sep="\t")
blastp.ref <- blastp.ref[-c(3:10,12)]
blastp.ref <- unique(blastp.ref)
colnames(blastp.ref) <- c("prot_id","db_id","evalue")

# Combine 
blastp.all <- rbind(blastp_edit,blastp.ref)
blastp.all.sorted <- blastp.all[order(blastp.all$prot_id),]

####################################
# Annotation data - Uniprot & NCBI
####################################
# Pull in UniProt data for swissprot and refseq dbs
uniprot <- read.csv("uniprot_data.txt", sep="\t", header=TRUE, as.is=TRUE, na.strings=c("","NA"), fill=TRUE)
uniprot <- uniprot[-c(3,10)]
colnames(uniprot)[7:8] <- c('Refseq','GO.terms')

uniprot.ref <- read.csv("refseq_uniprot_data.txt", sep="\t", header=TRUE, as.is=TRUE, na.strings=c("","NA"), fill=TRUE)
uniprot.ref <- uniprot.ref[-c(3,10:11)]
colnames(uniprot.ref)[7:8] <- c('Refseq','GO.terms')

####################################
# Create new annotation dataset
####################################
# Create new database
blastp_full <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(blastp_full) <- c(colnames(blastp_edit),"Protein.names","Gene.names","Organism","Length","Refseq","GO.terms")

missing.db <- vector()
missing.prot <- vector()

for (i in 1:nrow(blastp.all.sorted)) {
  r <- blastp.all.sorted[i,]
  p <- paste("\\b", as.character(blastp.all.sorted$prot_id)[i], "\\b", sep="")
  t <- blastp.all.sorted$db_id[i]
  if (sum(grepl(p, as.character(blastp_full$prot_id)))>0){
    next
  } else {
    s <- subset(uniprot, grepl(t, uniprot$Entry, fixed=TRUE))
    if (nrow(s)>0) {
      full.row <- cbind(r,s)
      colnames(full.row)
      blastp_full <- rbind(blastp_full, full.row)
    } else {
      s <- subset(uniprot.ref, grepl(t, uniprot.ref$Refseq))
      if (nrow(s)>0) {
        full.row <- cbind(r,s)
        colnames(full.row)
        blastp_full <- rbind(blastp_full, full.row)
      } else {
        missing.db <- append(missing.db, t)
        missing.prot <- append(missing.prot, p)
      }
    }
  }
}


# Write out list of blast hit information to file
write.table(blastp_full, paste0(outputPrefix, "-blasthit-uniprot-info.txt"), sep="\t", quote=FALSE, row.names=FALSE)

# Write out list of proteins without blast hits to file
write.table(missing.db, paste0(outputPrefix, "-nouniprothits-info.txt"), sep="\t", quote=FALSE, row.names=FALSE)

# Take Bs_protein_blast-prots-nouniprothits-info.txt file and split into files with ~2000 IDs using 'split -l 2000'
# Upload these files to Batch Entrez (https://www.ncbi.nlm.nih.gov/sites/batchentrez) and search the protein DB
# Download the summary files (as protein_result_[A-Z].txt and combine with 'cat protein_result_*.txt protein_result_allRefSeq.txt')
# Process this file with the script process_batch_entrez.sh, which will produce 'protein_result_allRefSeq_edited_organized_nouncharacterized.txt'
# This file has info for refseq hits that did not match to uniprot IDs and are not considered "uncharacterized"

####################################
# Add RefSeq Annotations
####################################

# Pull in RefSeq data
refseq <- read.csv("protein_result_allRefSeq_edited_organized_nouncharacterized.txt", sep="\t", header=FALSE, as.is=TRUE, na.strings=c("","NA"), fill=TRUE)
colnames(refseq) <- c('Refseq','Gene.names',"Organism","Length")

# Add RefSeq data to blast DB
blastp_refseq <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(blastp_refseq) <- c('Refseq','Gene.names','Organism','Length','prot_id')

for (i in 1:length(missing.db)) {
  p <- paste("\\b", missing.prot[i], "\\b", sep="")
  t <- missing.db[i]
  refs <- subset(blastp.ref, grepl(t, blastp.ref$db_id))
  if (sum(grepl(p, as.character(blastp_full$prot_id)))>0){
    next
  } else if (sum(grepl(p, as.character(blastp_refseq$prot_id)))>0){
    next
  } else {
    s <- subset(refseq, grepl(t, refseq$Refseq))
    if (nrow(s)>0) {
      l <- subset(refs, grepl(p, refs$prot_id))
      sl <- cbind(l[1,],s)
      blastp_refseq <- rbind(blastp_refseq,sl)
    } 
  }
}

# Create full table plus refseq
blastp_full_refseq <- bind_rows(blastp_full,blastp_refseq)

# Write out list of blast hit information to file
write.table(blastp_full_refseq, paste0(outputPrefix, "-blasthit-ALL-info.txt"), sep="\t", quote=FALSE, row.names=FALSE)
