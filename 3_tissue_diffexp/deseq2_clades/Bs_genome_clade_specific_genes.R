####################################
# Bs_genome_clade_specific_genes.R
# Written by: Jessica A. Goodheart
# Last Updated: 26 July 2022
# Purpose: To analyze the clade and tissue distribution of genes in Berghia
# Inputs used: Counts from htseq-counts, Orthofinder/kinfin results, and functional annotations
####################################

####################################
# Initial setup
####################################
# Packages
require(ggplot2)
require(dplyr)
require(stringr)

# Set the working directory
directory = "/Users/jessicagoodheart/Dropbox/Research/2_Projects/berghia_genome/April_2021/annotations_analysis/tissue_diffexp/deseq2/"
setwd(directory)

# Pull in counts data for all genes
all_genes = read.csv("Bs_tissues_DESeq2-full-normalized-counts.csv")
all_genes$gene = as.character(all_genes$gene)

# Pull in clade-specific lists of genes
berghia = scan("berghia/Berghia_stephanieae_genes.txt", what="", sep="\n")
aeolidina = scan("aeolidina/Aeolidina_only_genes.txt", what="", sep="\n")
nudibranchia = scan("nudibranchia/Nudibranchia_only_genes.txt", what="", sep="\n")
gastropoda = scan("gastropoda/Gastropoda_only_genes.txt", what="", sep="\n")
mollusca = scan("mollusca/Mollusca_only_genes.txt", what="", sep="\n")
other = scan("all_others//All_other_genes.txt", what="", sep="\n")

# Pull in tissue-specific lists of genes
BR_only = scan("all_genes/gene_annotations/Bs_tissues_allgenes_DESeq2-BR-only-genelist.txt", what="", sep="\n")
OT_only = scan("all_genes/gene_annotations/Bs_tissues_allgenes_DESeq2-OT-only-genelist.txt", what="", sep="\n")
RH_only = scan("all_genes/gene_annotations/Bs_tissues_allgenes_DESeq2-RH-only-genelist.txt", what="", sep="\n")
DC_only = scan("all_genes/gene_annotations/Bs_tissues_allgenes_DESeq2-DC-only-genelist.txt", what="", sep="\n")
PC_only = scan("all_genes/gene_annotations/Bs_tissues_allgenes_DESeq2-PC-only-genelist.txt", what="", sep="\n")
FO_only = scan("all_genes/gene_annotations/Bs_tissues_allgenes_DESeq2-FO-only-genelist.txt", what="", sep="\n")
TA_only = scan("all_genes/gene_annotations/Bs_tissues_allgenes_DESeq2-TA-only-genelist.txt", what="", sep="\n")

# Pull in lists of genes upregulated in each tissue
BR_up = read.csv("../deseq2_tissues/brain/Bs_tissues_brain_DESeq2-upregulated-stats-filtered.csv", header=TRUE)[,1]
OT_up = read.csv("../deseq2_tissues/oraltentacle/Bs_tissues_oraltentacle_DESeq2-upregulated-stats-filtered.csv", header=TRUE)[,1]
RH_up = read.csv("../deseq2_tissues/rhinophore/Bs_tissues_rhinophore_DESeq2-upregulated-stats-filtered.csv", header=TRUE)[,1]
DC_up = read.csv("../deseq2_tissues/distceras/Bs_tissues_distceras_DESeq2-upregulated-stats-filtered.csv", header=TRUE)[,1]
PC_up = read.csv("../deseq2_tissues/proximalceras/Bs_tissues_proximalceras_DESeq2-upregulated-stats-filtered.csv", header=TRUE)[,1]
FO_up = read.csv("../deseq2_tissues/foot/Bs_tissues_foot_DESeq2-upregulated-stats-filtered.csv", header=TRUE)[,1]
TA_up = read.csv("../deseq2_tissues/tail/Bs_tissues_tail_DESeq2-upregulated-stats-filtered.csv", header=TRUE)[,1]

# Pull in lists of genes downregulated in each tissue
BR_down = read.csv("../deseq2_tissues/brain/Bs_tissues_brain_DESeq2-downregulated-stats-filtered.csv", header=TRUE)[,1]
OT_down = read.csv("../deseq2_tissues/oraltentacle/Bs_tissues_oraltentacle_DESeq2-downregulated-stats-filtered.csv", header=TRUE)[,1]
RH_down = read.csv("../deseq2_tissues/rhinophore/Bs_tissues_rhinophore_DESeq2-downregulated-stats-filtered.csv", header=TRUE)[,1]
DC_down = read.csv("../deseq2_tissues/distceras/Bs_tissues_distceras_DESeq2-downregulated-stats-filtered.csv", header=TRUE)[,1]
PC_down = read.csv("../deseq2_tissues/proximalceras/Bs_tissues_proximalceras_DESeq2-downregulated-stats-filtered.csv", header=TRUE)[,1]
FO_down = read.csv("../deseq2_tissues/foot/Bs_tissues_foot_DESeq2-downregulated-stats-filtered.csv", header=TRUE)[,1]
TA_down = read.csv("../deseq2_tissues/tail/Bs_tissues_tail_DESeq2-downregulated-stats-filtered.csv", header=TRUE)[,1]

# Pull in interproscan database accessions
annot_db = read.table("../functional_annotation.txt", header=TRUE, sep="\t")
annot_db$protein_id = gsub(".t[0-9]*", "", as.character(annot_db$protein_id))
annot.db.sub = distinct(annot_db, protein_id, .keep_all=TRUE)

# Pull in blastp annotations
blastp = read.csv("../../genome_annot_nov2021/protein_annotations/filtered/Bs_protein-blasthit-ALL-info.txt", sep="\t", 
                   header=TRUE, as.is=TRUE, na.strings=c("","NA"), fill=TRUE)
blastp$prot_id = gsub(".t[0-9]*", "", as.character(blastp$prot_id))
blastp.sub = distinct(blastp, prot_id, .keep_all=TRUE)
blastp.sub = blastp.sub[,c(1:4,6,10:11)]

####################################
# Create Full Table
####################################

### Clade-specific data

# Create data frame that labels each gene with specific clade or "All Others"
clades = c()
for (i in all_genes$gene) {
  ie = paste("\\b", i, "\\b", sep="")
  index = grep(ie, all_genes$gene)
  if (i %in% berghia) { 
    clades[index]="Berghia"
  } else if (i %in% aeolidina) { 
    clades[index]="Aeolidina"
  } else if (i %in% nudibranchia) { 
    clades[index]="Nudibranchia"
  } else if (i %in% gastropoda) { 
    clades[index]="Gastropoda"
  } else if (i %in% mollusca) { 
    clades[index]="Mollusca"
  } else { 
    clades[index]="Other"
  }
}

# Build full data frame with all clade-specific labels
clade.spec = data.frame(all_genes$gene,clades)
colnames(clade.spec) = c("genes", "clade")
clade.spec$genes = as.character(clade.spec$genes)

### Tissue-specific data

# Create vector with tissue-specific labels or "other"
tissues = c()
for (i in clade.spec$genes) {
  ie = paste("\\b", i, "\\b", sep="")
  index = grep(ie, clade.spec$genes)
  if (i %in% BR_only) { 
    tissues[index]="brain"
  } else if (i %in% OT_only) { 
    tissues[index]="oral tentacle"
  } else if (i %in% RH_only) { 
    tissues[index]="rhinophore"
  } else if (i %in% DC_only) { 
    tissues[index]="distal ceras"
  } else if (i %in% PC_only) { 
    tissues[index]="proximal ceras"
  } else if (i %in% FO_only) { 
    tissues[index]="foot"
  } else if (i %in% TA_only) { 
    tissues[index]="tail" 
  } else { 
    tissues[index]="other"
  }
}

# Add tissue-specific labels to data frame
full.df = data.frame(clade.spec, tissues)
colnames(full.df)[3] <- "tissue.specificity"

### Upregulated in which tissue(s) data

# Create vector with tissue-specific labels or "other"
upreg.df = data.frame(matrix(ncol = 2, nrow = 0))
colnames(upreg.df) = c("genes", "upreg.tissues")
for (i in full.df$genes) {
  ie = paste("\\b", i, "\\b", sep="")
  index = grep(ie, full.df$genes)
  tissues = c()
  if (i %in% BR_up) { 
    tissues = append(tissues, "brain")
  }
  if (i %in% OT_up) { 
    tissues = append(tissues, "oral tentacle")
  } 
  if (i %in% RH_up) { 
    tissues = append(tissues, "rhinophore")
  } 
  if (i %in% DC_up) { 
    tissues = append(tissues, "distal ceras")
  }
  if (i %in% PC_up) { 
    tissues = append(tissues, "proximal ceras")
  }
  if (i %in% FO_up) { 
    tissues = append(tissues, "foot")
  } 
  if (i %in% TA_up) { 
    tissues = append(tissues, "tail")
  } 
  if (length(tissues)==0) { 
    tissues = "none"
  }
  tissues.df = data.frame(i, paste(tissues, collapse="; "))
  colnames(tissues.df) = c("genes", "upreg.tissues")
  upreg.df = rbind(upreg.df, tissues.df)
}

# Add tissue-specific labels to data frame
full.df = data.frame(full.df, upreg.df$upreg.tissues)
colnames(full.df)[4] = "upreg.tissues"

### Upregulated in which tissue(s) data

# Create vector with tissue-specific labels or "other"
downreg.df = data.frame(matrix(ncol = 2, nrow = 0))
colnames(downreg.df) = c("genes", "downreg.tissues")
for (i in full.df$genes) {
  ie = paste("\\b", i, "\\b", sep="")
  index = grep(ie, full.df$genes)
  tissues = c()
  if (i %in% BR_down) { 
    tissues = append(tissues, "brain")
  }
  if (i %in% OT_down) { 
    tissues = append(tissues, "oral tentacle")
  } 
  if (i %in% RH_down) { 
    tissues = append(tissues, "rhinophore")
  } 
  if (i %in% DC_down) { 
    tissues = append(tissues, "distal ceras")
  }
  if (i %in% PC_down) { 
    tissues = append(tissues, "proximal ceras")
  }
  if (i %in% FO_down) { 
    tissues = append(tissues, "foot")
  } 
  if (i %in% TA_down) { 
    tissues = append(tissues, "tail")
  } 
  if (length(tissues)==0) { 
    tissues = "none"
  }
  tissues.df = data.frame(i, paste(tissues, collapse="; "))
  colnames(tissues.df) = c("genes", "downreg.tissues")
  downreg.df = rbind(downreg.df, tissues.df)
}

# Add tissue-specific labels to data frame
full.df = data.frame(full.df, downreg.df$downreg.tissues)
colnames(full.df)[5] = "downreg.tissues"

### Annotations

# Create matching data frames for annotations - Interproscan accessions
ipr.annot = data.frame(matrix(ncol = ncol(annot.db.sub), nrow = 0))
colnames(ipr.annot) = colnames(annot.db.sub)
for (i in full.df$genes) {
  ie = paste("\\b", i, "\\b", sep="")
  s = subset(annot.db.sub, grepl(ie, annot.db.sub$protein_id))
  if(nrow(s)==0) 
    s = data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  ipr.annot = rbind(ipr.annot, s)
}

# Create matching data frames for annotations - BLASTP results
blastp.annot = data.frame(matrix(ncol = ncol(blastp.sub), nrow = 0))
colnames(blastp.annot) = colnames(blastp.sub)
for (i in full.df$genes) {
  ie = paste("\\b", i, "\\b", sep="")
  s = subset(blastp.sub, grepl(ie, blastp.sub$prot_id))
  if(nrow(s)==0) 
    s = data.frame("prot_id"=i, "db_id"=NA, "evalue"=NA, "Entry"=NA,
                   "Protein.names"=NA, "Refseq"=NA, "GO.terms"=NA,stringsAsFactors = FALSE)
  blastp.annot = rbind(blastp.annot, s)
}

# Combine datasets into full data frame
full.df = data.frame(full.df, blastp.annot, ipr.annot)
full.df = full.df[,c(1:5,7:12,14:17)]

# Save table
write.table(full.df, "gene-specificity-annotation-summary.txt", sep = '\t', quote=FALSE, row.names=FALSE, na="NA")
