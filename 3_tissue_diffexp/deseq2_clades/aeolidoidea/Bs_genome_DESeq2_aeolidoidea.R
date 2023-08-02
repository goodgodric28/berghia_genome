####################################
# Bs_genome_DESeq2_aeolidoidea.R
# Written by: Jessica A. Goodheart
# Last Updated: 30 May 2023
# Purpose: To analyze differential expression data from Berghia tissues in genes specific to Aeolidoidea.
# Inputs used: Expression counts from htseq-count and Orthofinder/kinfin orthogroup results
####################################

####################################
# Initial setup
####################################

# Set the working directory
directory <- "[PATH TO]/tissue_diffexp/deseq2_clades/aeolidoidea/"
read_dir <- "[PATH TO]/tissue_diffexp/read_files/"
setwd(directory)

# Read in environment from DESeq2 analysis
library("DESeq2")
load("../Bs_genome_DESeq2.Rdata")

# Set the prefix for each output file name
outputPrefix <- "Bs_tissues_aeolidoidea_DESeq2"

# Pull in Berghia gene ids for Aeolidoidea-specific genes from kinfin analysis
aeolidoidea_novel_genes <- scan("Aeolidioidea_only_genes.txt", what="", sep="\n")

# Subset results to only include Aeolidoidea-specific genes and order results by padj value (most significant to least)
res <- resdata[resdata$gene %in% aeolidoidea_novel_genes, ]
res <- res[order(res$padj),]

# Send subsetted and normalized counts to tab delimited file
dds_filt <- dds[aeolidoidea_novel_genes, ]

# Save filtered data set and counts
genes <- as.vector(rownames(counts(dds_filt,normalized=T)))
write.table(cbind(as.data.frame(genes), as.data.frame(counts(dds_filt,normalized=T))), 
            file = paste0(outputPrefix, "-filtered-normalized-counts.txt"), sep = '\t', quote=FALSE)

# Clustering analysis
# Excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
library("RColorBrewer")
library("gplots")
distsRL <- dist(t(assay(rld[aeolidoidea_novel_genes, ])))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds_filt),sampleNames)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
png(paste0(outputPrefix, "-clustering.png"), width=1100, height=1000)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13), cexRow = 1.7, cexCol = 1.7)
dev.off()

# Principal components plot shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(rld[aeolidoidea_novel_genes, ]))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(rld)[select,]))

condition <- treatments
sampleCondition2 <- c("Brain", "Brain",
                      "Distal Ceras","Distal Ceras","Distal Ceras",
                      "Foot","Foot",
                      "Oral Tentacle","Oral Tentacle","Oral Tentacle",
                      "Proximal Ceras","Proximal Ceras","Proximal Ceras",
                      "Rhinophore","Rhinophore","Rhinophore",
                      "Tail","Tail")
scores <- data.frame(pc$x, sampleCondition2)

# PCA 1 vs PCA 2
(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(sampleCondition2))))
  + geom_point(size = 5)
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA),
    text = element_text(size = 20)
  ))

ggsave(pcaplot,file=paste0(outputPrefix, "-pcaplot-1v2.pdf"))
dev.off()

# PCA 2 vs PCA 3
(pcaplot2 <- ggplot(scores, aes(x = PC2, y = PC3, col = (factor(sampleCondition2))))
  + geom_point(size = 5)
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA),
    text = element_text(size = 20)
  ))

ggsave(pcaplot2,file=paste0(outputPrefix, "-pcaplot-2v3.pdf"))
dev.off()

# PCA 1 vs PCA 3
(pcaplot3 <- ggplot(scores, aes(x = PC1, y = PC3, col = (factor(sampleCondition2))))
  + geom_point(size = 5)
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA),
    text = element_text(size = 20)
  ))

ggsave(pcaplot3,file=paste0(outputPrefix, "-pcaplot-1v3.pdf"))
dev.off()

# Data heatmaps
library("RColorBrewer")
library("gplots")

# Clade-specific expression with heatmap.2
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
png(file=paste0(outputPrefix, "-HEATMAP-clade-expression.png"), width=1000, height=1000)
heatmap.2(assay(rld[aeolidoidea_novel_genes, ]), col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="Clade-specific Expression Heatmap")
dev.off()

################################################
### Venn Diagram from presence/absence data ####
################################################
library(ggpolypath)
library(RColorBrewer)
library(venn)

# Create presence/absence data frames for venn diagram construction
dds_filt_df <- as.data.frame(counts(dds_filt,normalized=T))
dds_filt_df.2 <- dds_filt_df[rowSums(dds_filt_df) > 10.0,]

# Create empty venn diagram data frame
venn_df <- data.frame(matrix(ncol = 0, nrow = nrow(dds_filt_df.2)))
rownames(venn_df) <- rownames(dds_filt_df.2)

# Calculate average expression for each tissue type in new data frame
venn_df$DC <- (dds_filt_df.2$`Distal Ceras 1` + dds_filt_df.2$`Distal Ceras 2` + dds_filt_df.2$`Distal Ceras 3`)/3
venn_df$PC <- (dds_filt_df.2$`Proximal Ceras 1` + dds_filt_df.2$`Proximal Ceras 2` + dds_filt_df.2$`Proximal Ceras 3`)/3
venn_df$OT <- (dds_filt_df.2$`Oral Tentacle 1` + dds_filt_df.2$`Oral Tentacle 2` + dds_filt_df.2$`Oral Tentacle 3`)/3
venn_df$RH <- (dds_filt_df.2$`Rhinophore 1` + dds_filt_df.2$`Rhinophore 2` + dds_filt_df.2$`Rhinophore 3`)/3
venn_df$BR <- (dds_filt_df.2$`Brain 7` + dds_filt_df.2$`Brain 8`)/2
venn_df$TL <- (dds_filt_df.2$`Tail 1` + dds_filt_df.2$`Tail 2`)/2
venn_df$FO <- (dds_filt_df.2$`Foot 1` + dds_filt_df.2$`Foot 2`)/2

# Determine presence absence using average expression: Normalized counts > 0.5 considered present
venn_df[venn_df < 0.5] <- as.integer(0)
venn_df[venn_df >= 0.5] <- as.integer(1) 
colnames(venn_df)<-c("Distal Ceras",  "Proximal Ceras",
                     "Oral Tentacle", "Rhinophore", "Brain",
                    "Tail", "Foot")

# Plot venn diagram
cols = brewer.pal(7, "Set2")

png(file=paste0(outputPrefix, "-venn-diagram.png"), width=1000, height=1000)
venn::venn(venn_df, zcolor=cols, ilcs = 2, sncs = 2, plotsize=100)
dev.off()


################################################
##### Gene list from presence/absence data #####
################################################
library(dplyr)
library(reshape2)
library(R.utils)
library(ProtDomSeq)
library(UniprotR)

# Change presence/absence data set column names
colnames(venn_df) <- c("DC",  "PC", "OT", "RH", "BR", "TA", "FO")

# Pull out gene names for genes only expressed in certain tissues and save
dc_only <- rownames(filter(venn_df, DC=="1" & PC=="0" & OT=="0" & RH=="0" & BR=="0" & TA=="0" & FO=="0"))
ot_only <- rownames(filter(venn_df, DC=="0" & PC=="0" & OT=="1" & RH=="0" & BR=="0" & TA=="0" & FO=="0"))
rh_only <- rownames(filter(venn_df, DC=="0" & PC=="0" & OT=="0" & RH=="1" & BR=="0" & TA=="0" & FO=="0"))
br_only <- rownames(filter(venn_df, DC=="0" & PC=="0" & OT=="0" & RH=="0" & BR=="1" & TA=="0" & FO=="0"))
ta_only <- rownames(filter(venn_df, DC=="0" & PC=="0" & OT=="0" & RH=="0" & BR=="0" & TA=="1" & FO=="0"))
fo_only <- rownames(filter(venn_df, DC=="0" & PC=="0" & OT=="0" & RH=="0" & BR=="0" & TA=="0" & FO=="1"))
pc_only <- rownames(filter(venn_df, DC=="0" & PC=="1" & OT=="0" & RH=="0" & BR=="0" & TA=="0" & FO=="0"))

# Write out tables with genes only expressed in (a) certain tissue(s)
write.table(dc_only, paste0("gene_annotations/", outputPrefix, "-DC-only-genelist.txt"), quote=FALSE, row.names = FALSE, sep="\t", col.names = FALSE)
write.table(ot_only, paste0("gene_annotations/", outputPrefix, "-OT-only-genelist.txt"), quote=FALSE, row.names = FALSE, sep="\t", col.names = FALSE) 
write.table(rh_only, paste0("gene_annotations/", outputPrefix, "-RH-only-genelist.txt"), quote=FALSE, row.names = FALSE, sep="\t", col.names = FALSE) 
write.table(br_only, paste0("gene_annotations/", outputPrefix, "-BR-only-genelist.txt"), quote=FALSE, row.names = FALSE, sep="\t", col.names = FALSE) 
write.table(ta_only, paste0("gene_annotations/", outputPrefix, "-TA-only-genelist.txt"), quote=FALSE, row.names = FALSE, sep="\t", col.names = FALSE) 
write.table(fo_only, paste0("gene_annotations/", outputPrefix, "-FO-only-genelist.txt"), quote=FALSE, row.names = FALSE, sep="\t", col.names = FALSE)
write.table(pc_only, paste0("gene_annotations/", outputPrefix, "-PC-only-genelist.txt"), quote=FALSE, row.names = FALSE, sep="\t", col.names = FALSE)

# Pull in interproscan annotation data 
annot <- read.csv("../inputs/berghia_RM_iso_anysupport_ipscan_2021_11.tsv", sep="\t")
annot_edit <- annot[-c(2:3,7:11)]
colnames(annot_edit) <- c("prot_id", "source", "type", "label", "ipr_id", "go_term")

# Pull in interproscan database accessions
annot_db <- read.table("../../inputs/functional_annotation.txt", header=TRUE, sep="\t")

# Pull in blastp hits separately
blastp <- read.csv("../../inputs/Bs_protein-blasthit-ALL-info.txt", sep="\t", 
                     header=TRUE, as.is=TRUE, na.strings=c("","NA"), fill=TRUE)


# Pull out annotation data for tissue-specific Aeolidoidea genes
### Distal Ceras Only ###
## Interpro Scan annotations ##
# Pull out annotations into separate data frame
dc_only_annot <- data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in dc_only) {
  colnames(dc_only_annot) <- colnames(annot_edit)
  t <- paste(i, ".", sep="")
  s <- subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  dc_only_annot <- rbind(dc_only_annot, s)
}

# Sort annotations by number of times present
dc_label_summary <- sort(summary(dc_only_annot$label), decreasing=TRUE)
dc_label_summary <- dc_label_summary[!dc_label_summary %in% 0]

# Write out list of annotations to files
if (nrow(dc_only_annot)>0) {
  write.table(dc_only_annot, paste0("gene_annotations/", outputPrefix, "-DC-only-annotations.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

if (!is.null(nrow(as.table(dc_label_summary)))) {
  write.table(dc_label_summary, paste0("gene_annotations/", outputPrefix, "-DC-label-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

# Sort go terms by number of times present
dc_goterm_summary <- sort(summary(dc_only_annot$go_term), decreasing=TRUE)
dc_goterm_summary <- dc_goterm_summary[!dc_goterm_summary %in% 0]

# Write out list of go terms to file
if (!is.null(nrow(as.table(dc_goterm_summary)))) {
  write.table(dc_goterm_summary, paste0("gene_annotations/", outputPrefix, "-DC-goterms-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
dc_only_acc <- data.frame("protein_id"=character(0), "GO"=character(0), "IPR"=character(0), 
                          "SignalP_EUK"=character(0), "Pfam"=character(0), stringsAsFactors = FALSE)
for (i in dc_only) {
  t <- paste(i, ".", sep="")
  s <- subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s <- data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  dc_only_acc <- rbind(dc_only_acc, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(dc_only_acc)>0) {
  write.table(dc_only_acc, paste0("gene_annotations/", outputPrefix, "-DC-accessions-summary.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

## BLASTP annotations ##
# Pull out annotations into separate data frame
dc_only_blastp <- data.frame(matrix(ncol = ncol(blastp), nrow = 0))
colnames(dc_only_blastp) <- colnames(blastp)
for (i in dc_only) {
  colnames(dc_only_blastp) <- colnames(blastp)
  t <- paste(i, ".", sep="")
  s <- subset(blastp, grepl(t, blastp$prot_id, fixed=TRUE))
  dc_only_blastp <- rbind(dc_only_blastp, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(dc_only_blastp)>0) {
  write.table(dc_only_blastp, paste0("gene_annotations/", outputPrefix, "-DC-blastp-hits.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

### Oral Tentacle Only ###
## Interpro Scan annotations ##
# Pull out annotations into separate data frame
ot_only_annot <- data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in ot_only) {
  colnames(ot_only_annot) <- colnames(annot_edit)
  t <- paste(i, ".", sep="")
  s <- subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  ot_only_annot <- rbind(ot_only_annot, s)
}

# Sort annotations by number of times present
ot_label_summary <- sort(summary(ot_only_annot$label), decreasing=TRUE)
ot_label_summary <- ot_label_summary[!ot_label_summary %in% 0]

# Write out list of annotations to files
if (nrow(ot_only_annot)>0) {
  write.table(ot_only_annot, paste0("gene_annotations/", outputPrefix, "-OT-only-annotations.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

if (!is.null(nrow(as.table(ot_label_summary)))) {
  write.table(ot_label_summary, paste0("gene_annotations/", outputPrefix, "-OT-label-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

# Sort go terms by number of times present
ot_goterm_summary <- sort(summary(ot_only_annot$go_term), decreasing=TRUE)
ot_goterm_summary <- ot_goterm_summary[!ot_goterm_summary %in% 0]

# Write out list of go terms to file
if (!is.null(nrow(as.table(ot_goterm_summary)))) {
  write.table(ot_goterm_summary, paste0("gene_annotations/", outputPrefix, "-OT-goterms-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
ot_only_acc <- data.frame(matrix(ncol = ncol(annot_db), nrow = 0))
for (i in ot_only) {
  colnames(ot_only_acc) <- colnames(annot_db)
  t <- paste(i, ".", sep="")
  s <- subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s <- data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  ot_only_acc <- rbind(ot_only_acc, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(ot_only_acc)>0) {
  write.table(ot_only_acc, paste0("gene_annotations/", outputPrefix, "-OT-accessions-summary.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

## BLASTP annotations ##
# Pull out annotations into separate data frame
ot_only_blastp <- data.frame(matrix(ncol = ncol(blastp), nrow = 0))
colnames(ot_only_blastp) <- colnames(blastp)
for (i in ot_only) {
  colnames(ot_only_blastp) <- colnames(blastp)
  t <- paste(i, ".", sep="")
  s <- subset(blastp, grepl(t, blastp$prot_id, fixed=TRUE))
  ot_only_blastp <- rbind(ot_only_blastp, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(ot_only_blastp)>0) {
  write.table(ot_only_blastp, paste0("gene_annotations/", outputPrefix, "-OT-blastp-hits.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

### Rhinophore Only ###
## Interpro Scan annotations ##
# Pull out annotations into separate data frame
rh_only_annot <- data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in rh_only) {
  colnames(rh_only_annot) <- colnames(annot_edit)
  t <- paste(i, ".", sep="")
  s <- subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  rh_only_annot <- rbind(rh_only_annot, s)
}

# Sort annotations by number of times present
rh_label_summary <- sort(summary(rh_only_annot$label), decreasing=TRUE)
rh_label_summary <- rh_label_summary[!rh_label_summary %in% 0]

# Write out list of annotations to files
if (nrow(rh_only_annot)>0) {
  write.table(rh_only_annot, paste0("gene_annotations/", outputPrefix, "-RH-only-annotations.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

if (!is.null(nrow(as.table(rh_label_summary)))) {
  write.table(rh_label_summary, paste0("gene_annotations/", outputPrefix, "-RH-label-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

# Sort go terms by number of times present
rh_goterm_summary <- sort(summary(rh_only_annot$go_term), decreasing=TRUE)
rh_goterm_summary <- rh_goterm_summary[!rh_goterm_summary %in% 0]

# Write out list of go terms to file
if (!is.null(nrow(as.table(rh_goterm_summary)))) {
  write.table(rh_goterm_summary, paste0("gene_annotations/", outputPrefix, "-RH-goterms-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
rh_only_acc <- data.frame(matrix(ncol = ncol(annot_db), nrow = 0))
for (i in rh_only) {
  colnames(rh_only_acc) <- colnames(annot_db)
  t <- paste(i, ".", sep="")
  s <- subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s <- data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  rh_only_acc <- rbind(rh_only_acc, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(rh_only_acc)>0) {
  write.table(rh_only_acc, paste0("gene_annotations/", outputPrefix, "-RH-accessions-summary.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

## BLASTP annotations ##
# Pull out annotations into separate data frame
rh_only_blastp <- data.frame(matrix(ncol = ncol(blastp), nrow = 0))
colnames(rh_only_blastp) <- colnames(blastp)
for (i in rh_only) {
  colnames(rh_only_blastp) <- colnames(blastp)
  t <- paste(i, ".", sep="")
  s <- subset(blastp, grepl(t, blastp$prot_id, fixed=TRUE))
  rh_only_blastp <- rbind(rh_only_blastp, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(rh_only_blastp)>0) {
  write.table(rh_only_blastp, paste0("gene_annotations/", outputPrefix, "-RH-blastp-hits.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

### Brain Only ###
## Interpro Scan annotations ##
# Pull out annotations into separate data frame
br_only_annot <- data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in br_only) {
  colnames(br_only_annot) <- colnames(annot_edit)
  t <- paste(i, ".", sep="")
  s <- subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  br_only_annot <- rbind(br_only_annot, s)
}

# Sort annotations by number of times present
br_label_summary <- sort(summary(br_only_annot$label), decreasing=TRUE)
br_label_summary <- br_label_summary[!br_label_summary %in% 0]

# Write out list of annotations to files
if (nrow(br_only_annot)>0) {
  write.table(br_only_annot, paste0("gene_annotations/", outputPrefix, "-BR-only-annotations.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

if (!is.null(nrow(as.table(br_label_summary)))){
  write.table(br_label_summary, paste0("gene_annotations/", outputPrefix, "-BR-label-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

# Sort go terms by number of times present
br_goterm_summary <- sort(summary(br_only_annot$go_term), decreasing=TRUE)
br_goterm_summary <- br_goterm_summary[!br_goterm_summary %in% 0]

# Write out list of go terms to file
if (!is.null(nrow(as.table(br_goterm_summary)))) {
  write.table(br_goterm_summary, paste0("gene_annotations/", outputPrefix, "-BR-goterms-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
br_only_acc <- data.frame(matrix(ncol = ncol(annot_db), nrow = 0))
for (i in br_only) {
  colnames(br_only_acc) <- colnames(annot_db)
  t <- paste(i, ".", sep="")
  s <- subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s <- data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  br_only_acc <- rbind(br_only_acc, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(br_only_acc)>0) {
  write.table(br_only_acc, paste0("gene_annotations/", outputPrefix, "-BR-accessions-summary.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

## BLASTP annotations ##
# Pull out annotations into separate data frame
br_only_blastp <- data.frame(matrix(ncol = ncol(blastp), nrow = 0))
colnames(br_only_blastp) <- colnames(blastp)
for (i in br_only) {
  colnames(br_only_blastp) <- colnames(blastp)
  t <- paste(i, ".", sep="")
  s <- subset(blastp, grepl(t, blastp$prot_id, fixed=TRUE))
  br_only_blastp <- rbind(br_only_blastp, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(br_only_blastp)>0) {
  write.table(br_only_blastp, paste0("gene_annotations/", outputPrefix, "-BR-blastp-hits.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

### Foot Only ###
## Interpro Scan annotations ##
# Pull out annotations into separate data frame
fo_only_annot <- data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in fo_only) {
  colnames(fo_only_annot) <- colnames(annot_edit)
  t <- paste(i, ".", sep="")
  s <- subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  fo_only_annot <- rbind(fo_only_annot, s)
}

# Sort annotations by number of times present
fo_label_summary <- sort(summary(fo_only_annot$label), decreasing=TRUE)
fo_label_summary <- fo_label_summary[!fo_label_summary %in% 0]

# Write out list of annotations to files
if (nrow(fo_only_annot)>0) {
  write.table(fo_only_annot, paste0("gene_annotations/", outputPrefix, "-FO-only-annotations.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

if (!is.null(nrow(as.table(fo_label_summary)))) {
  write.table(fo_label_summary, paste0("gene_annotations/", outputPrefix, "-FO-label-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

# Sort go terms by number of times present
fo_goterm_summary <- sort(summary(fo_only_annot$go_term), decreasing=TRUE)
fo_goterm_summary <- fo_goterm_summary[!fo_goterm_summary %in% 0]

# Write out list of go terms to file
if (!is.null(nrow(as.table(fo_goterm_summary)))) {
  write.table(fo_goterm_summary, paste0("gene_annotations/", outputPrefix, "-FO-goterms-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
fo_only_acc <- data.frame(matrix(ncol = ncol(annot_db), nrow = 0))
for (i in fo_only) {
  colnames(fo_only_acc) <- colnames(annot_db)
  t <- paste(i, ".", sep="")
  s <- subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s <- data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  fo_only_acc <- rbind(fo_only_acc, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(fo_only_acc)>0) {
  write.table(fo_only_acc, paste0("gene_annotations/", outputPrefix, "-FO-accessions-summary.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

## BLASTP annotations ##
# Pull out annotations into separate data frame
fo_only_blastp <- data.frame(matrix(ncol = ncol(blastp), nrow = 0))
colnames(fo_only_blastp) <- colnames(blastp)
for (i in fo_only) {
  colnames(fo_only_blastp) <- colnames(blastp)
  t <- paste(i, ".", sep="")
  s <- subset(blastp, grepl(t, blastp$prot_id, fixed=TRUE))
  fo_only_blastp <- rbind(fo_only_blastp, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(fo_only_blastp)>0) {
  write.table(fo_only_blastp, paste0("gene_annotations/", outputPrefix, "-FO-blastp-hits.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

### Tail Only ###
## Interpro Scan annotations ##
# Pull out annotations into separate data frame
ta_only_annot <- data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in ta_only) {
  colnames(ta_only_annot) <- colnames(annot_edit)
  t <- paste(i, ".", sep="")
  s <- subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  ta_only_annot <- rbind(ta_only_annot, s)
}

# Sort annotations by number of times present
ta_label_summary <- sort(summary(ta_only_annot$label), decreasing=TRUE)
ta_label_summary <- ta_label_summary[!ta_label_summary %in% 0]

# Write out list of annotations to files
if (nrow(ta_only_annot)>0) {
  write.table(ta_only_annot, paste0("gene_annotations/", outputPrefix, "-TA-only-annotations.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

if (!is.null(nrow(as.table(ta_label_summary)))) {
  write.table(ta_label_summary, paste0("gene_annotations/", outputPrefix, "-TA-label-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

# Sort go terms by number of times present
ta_goterm_summary <- sort(summary(ta_only_annot$go_term), decreasing=TRUE)
ta_goterm_summary <- ta_goterm_summary[!ta_goterm_summary %in% 0]

# Write out list of go terms to file
if (!is.null(nrow(as.table(ta_goterm_summary)))) {
  write.table(ta_goterm_summary, paste0("gene_annotations/", outputPrefix, "-TA-goterms-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
ta_only_acc <- data.frame(matrix(ncol = ncol(annot_db), nrow = 0))
for (i in ta_only) {
  colnames(ta_only_acc) <- colnames(annot_db)
  t <- paste(i, ".", sep="")
  s <- subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s <- data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  ta_only_acc <- rbind(ta_only_acc, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(ta_only_acc)>0) {
  write.table(ta_only_acc, paste0("gene_annotations/", outputPrefix, "-TA-accessions-summary.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

## BLASTP annotations ##
# Pull out annotations into separate data frame
ta_only_blastp <- data.frame(matrix(ncol = ncol(blastp), nrow = 0))
colnames(ta_only_blastp) <- colnames(blastp)
for (i in ta_only) {
  colnames(ta_only_blastp) <- colnames(blastp)
  t <- paste(i, ".", sep="")
  s <- subset(blastp, grepl(t, blastp$prot_id, fixed=TRUE))
  ta_only_blastp <- rbind(ta_only_blastp, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(ta_only_blastp)>0) {
  write.table(ta_only_blastp, paste0("gene_annotations/", outputPrefix, "-TA-blastp-hits.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

### Proximal Ceras Only ###
## Interpro Scan annotations ##
# Pull out annotations into separate data frame
pc_only_annot <- data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in pc_only) {
  colnames(pc_only_annot) <- colnames(annot_edit)
  t <- paste(i, ".", sep="")
  s <- subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  pc_only_annot <- rbind(pc_only_annot, s)
}

# Sort annotations by number of times present
pc_label_summary <- sort(summary(pc_only_annot$label), decreasing=TRUE)
pc_label_summary <- pc_label_summary[!pc_label_summary %in% 0]

# Write out list of annotations to files
if (nrow(pc_only_annot)>0) {
  write.table(pc_only_annot, paste0("gene_annotations/", outputPrefix, "-PC-only-annotations.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

if (!is.null(nrow(as.table(pc_label_summary)))) {
  write.table(pc_label_summary, paste0("gene_annotations/", outputPrefix, "-PC-label-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

# Sort go terms by number of times present
pc_goterm_summary <- sort(summary(pc_only_annot$go_term), decreasing=TRUE)
pc_goterm_summary <- pc_goterm_summary[!pc_goterm_summary %in% 0]

# Write out list of go terms to file
if (!is.null(nrow(as.table(pc_goterm_summary)))) {
  write.table(pc_goterm_summary, paste0("gene_annotations/", outputPrefix, "-PC-goterms-summary.txt"), quote=FALSE, col.names = FALSE, sep="\t")
}

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
pc_only_acc <- data.frame(matrix(ncol = ncol(annot_db), nrow = 0))
for (i in pc_only) {
  colnames(pc_only_acc) <- colnames(annot_db)
  t <- paste(i, ".", sep="")
  s <- subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s <- data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  pc_only_acc <- rbind(pc_only_acc, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(pc_only_acc)>0) {
  write.table(pc_only_acc, paste0("gene_annotations/", outputPrefix, "-PC-accessions-summary.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}

## BLASTP annotations ##
# Pull out annotations into separate data frame
pc_only_blastp <- data.frame(matrix(ncol = ncol(blastp), nrow = 0))
colnames(pc_only_blastp) <- colnames(blastp)
for (i in pc_only) {
  colnames(pc_only_blastp) <- colnames(blastp)
  t <- paste(i, ".", sep="")
  s <- subset(blastp, grepl(t, blastp$prot_id, fixed=TRUE))
  pc_only_blastp <- rbind(pc_only_blastp, s)
}

# Write out annotation summary for tissue-specific genes to file
if (nrow(pc_only_blastp)>0) {
  write.table(pc_only_blastp, paste0("gene_annotations/", outputPrefix, "-PC-blastp-hits.txt"), quote=FALSE, row.names = FALSE, sep="\t")
}
