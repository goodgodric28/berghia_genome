####################################
# Bs_genome_DESeq2_distceras.R
# Written by: Jessica A. Goodheart
# Last Updated: 3 August 2023
# Purpose: To analyze the clade distribution and expression of distceras genes in Berghia
# Inputs used: Expression counts from htseq-count and Orthofinder/kinfin orthogroup results
####################################

####################################
# Initial setup
####################################

# Set the working directory
directory = "[PATH_TO]/tissue_diffexp/deseq2_tissues/distceras/"
read_dir = "[PATH_TO]/tissue_diffexp/read_files/"
setwd(directory)

# Set the prefix for each output file name
outputPrefix = "Bs_tissues_distceras_DESeq2"

####################################
# DESeq2 Analyses
# Modified from: https://benchtobioinformatics.wordpress.com/category/dexseq/
####################################

# Load DESeq2 library for differential expression
library("DESeq2")

# Read in counts files for each tissue
sampleFiles= c("Bb7-BsV1_genome.counts",
                "Bb8-BsV1_genome.counts",
                "D1-BsV1_genome.counts",
                "D2-BsV1_genome.counts",
                "D4-BsV1_genome.counts",
                "F1-BsV1_genome.counts",
                "F4-BsV1_genome.counts",
                "OT1-BsV1_genome.counts",
                "OT2-BsV1_genome.counts",
                "OT3-BsV1_genome.counts",
                "P1-BsV1_genome.counts",
                "P2-BsV1_genome.counts",
                "P3-BsV1_genome.counts",
                "R1-BsV1_genome.counts",
                "R2-BsV1_genome.counts",
                "R4-BsV1_genome.counts",
                "TL3-BsV1_genome.counts",
                "TL4-BsV1_genome.counts")

# Set sample names for counts
sampleNames = c("Brain 7", "Brain 8",
                 "Distal Ceras 1", "Distal Ceras 2", "Distal Ceras 3",
                 "Foot 1", "Foot 2",
                 "Oral Tentacle 1", "Oral Tentacle 2", "Oral Tentacle 3",
                 "Proximal Ceras 1", "Proximal Ceras 2", "Proximal Ceras 3",
                 "Rhinophore 1", "Rhinophore 2", "Rhinophore 3",
                 "Tail 1", "Tail 2")

# Set sample conditions for counts
sampleCondition = c("brain", "brain",
                     "distCeras","distCeras","distCeras",
                     "foot","foot",
                     "oralTent","oralTent","oralTent",
                     "proxCeras","proxCeras","proxCeras",
                     "rhinophore","rhinophore","rhinophore",
                     "tail","tail")

# Create DESeq2 sample table
sampleTable = data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

# Set treatment levels for DESeq2 analysis
treatments = c("brain","distCeras","foot","oralTent","proxCeras","rhinophore","tail")

# Create DESeq2 data set from counts
ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = read_dir,
                                       design = ~1 + condition)
colData(ddsHTSeq)$condition = factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

# Calculate differential expression based on the negative binomial (i.e., Gamma-Poisson) distribution
dds = DESeq(ddsHTSeq)
res = results(dds, contrast=c(-1/6,1,-1/6,-1/6,-1/6,-1/6,-1/6))

# Save data results and normalized reads to csv
resdata = merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] = 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"), quote=FALSE)

# order results by padj value (most significant to least)
res = res[order(res$padj),]

# Produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.csv"))

# Replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. Recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates. 

ddsClean = replaceOutliersWithTrimmedMean(dds)
ddsClean = DESeq(ddsClean)
tab = table(initial = results(dds)$padj < 0.05,
             cleaned = results(ddsClean)$padj < 0.05)
addmargins(tab)
write.csv(as.data.frame(tab),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean = results(ddsClean)
resClean = subset(res, padj<0.05)
resClean = resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# Subset results to only include upregulated distceras genes 
resClean.up = as.data.frame(subset(resClean, log2FoldChange>2))
counts = as.data.frame(counts(dds,normalized =TRUE))
counts.up = subset(counts, subset= rownames(counts) %in% rownames(resClean.up))
resClean.up.merged = merge(resClean.up, counts, by = 'row.names', sort = FALSE)

write.csv(resClean.up.merged,file = paste0(outputPrefix, "-upregulated-stats-filtered.csv"), quote=FALSE, row.names=FALSE)

# Subset results to only include downregulated distceras genes 
resClean.down = as.data.frame(subset(resClean, log2FoldChange < -2))
counts.down = subset(counts, subset= rownames(counts) %in% rownames(resClean.down))
resClean.down.merged = merge(resClean.down, counts, by = 'row.names', sort = FALSE)

write.csv(resClean.down.merged,file = paste0(outputPrefix, "-downregulated-stats-filtered.csv"), quote=FALSE, row.names=FALSE)

####################################################################################
# Exploratory data analysis of RNAseq data with DESeq2
#
# these next R scripts are for a variety of visualization, QC and other plots to
# get a sense of what the RNAseq data looks like based on DESEq2 analysis
#
# 1) MA plot
# 2) rlog stabilization and variance stabiliazation
# 3) variance stabilization plot
# 4) heatmap of clustering analysis
# 5) PCA plot
####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
plotMA(ddsClean, ylim=c(-8,8))
dev.copy(png, paste0(outputPrefix, "-MAplot_initial_analysis.png"))
dev.off()

# Transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld = rlogTransformation(ddsClean, blind=T)
vsd = varianceStabilizingTransformation(ddsClean, blind=T)

# Save normalized values
write.table(as.data.frame(assay(rld)),file = paste0(outputPrefix, "-rlog-transformed-counts.txt"), sep = '\t')
write.table(as.data.frame(assay(vsd)),file = paste0(outputPrefix, "-vst-transformed-counts.txt"), sep = '\t')

# Top expressed genes overall with heatmap.2
library("RColorBrewer")
library("gplots")
library("genefilter")
library("ggplot2")
library("grDevices")

select = order(rowMeans(counts(ddsClean,normalized=T)),decreasing=T)[1:1000]
my_palette = colorRampPalette(c('blue','white',"red"))(n=1000)

png(file=paste0(outputPrefix, "-HEATMAP-top1000.png"), width=1000, height=1000)
heatmap.2(assay(rld)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="Top 1000 Expressed Genes Heatmap")
dev.off()

# Top expressed genes upregulated
png(file=paste0(outputPrefix, "-HEATMAP-top-expressed-tissue.png"), width=1000, height=1000)
heatmap.2(assay(rld)[rownames(resClean.up),], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="Top Differentially Expressed Genes Heatmap")
dev.off()

##########################################
### Venn diagram of upregulated genes ####
##########################################
library(ggpolypath)
library(RColorBrewer)
library(venn)

# Create presence/absence data frame for venn diagram
venn_df = data.frame(matrix(ncol = 0, nrow = nrow(counts.up)))
rownames(venn_df) = rownames(counts.up)

# Calculate average expression for each tissue type in new data frame
venn_df$DC = (counts.up$`Distal Ceras 1` + counts.up$`Distal Ceras 2` + counts.up$`Distal Ceras 3`)/3
venn_df$PC = (counts.up$`Proximal Ceras 1` + counts.up$`Proximal Ceras 2` + counts.up$`Proximal Ceras 3`)/3
venn_df$OT = (counts.up$`Oral Tentacle 1` + counts.up$`Oral Tentacle 2` + counts.up$`Oral Tentacle 3`)/3
venn_df$RH = (counts.up$`Rhinophore 1` + counts.up$`Rhinophore 2` + counts.up$`Rhinophore 3`)/3
venn_df$BR = (counts.up$`Brain 7` + counts.up$`Brain 8`)/2
venn_df$TL = (counts.up$`Tail 1` + counts.up$`Tail 2`)/2
venn_df$FO = (counts.up$`Foot 1` + counts.up$`Foot 2`)/2

# Determine presence absence using average expression: Normalized counts > 0.5 considered present
venn_df[venn_df > 0.5] = as.integer(1) 
venn_df[venn_df <= 0.5] = as.integer(0) 
colnames(venn_df)=c("Distal Ceras",  "Proximal Ceras",
                     "Oral Tentacle", "Rhinophore", "Brain",
                     "Tail", "Foot")

# Plot venn diagram
cols = brewer.pal(7, "Set2")

png(file=paste0(outputPrefix, "-venn-diagram-up.png"), width=1000, height=1000)
venn::venn(venn_df, zcolor=cols, ilcs = 2, sncs = 2, plotsize=100)
dev.off()

####################################
# Gene annotation analysis - upreg
####################################
library(dplyr)
library(R.utils)
library(ProtDomSeq)

# Pull in interproscan annotation data 
annot = read.csv("../../../genome_annot_nov2021/protein_annotations/filtered/berghia_RM_iso_anysupport_ipscan_2021_11.tsv", sep="\t")
annot_edit = annot[-c(2:3,7:11)]
colnames(annot_edit) = c("prot_id", "source", "type", "label", "ipr_id", "go_term")

# Pull in interproscan database accessions
annot_db = read.table("../../functional_annotation.txt", header=TRUE, sep="\t")

# Pull in blastp hits separately
blastp = read.csv("../../../genome_annot_nov2021/protein_annotations/filtered/Bs_protein-blasthit-ALL-info.txt", sep="\t", 
                   header=TRUE, as.is=TRUE, na.strings=c("","NA"), fill=TRUE)

# Pull out upregulated gene annotations into separate data frame
tissue_annot = data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in rownames(resClean.up)) {
  colnames(tissue_annot) = colnames(annot_edit)
  t = paste(i, ".", sep="")
  s = subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  tissue_annot = rbind(tissue_annot, s)
}

# Sort annotations by number of times present
label_summary = sort(summary(as.factor(tissue_annot$label)), decreasing=TRUE)
label_summary = label_summary[!label_summary %in% 0]

# Write out list of annotations to files
write.table(tissue_annot, paste0(outputPrefix, "-annotations-up.txt"), quote=FALSE, row.names = FALSE)
write.table(label_summary, paste0(outputPrefix, "-label-summary-up.txt"), quote=FALSE, col.names = FALSE)

# Sort go terms by number of times present
goterm_summary = sort(summary(as.factor(tissue_annot$go_term)), decreasing=TRUE)
goterm_summary = goterm_summary[!goterm_summary %in% 0]

# Write out list of go terms to file
write.table(goterm_summary, paste0(outputPrefix, "-goterms-summary-up.txt"), quote=FALSE, col.names = FALSE)

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
acc = data.frame("protein_id"=character(0), "GO"=character(0), "IPR"=character(0), 
                          "SignalP_EUK"=character(0), "Pfam"=character(0), stringsAsFactors = FALSE)
for (i in rownames(resClean.up)) {
  t = paste(i, ".", sep="")
  s = subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s = data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                    "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  acc = rbind(acc, s)
}

# Write out annotation summary for tissue-specific genes to file
write.table(acc, paste0(outputPrefix, "-accessions-summary-up.txt"), quote=FALSE, row.names = FALSE)

## BLASTP annotations ##
# Pull out annotations into separate data frame
tissue_blastp = data.frame(matrix(ncol = ncol(blastp), nrow = 0))
colnames(tissue_blastp) = colnames(blastp)
for (i in rownames(resClean.up)) {
  t = paste(i, ".", sep="")
  s = subset(blastp, grepl(t, blastp$prot_id, fixed=TRUE))
  tissue_blastp = rbind(tissue_blastp, s)
}

# Write out annotation summary for tissue-specific genes to file
write.table(tissue_blastp, paste0(outputPrefix, "-blastp-hits-up.txt"), quote=FALSE, row.names = FALSE, sep="\t")

####################################
# GO term analysis - upregulated
####################################

# Library
require("goseq")
require("dplyr")
require("tidyr")
require("GO.db")

# Annotation data
annot.go = data.frame(gsub(".t[0-9]*", "", annot_db$protein_id), annot_db$GO)
colnames(annot.go) = c("gene_id", "go.terms")
annot.go.edit = as.data.frame(annot.go %>% mutate(go.terms = strsplit(as.character(go.terms), ";")) %>% unnest(go.terms))

# Gene lengths data
gene.lengths.df = read.table("../gene_lengths.txt", sep = " ")
colnames(gene.lengths.df) = c("gene", "length")
gene.lengths = unlist(gene.lengths.df$length)
names(gene.lengths) = gene.lengths.df$gene

# Get GO data for all genes
de.genes = as.integer(names(gene.lengths) %in% rownames(resClean.up))
names(de.genes) = names(gene.lengths)

# Fit probability weight function
pwf = nullp(de.genes, bias.data=gene.lengths)
head(pwf)

# Goseq analysis
GO.wall = goseq(pwf, gene2cat = annot.go.edit)
GO.samp = goseq(pwf, gene2cat = annot.go.edit, method="Sampling", repcnt=1000)
GO.nobias = goseq(pwf, gene2cat = annot.go.edit, method="Hypergeometric")

# Write out significant go.term summary for tissue-specific genes to file
GO.samp.up = subset(GO.samp, over_represented_pvalue<0.05|over_represented_pvalue<0.05)
write.table(GO.samp.up, paste0(outputPrefix, "-goterm-stats-up.txt"), quote=FALSE, row.names = FALSE, sep="\t")

##########################################
### Venn diagram of downregulated genes ##
##########################################
library(ggpolypath)
library(RColorBrewer)
library(venn)

# Create presence/absence data frame for venn diagram
venn_df = data.frame(matrix(ncol = 0, nrow = nrow(counts.down)))
rownames(venn_df) = rownames(counts.down)

# Calculate average expression for each tissue type in new data frame
venn_df$DC = (counts.down$`Distal Ceras 1` + counts.down$`Distal Ceras 2` + counts.down$`Distal Ceras 3`)/3
venn_df$PC = (counts.down$`Proximal Ceras 1` + counts.down$`Proximal Ceras 2` + counts.down$`Proximal Ceras 3`)/3
venn_df$OT = (counts.down$`Oral Tentacle 1` + counts.down$`Oral Tentacle 2` + counts.down$`Oral Tentacle 3`)/3
venn_df$RH = (counts.down$`Rhinophore 1` + counts.down$`Rhinophore 2` + counts.down$`Rhinophore 3`)/3
venn_df$BR = (counts.down$`Brain 7` + counts.down$`Brain 8`)/2
venn_df$TL = (counts.down$`Tail 1` + counts.down$`Tail 2`)/2
venn_df$FO = (counts.down$`Foot 1` + counts.down$`Foot 2`)/2

# Determine presence absence using average expression: Normalized counts > 0.5 considered present
venn_df[venn_df > 0.5] = as.integer(1) 
venn_df[venn_df <= 0.5] = as.integer(0) 
colnames(venn_df)=c("Distal Ceras",  "Proximal Ceras",
                    "Oral Tentacle", "Rhinophore", "Brain",
                    "Tail", "Foot")

# Plot venn diagram
cols = brewer.pal(7, "Set2")

png(file=paste0(outputPrefix, "-venn-diagram-down.png"), width=1000, height=1000)
venn::venn(venn_df, zcolor=cols, ilcs = 2, sncs = 2, plotsize=100)
dev.off()

####################################
# Gene annotation analysis - downreg
####################################
library(dplyr)
library(R.utils)
library(ProtDomSeq)

# Pull in interproscan annotation data 
annot = read.csv("../../../genome_annot_nov2021/protein_annotations/filtered/berghia_RM_iso_anysupport_ipscan_2021_11.tsv", sep="\t")
annot_edit = annot[-c(2:3,7:11)]
colnames(annot_edit) = c("prot_id", "source", "type", "label", "ipr_id", "go_term")

# Pull in interproscan database accessions
annot_db = read.table("../../functional_annotation.txt", header=TRUE, sep="\t")

# Pull in blastp hits separately
blastp = read.csv("../../../genome_annot_nov2021/protein_annotations/filtered/Bs_protein-blasthit-ALL-info.txt", sep="\t", 
                  header=TRUE, as.is=TRUE, na.strings=c("","NA"), fill=TRUE)

# Pull out downregulated gene annotations into separate data frame
tissue_annot = data.frame(matrix(ncol = ncol(annot_edit), nrow = 0))
for (i in rownames(resClean.down)) {
  colnames(tissue_annot) = colnames(annot_edit)
  t = paste(i, ".", sep="")
  s = subset(annot_edit, grepl(t, annot_edit$prot_id, fixed=TRUE))
  tissue_annot = rbind(tissue_annot, s)
}

# Sort annotations by number of times present
label_summary = sort(summary(as.factor(tissue_annot$label)), decreasing=TRUE)
label_summary = label_summary[!label_summary %in% 0]

# Write out list of annotations to files
write.table(tissue_annot, paste0(outputPrefix, "-annotations-down.txt"), quote=FALSE, row.names = FALSE)
write.table(label_summary, paste0(outputPrefix, "-label-summary-down.txt"), quote=FALSE, col.names = FALSE)

# Sort go terms by number of times present
goterm_summary = sort(summary(as.factor(tissue_annot$go_term)), decreasing=TRUE)
goterm_summary = goterm_summary[!goterm_summary %in% 0]

# Write out list of go terms to file
write.table(goterm_summary, paste0(outputPrefix, "-goterms-summary-down.txt"), quote=FALSE, col.names = FALSE)

## Annotation accessions for different interproscan databases ##
# Pull out annotations into separate data frame
acc = data.frame("protein_id"=character(0), "GO"=character(0), "IPR"=character(0), 
                 "SignalP_EUK"=character(0), "Pfam"=character(0), stringsAsFactors = FALSE)
for (i in rownames(resClean.down)) {
  t = paste(i, ".", sep="")
  s = subset(annot_db, grepl(t, annot_db$protein_id, fixed=TRUE))
  if(nrow(s)==0) 
    s = data.frame("protein_id"=i, "GO"="None", "IPR"="None", 
                   "SignalP_EUK"="None", "Pfam"="None", stringsAsFactors = FALSE)
  acc = rbind(acc, s)
}

# Write out annotation summary for tissue-specific genes to file
write.table(acc, paste0(outputPrefix, "-accessions-summary-down.txt"), quote=FALSE, row.names = FALSE)

## BLASTP annotations ##
# Pull out annotations into separate data frame
tissue_blastp = data.frame(matrix(ncol = ncol(blastp), nrow = 0))
colnames(tissue_blastp) = colnames(blastp)
for (i in rownames(resClean.down)) {
  t = paste(i, ".", sep="")
  s = subset(blastp, grepl(t, blastp$prot_id, fixed=TRUE))
  tissue_blastp = rbind(tissue_blastp, s)
}

# Write out annotation summary for tissue-specific genes to file
write.table(tissue_blastp, paste0(outputPrefix, "-blastp-hits-down.txt"), quote=FALSE, row.names = FALSE, sep="\t")

####################################
# GO term analysis - downreg
####################################

# Library
require("goseq")
require("dplyr")
require("tidyr")
require("GO.db")

# Annotation data
annot.go = data.frame(gsub(".t[0-9]*", "", annot_db$protein_id), annot_db$GO)
colnames(annot.go) = c("gene_id", "go.terms")
annot.go.edit = as.data.frame(annot.go %>% mutate(go.terms = strsplit(as.character(go.terms), ";")) %>% unnest(go.terms))

# Gene lengths data
gene.lengths.df = read.table("../gene_lengths.txt", sep = " ")
colnames(gene.lengths.df) = c("gene", "length")
gene.lengths = unlist(gene.lengths.df$length)
names(gene.lengths) = gene.lengths.df$gene

# Get GO data for all genes
de.genes = as.integer(names(gene.lengths) %in% rownames(resClean.down))
names(de.genes) = names(gene.lengths)

# Fit probability weight function
pwf = nullp(de.genes, bias.data=gene.lengths)
head(pwf)

# Goseq analysis
GO.wall = goseq(pwf, gene2cat = annot.go.edit)
GO.samp = goseq(pwf, gene2cat = annot.go.edit, method="Sampling", repcnt=1000)
GO.nobias = goseq(pwf, gene2cat = annot.go.edit, method="Hypergeometric")

# Write out significant go.term summary for tissue-specific genes to file
GO.samp.down = subset(GO.samp, over_represented_pvalue<0.05|over_represented_pvalue<0.05)
write.table(GO.samp.down, paste0(outputPrefix, "-goterm-stats-down.txt"), quote=FALSE, row.names = FALSE, sep="\t")


####################################
# Lineage-specific gene analysis
####################################
library(purrr)
library(dplyr)

# Pull in normalized counts data for all genes (from all_genes diff exp)
diff_exp = read.csv("[PATH_TO]/tissue_diffexp/deseq2/Bs_tissues_DESeq2-full-normalized-counts.csv")
diff_exp$X = NULL

# Pull in Berghia gene ids for lineage-specific genes from kinfin analysis
other_genes = scan("[PATH_TO]/tissue_diffexp/deseq2/all_others/gene_annotations/Bs_tissues_all_others_DESeq2-DC-only-genelist.txt", what="", sep="\n")
mollusca_genes = scan("[PATH_TO]/tissue_diffexp/deseq2/mollusca/gene_annotations/Bs_tissues_mollusca_DESeq2-DC-only-genelist.txt", what="", sep="\n")
gastropoda_genes = scan("[PATH_TO]/tissue_diffexp/deseq2/gastropoda/gene_annotations/Bs_tissues_gastropoda_DESeq2-DC-only-genelist.txt", what="", sep="\n")
nudibranchia_genes = scan("[PATH_TO]/tissue_diffexp/deseq2/nudibranchia/gene_annotations/Bs_tissues_nudibranchia_DESeq2-DC-only-genelist.txt", what="", sep="\n")
aeolidina_genes = scan("[PATH_TO]/tissue_diffexp/deseq2/aeolidina/gene_annotations/Bs_tissues_aeolidina_DESeq2-DC-only-genelist.txt", what="", sep="\n")
berghia_genes = scan("[PATH_TO]/tissue_diffexp/deseq2/berghia/gene_annotations/Bs_tissues_berghia_DESeq2-DC-only-genelist.txt", what="", sep="\n")

# Subset results to only include clade-specific genes
diff_exp.other = diff_exp[diff_exp$gene %in% other_genes, ]
diff_exp.mollusca = diff_exp[diff_exp$gene %in% mollusca_genes, ]
diff_exp.gastropoda = diff_exp[diff_exp$gene %in% gastropoda_genes, ]
diff_exp.nudibranchia = diff_exp[diff_exp$gene %in% nudibranchia_genes, ]
diff_exp.aeolidina = diff_exp[diff_exp$gene %in% aeolidina_genes, ]
diff_exp.berghia = diff_exp[diff_exp$gene %in% berghia_genes, ]

# Data frames with lineage-specific categories for each grouping
diff_exp.other$clade = rep("Other", nrow(diff_exp.other))
diff_exp.mollusca$clade = rep("Mollusca", nrow(diff_exp.mollusca))
diff_exp.gastropoda$clade = rep("Gastropoda", nrow(diff_exp.gastropoda))
diff_exp.nudibranchia$clade = rep("Nudibranchia", nrow(diff_exp.nudibranchia))
diff_exp.aeolidina$clade = rep("Aeolidina", nrow(diff_exp.aeolidina))
diff_exp.berghia$clade = rep("Berghia", nrow(diff_exp.berghia))

diff_exp = do.call("rbind", list(diff_exp.mollusca,diff_exp.gastropoda,diff_exp.nudibranchia,diff_exp.aeolidina,diff_exp.berghia,diff_exp.other))
write.table(diff_exp,file = paste0(outputPrefix, "-lin-spec-expression.txt"), sep = '\t', row.names=FALSE)

# Create tidy data frame with lists of genes
lin.spec.genes = as.data.frame(diff_exp$gene)
lin.spec.genes$clade = diff_exp$clade
colnames(lin.spec.genes) = c("gene","clade")
write.table(as.data.frame(lin.spec.genes),file = paste0(outputPrefix, "-lin-spec-genes.txt"), sep = '\t', quote=FALSE, row.names=FALSE)

# Summary table
lin.spec.genes.sum = data.frame(names(summary(as.factor(lin.spec.genes$clade))),unname(summary(as.factor(lin.spec.genes$clade))),rep("Clade",length(summary(as.factor(lin.spec.genes$clade)))))
colnames(lin.spec.genes.sum) = c("Clade","lin_spec_genes","clade")

levs = c("Berghia","Aeolidina","Nudibranchia","Gastropoda","Mollusca","Other")
for (i in 1:length(levs)) {
  n = levs[i]
  if (n %in% as.character(lin.spec.genes.sum$Clade)) {
    next
  } else {
    s = data.frame(n, 0, "Clade")
    colnames(s) = colnames(lin.spec.genes.sum)
    lin.spec.genes.sum = rbind(lin.spec.genes.sum,s)
  }
}

lin.spec.genes.sum$Clade = factor(lin.spec.genes.sum$Clade,levels = levs)
write.table(as.data.frame(lin.spec.genes.sum),file = paste0(outputPrefix, "-lin-spec-genes-summary.txt"), sep = '\t', quote=FALSE, row.names=FALSE)

# Plot figure
cols = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3", "darkgrey")
dev.off()

c = ggplot(lin.spec.genes.sum, aes(x = clade, y = as.integer(lin_spec_genes), fill = Clade)) +
  geom_col() + 
  ylab("Number of clade-specific genes") + 
  scale_fill_manual(values=cols) +
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40))

png(file=paste0(outputPrefix, "-clade-distribution-lin-spec-genes.png"), width=1000, height=1000)
c 
dev.off()

###################################################
### Venn diagram of specific/upregulated genes ####
###################################################

# Create presence/absence data frame for venn diagram
only.exp = sort(as.character(lin.spec.genes$gene))
upreg = (rownames(resClean.up))
dat = union(as.character(lin.spec.genes$gene), rownames(resClean.up))
out = data.frame(only.exp = as.integer(dat %in% as.character(lin.spec.genes$gene)),
                  sig.upreg = as.integer(dat %in% rownames(resClean.up)))
row.names(out) = dat
colnames(out) = c("Only expressed in distal ceras", "Upregulated in distal ceras")

# Plot venn diagram
cols = brewer.pal(7, "Set2")

png(file=paste0(outputPrefix, "-venn-diagram-diffexp-specific.png"), width=1000, height=1000)
venn::venn(out, zcolor=cols, ilcs = 2, sncs = 2, plotsize=100)
dev.off()

# Pull out genes both upregulated and exclusively expressed in the distceras
colnames(out) <- c("OB", "UB")
both <- rownames(filter(out, c(OB=="1" & UB=="1")))

both.annot = data.frame(matrix(ncol = ncol(blastp), nrow = 0))
colnames(both.annot) = colnames(blastp)
for (i in both) {
  t = paste(i, ".", sep="")
  s = subset(blastp, grepl(t, blastp$prot_id, fixed=TRUE))
  if(nrow(s)==0) 
    next
  both.annot = rbind(both.annot, s)
}

# Write out annotation summary for tissue-specific genes to file
write.table(both.annot, paste0(outputPrefix, "-blastp-only-expressed-and-upregulated.txt"), quote=FALSE, row.names = FALSE, sep="\t")

# Get GO data for all genes
both.genes = as.integer(names(gene.lengths) %in% both)
names(both.genes) = names(gene.lengths)

# Fit probability weight function
pwf = nullp(both.genes, bias.data=gene.lengths)
head(pwf)

# Goseq analysis
GO.wall = goseq(pwf, gene2cat = annot.go.edit)
GO.samp = goseq(pwf, gene2cat = annot.go.edit, method="Sampling", repcnt=1000)
GO.nobias = goseq(pwf, gene2cat = annot.go.edit, method="Hypergeometric")

# Write out significant go.term summary for tissue-specific genes to file
GO.samp.both = subset(GO.samp, over_represented_pvalue<0.05|over_represented_pvalue<0.05)
write.table(GO.samp.both, paste0(outputPrefix, "-goterm-stats-only-expressed-and-upregulated.txt"), quote=FALSE, row.names = FALSE, sep="\t")


