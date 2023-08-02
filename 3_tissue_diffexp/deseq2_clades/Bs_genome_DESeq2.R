####################################
# Bs_genome_DESeq2.R
# Written by: Jessica A. Goodheart
# Last Updated: 30 May 2023
# Purpose: To analyze differential expression data from Berghia tissues in all genes.
# Inputs used: Expression counts from htseq-count
####################################

####################################
# Initial setup
####################################

# Set the working directory
directory <- "[PATH TO]/tissue_diffexp/deseq2_clades/"
read_dir <- "[PATH TO]/tissue_diffexp/read_files/"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "Bs_tissues_DESeq2"

####################################
# DESeq2 Analyses
# Modified from: https://benchtobioinformatics.wordpress.com/category/dexseq/
####################################

# Load DESeq2 library for differential expression
library("DESeq2")

# Read in counts files for each tissue
sampleFiles<- c("Bb7-BsV1_genome.counts",
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
sampleNames <- c("Brain 7", "Brain 8",
                 "Distal Ceras 1", "Distal Ceras 2", "Distal Ceras 3",
                 "Foot 1", "Foot 2",
                 "Oral Tentacle 1", "Oral Tentacle 2", "Oral Tentacle 3",
                 "Proximal Ceras 1", "Proximal Ceras 2", "Proximal Ceras 3",
                 "Rhinophore 1", "Rhinophore 2", "Rhinophore 3",
                 "Tail 1", "Tail 2")

# Set sample conditions for counts
sampleCondition <- c("brain", "brain",
                     "distCeras","distCeras","distCeras",
                     "foot","foot",
                     "oralTent","oralTent","oralTent",
                     "proxCeras","proxCeras","proxCeras",
                     "rhinophore","rhinophore","rhinophore",
                     "tail","tail")

# Create DESeq2 sample table
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

# Set treatment levels for DESeq2 analysis
treatments = c("brain","distCeras","foot","oralTent","proxCeras","rhinophore","tail")

# Create DESeq2 data set from counts
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = read_dir,
                                       design = ~1 + condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

# Calculate differential expression based on the negative binomial (i.e., Gamma-Poisson) distribution
dds <- DESeq(ddsHTSeq)

write.table(assay(as.data.frame(dds)), 
            file = paste0(outputPrefix, "-deseq2-dds-object.txt"), sep = '\t', quote=FALSE, row.names = FALSE)


# Calculate results
res <- results(dds)

# Stats
res_nocounts <- subset(res,baseMean==0.000000)
num_nocounts <- nrow(res_nocounts)

statsfile<-file("summary_stats.txt")
writeLines(c(print(paste0("Number of genes not expressed: ", num_nocounts)),
             print(paste0("Proportion of genes not expressed: ", (num_nocounts/nrow(res)))),
             print(paste0("Proportion of genes not expressed: ", ((nrow(res)-num_nocounts)/nrow(res))))), statsfile)
close(statsfile)

# Save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-full-normalized-counts.csv"), quote=FALSE)

# order results by padj value (most significant to least)
res <- res[order(res$padj),]

# Produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.csv"), quote=FALSE)

# Replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. Recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates. 

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.05,
             cleaned = results(ddsClean)$padj < 0.05)
addmargins(tab)
write.csv(as.data.frame(tab),file = paste0(outputPrefix, "-replaceoutliers.csv"), quote=FALSE)

resClean <- results(ddsClean)
resClean = subset(res, padj<0.05)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "-replaceoutliers-results.csv"), quote=FALSE)

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
plotMA(res, ylim=c(-8,8))
dev.copy(png, paste0(outputPrefix, "-MAplot_initial_analysis.png"))
dev.off()

# Transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# Plot to show effect of transformation
# Axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,paste0(outputPrefix, "-variance_stabilizing.png"))
dev.off()

# Save normalized values
write.table(cbind(as.data.frame(rownames(assay(rld))), as.data.frame(assay(rld))), 
            file = paste0(outputPrefix, "-rlog-transformed-counts.txt"), sep = '\t', quote=FALSE, row.names = FALSE)
write.table(cbind(as.data.frame(rownames(assay(rld))), as.data.frame(assay(vsd))), 
            file = paste0(outputPrefix, "-vst-transformed-counts.txt"), sep = '\t', quote=FALSE, row.names = FALSE)

save.image("Bs_genome_DESeq2.Rdata")
