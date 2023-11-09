####################################
# counts_and_location.R
# Written by: Jessica A. Goodheart
# Last Updated: 27 October 2023
# Purpose: Takes lists of clade-specific genes generated using kinfin and produced count charts by
# clade and chromosome/scaffold.
####################################

####################################
# Initial setup
####################################

# Load required packages
require(ggplot2)
require(tidyr)
require(plyr)
require(dplyr)
require(stringr)
require(patchwork)

# Set working directory
setwd("[PATH_TO_DIR]")

####################################
# Novel gene counts for each clade
####################################

# Read in table, add column names, and set factor order for plot
counts = read.table(file="novel_gene_nums")
colnames(counts) = c("clade","novel_genes")
counts$Clades = c("Clades","Clades","Clades","Clades","Clades","Clades","Clades")
counts$clade = factor(counts$clade, counts$clade)

# Set levels
counts$clade <- factor(counts$clade,levels = rev(c("Berghia","Aeolidina","Nudibranchia","Heterobranchia","Gastropoda","Mollusca","Other")))

# Save counts table
counts.save <- counts[,1:2]
colnames(counts.save) <- c("Clade", "Number of Genes")
write.table(counts.save, "Berghia-number-genes-per-clade.txt", quote=FALSE, row.names = FALSE, sep="\t")

# Plot figure
cols <- rev(c("#F8766D","#BB9D00","#00B81F","#00C0B8","#00A5FF","#E76BF3","darkgrey"))
c <- ggplot(counts, aes(x = Clades, y = novel_genes, fill = clade)) +
  geom_col() + 
  scale_fill_manual(values=cols) +
  ylab("Number of clade-specific genes") + 
  theme_bw() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.position = "none", text = element_text(size = 40))
c  

# Save figure to file
png(file="bar_chart_counts.png", width=500, height=900)
c
dev.off()

####################################
# Novel gene counts for each chromosome by clade (NO NORMALIZATION)
####################################

# Read in list of chromosomes/scaffolds for each gene ID
# These files were built using with grep and awk on braker_annotations.gff3 (see get_locations.sh)
Be <- scan("clade_files/Berghia_stephanieae_genes.txt.locations", what="", sep="\n")
Ae <- scan("clade_files/Aeolidina_only_genes.txt.locations", what="", sep="\n")
Nu <- scan("clade_files/Nudibranchia_only_genes.txt.locations", what="", sep="\n")
He <- scan("clade_files/Heterobranchia_only_genes.txt.locations", what="", sep="\n")
Ga <- scan("clade_files/Gastropoda_only_genes.txt.locations", what="", sep="\n")
Mo <- scan("clade_files/Mollusca_only_genes.txt.locations", what="", sep="\n")

# Calculate frequencies and read in lists as data.frames
BeF <- as.data.frame(table(Be))
AeF <- as.data.frame(table(Ae))
NuF <- as.data.frame(table(Nu))
HeF <- as.data.frame(table(He))
GaF <- as.data.frame(table(Ga))
MoF <- as.data.frame(table(Mo))

# Set the same column names for each df
colnames(BeF) = colnames(AeF) = colnames(NuF) = colnames(HeF) =colnames(GaF) = colnames(MoF) = c("chromosome","frequency")

# Merge frequency lists and set column names
chrom_freqs <- join_all(list(BeF, AeF, NuF, HeF, GaF, MoF), by="chromosome") 
colnames(chrom_freqs) <- c("Chromosome","Berghia","Aeolidina","Nudibranchia","Heterobranchia","Gastropoda","Mollusca")

# Subset to only chromosomes, order by chromosome and change names to chromosome number
chrom_freqs <- chrom_freqs[1:15,]
chrom_freqs <- chrom_freqs[order(chrom_freqs$Chromosome),] 
chrom_freqs$Chromosome <- c("1","10","11","12","13","14","15","2","3","4","5","6","7","8","9")

# Create tidy df for dataset while setting each factor so the plot elements are ordered correctly
cf_tidy <- gather(chrom_freqs, "clade", "num_genes", c("Berghia","Aeolidina","Nudibranchia","Heterobranchia","Gastropoda","Mollusca"))
cf_tidy$clade <- factor(cf_tidy$clade,levels = rev(c("Berghia","Aeolidina","Nudibranchia","Heterobranchia","Gastropoda","Mollusca")))
cf_tidy$Chromosome <- factor(cf_tidy$Chromosome, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))

# Plot figure
p <- ggplot(cf_tidy, aes(x = Chromosome, y = num_genes, fill = clade)) +
  geom_col() +
  xlab("Chromosome") + 
  ylab("# Clade-specific Genes") + 
  scale_fill_discrete(name="Clade")  +
  theme_bw() +
  theme(text = element_text(size = 30))
p

# Save figure to file
png(file="genome_distribution.png", width=1200, height=700)
p
dev.off()

png(file="genome_distribution_narrow.png", width=800, height=700)
p
dev.off()

####################################
# Novel gene counts for each chromosome by clade (WITH NORMALIZATION)
####################################

# Pull in sequence lengths for each chromosome
seq_lengths <- read.csv("../../figures/2_genome_stats/barplot/seq_lengths_sorted.csv")
seq_lengths <- seq_lengths[1:15,]
seq_lengths$m <- seq_lengths$l/1000000

# Normalize frequency data by chromosome length
chrom_freqs_norm <- as.data.frame(chrom_freqs$Chromosome)
colnames(chrom_freqs_norm) <- "Chromosome"

chrom_freqs_norm$Berghia <- chrom_freqs$Berghia/seq_lengths$m
chrom_freqs_norm$Aeolidina <- chrom_freqs$Aeolidina/seq_lengths$m
chrom_freqs_norm$Nudibranchia <- chrom_freqs$Nudibranchia/seq_lengths$m
chrom_freqs_norm$Heterobranchia <- chrom_freqs$Heterobranchia/seq_lengths$m
chrom_freqs_norm$Gastropoda <- chrom_freqs$Gastropoda/seq_lengths$m
chrom_freqs_norm$Mollusca <- chrom_freqs$Mollusca/seq_lengths$m

# Create tidy df for dataset while setting each factor so the plot elements are ordered correctly
cf_tidy_norm <- gather(chrom_freqs_norm, "clade", "gpmb", c("Berghia","Aeolidina","Nudibranchia","Heterobranchia","Gastropoda","Mollusca"))
cf_tidy_norm$clade <- factor(cf_tidy_norm$clade,levels = rev(c("Berghia","Aeolidina","Nudibranchia","Heterobranchia","Gastropoda","Mollusca")))
cf_tidy_norm$Chromosome <- factor(cf_tidy_norm$Chromosome, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))

# Plot figure
n <- ggplot(cf_tidy_norm, aes(x = Chromosome, y = gpmb, fill = clade)) +
  geom_col() + 
  scale_fill_manual(values=cols[2:7]) +
  xlab("Chromosome") +
  ylab("# Genes per Mb") +
  theme_bw() +  
  theme(text = element_text(size = 30))
n

# Save figure to file
png(file="genome_distribution_normalized.png", width=1200, height=700)
n
dev.off()

png(file="genome_distribution_normalized_narrow.png", width=800, height=700)
n
dev.off()

####################################
# Gene counts for each chromosome by clade for ALL OTHER genes (NO NORMALIZATION)
####################################

# Read in list of chromosomes/scaffolds for each gene ID
# These files were built using with grep and awk on braker_annotations.gff3 (see get_locations.sh)
ao <- scan("clade_files/All_other_genes.txt.locations", what="", sep="\n")

# # Calculate frequencies and read in list as data.frames
aoF <- as.data.frame(table(ao))
aoF <- aoF[1:15,]

# Set the same column names for each df
colnames(aoF) = c("chromosome","frequency")
aoF$chromosome <- c("1","10","11","12","13","14","15","2","3","4","5","6","7","8","9")
aoF$chromosome <- factor(aoF$chromosome, levels=c("1","2","3","4","5","6","7","8",
                                                  "9","10","11","12","13","14","15"))


# Plot figure
a <- ggplot(aoF, aes(x = chromosome, y = frequency)) +
  geom_col() +
  xlab("Chromosome") + 
  ylab("# of Genes") + 
  theme_bw() +
  theme(text = element_text(size = 20))
a

# Save figure to file
png(file="genome_distribution_AllOtherGenes.png", width=1100, height=700)
a
dev.off()

png(file="genome_distribution_AllOtherGenes_narrow.png", width=800, height=700)
a
dev.off()

####################################
# Gene counts for each chromosome by clade for ALL OTHER genes (WITH NORMALIZATION)
####################################

# Normalize frequency data by chromosome length
aoF$gpmb <- aoF$frequency/seq_lengths$m

# Plot figure
b <- ggplot(aoF, aes(x = chromosome, y = gpmb)) +
  geom_col() +
  xlab("Chromosome") + 
  ylab("# Genes per Mb") + 
  theme_bw() +
  theme(text = element_text(size = 20))
b

# Save figure to file
png(file="genome_distribution_normalized_AllOtherGenes.png", width=1100, height=700)
b
dev.off()

png(file="genome_distribution_normalized_AllOtherGenes_narrow.png", width=800, height=700)
b
dev.off()

####################################
# Gene counts for each chromosome by clade for ALL genes (NO NORMALIZATION)
####################################

# Combine clade-specific data with all other genes
aoF.2 <- data.frame(aoF[,1:2],rep("Other",length(aoF)))
colnames(aoF.2) <- c("Chromosome","num_genes","clade")

cf_tidy_all <- rbind(cf_tidy, aoF.2)
colnames(cf_tidy_all) <- c("Chromosome","Clade","num_genes")
cf_tidy_all$Clade <- factor(cf_tidy_all$Clade,
                            levels = rev(c("Berghia","Aeolidina","Nudibranchia","Heterobranchia","Gastropoda","Mollusca","Other")))

# Plot figure
d <- ggplot(cf_tidy_all, aes(x = Chromosome, y = num_genes, fill=Clade)) +
  geom_col() +
  scale_fill_manual(values=cols) +
  xlab("Chromosome") + 
  ylab("# of Genes") + 
  theme_bw() +
  theme(text = element_text(size = 20))
d

# Save figure to file
png(file="genome_distribution_AllGenes.png", width=1100, height=700)
d
dev.off()

png(file="genome_distribution_AllGenes_narrow.png", width=700, height=700)
d
dev.off()

####################################
# Gene counts for each chromosome by clade for ALL genes (WITH NORMALIZATION)
####################################
# Combine clade-specific data with all other genes
aoF.3 <- cbind(aoF[,c(1,3)],rep("Other",nrow(aoF)))
colnames(aoF.3) <- c("Chromosome","gpmb","clade")

cf_tidy_all_norm <- rbind(cf_tidy_norm, aoF.3)
colnames(cf_tidy_all_norm) <- c("Chromosome","Clade","gpmb")
cf_tidy_all_norm$Clade <- factor(cf_tidy_all_norm$Clade,
                            levels = rev(c("Berghia","Aeolidina","Nudibranchia","Heterobranchia","Gastropoda","Mollusca","Other")))


# Plot figure
e <- ggplot(cf_tidy_all_norm, aes(x = Chromosome, y = gpmb, fill=Clade)) +
  geom_col() +
  scale_fill_manual(values=cols) +
  xlab("Chromosome") + 
  ylab("# of Genes per Mb") + 
  theme_bw() +
  theme(text = element_text(size = 20))
e

# Save figure to file
png(file="genome_distribution_AllGenes_normalized.png", width=1100, height=700)
e
dev.off()

png(file="genome_distribution_AllGenes_normalized_narrow.png", width=700, height=700)
e
dev.off()

####################################
# Plot proportions of novelty by chromosome
####################################

# Create tidy df for dataset while setting each factor so the plot elements are ordered correctly
cf_tidy_all <- cf_tidy_all %>% group_by(Chromosome) %>% mutate(prop=prop.table(num_genes))
cf_tidy_all$prop <- with(cf_tidy_all, ave(num_genes, Chromosome, FUN = prop.table))
   
# Plot figure
p <- ggplot(cf_tidy_all, aes(x = Chromosome, y = prop, fill = Clade)) +
  geom_col() + 
  scale_fill_manual(values=cols) +
  ylab("Proportion of genes") + 
  theme_bw() +
  theme(axis.ticks.x=element_blank(), text = element_text(size = 25))
p
 
# Save figures to file                                             
png(file="genome_distribution_proportions_ALL.png", width=1000, height=1000)
p
dev.off()

png(file="genome_distribution_proportions_ALL_narrow.png", width=700, height=700)
p
dev.off()

####################################
# Chi-square tests of frequency
####################################

# Combine clade-specific data with all other genes
cf_tidy_all <- rbind(cf_tidy, aoF.2)
colnames(cf_tidy_all) <- c("Chromosome","Clade","num_genes")
cf_tidy_all$Clade <- factor(cf_tidy_all$Clade,levels = rev(c("Berghia","Aeolidina","Nudibranchia","Heterobranchia","Gastropoda","Mollusca","Other")))

cf_all <- spread(cf_tidy_all, Chromosome, num_genes)

# Spread out tidy data into table and run chi-squared test
chisq_result <- chisq.test(cf_all[,2:16]) 

####################################
# Annotation data
####################################
# Pull in annotation and clade/tissue specificity data
annot.summary <- read.delim("../../../April_2021_data_analysis/annotations_analysis/tissue_diffexp/deseq2/gene-specificity-annotation-summary.txt", 
                            sep = "\t", header=TRUE)

# Pull in clade-specific lists of genes
berghia = scan("../../../April_2021_data_analysis/annotations_analysis/tissue_diffexp/deseq2/berghia/Berghia_stephanieae_genes.txt", what="", sep="\n")
aeolidina = scan("../../../April_2021_data_analysis/annotations_analysis/tissue_diffexp/deseq2/aeolidina/Aeolidina_only_genes.txt", what="", sep="\n")
nudibranchia = scan("../../../April_2021_data_analysis/annotations_analysis/tissue_diffexp/deseq2/nudibranchia/Nudibranchia_only_genes.txt", what="", sep="\n")
heterobranchia = scan("../../../April_2021_data_analysis/annotations_analysis/tissue_diffexp/deseq2/heterobranchia/Heterobranchia_only_genes.txt", what="", sep="\n")
gastropoda = scan("../../../April_2021_data_analysis/annotations_analysis/tissue_diffexp/deseq2/gastropoda/Gastropoda_only_genes.txt", what="", sep="\n")
mollusca = scan("../../../April_2021_data_analysis/annotations_analysis/tissue_diffexp/deseq2/mollusca/Mollusca_only_genes.txt", what="", sep="\n")
other = scan("../../../April_2021_data_analysis/annotations_analysis/tissue_diffexp/deseq2/all_others/All_other_genes.txt", what="", sep="\n")

# List of clades
strings2 <- c("Berghia","Aeolidina","Nudibranchia","Heterobranchia","Gastropoda","Mollusca","Other")

# Pull in clade-specific annotation data and create data frame
annotation.clade = data.frame(matrix(ncol = 3, nrow = 0))
colnames(annotation.clade) = c("gene", "clade", "annotation")
for (s in strings2) {
  sub <- gsub(" ", "", s)
  assign(paste0(sub, ".annot.clade"), annot.summary %>% filter(str_detect(clade, s), !is.na(db_id), 
                                                               !grepl("putative|hypothetical|uncharacterized",Protein.names, ignore.case=TRUE), 
                                                               !is.na(Protein.names)))
  assign(paste0(sub, ".annot.clade.g"), eval(as.name(paste0(sub, ".annot.clade")))[,1])
  assign(paste0(sub, ".annot.clade2"), annot.summary %>% filter(str_detect(clade, s), !is.na(db_id), 
                                                                grepl("putative|hypothetical|uncharacterized",Protein.names, ignore.case=TRUE), 
                                                                !is.na(Protein.names)))
  assign(paste0(sub, ".annot.clade.u"), eval(as.name(paste0(sub, ".annot.clade2")))[,1])
  for (i in as.character(eval(as.name(tolower(paste0(sub)))))) {
    if (i %in% as.character(eval(as.name(paste0(sub, ".annot.clade.g"))))) {
      r <- data.frame("gene"=i, "clade"=s, "annotation"="annotated")
      annotation.clade <- rbind(annotation.clade, r)
    } else if (i %in% as.character(eval(as.name(paste0(sub, ".annot.clade.u"))))) {
      r <- data.frame("gene"=i, "clade"=s, "annotation"="uncharacterized")
      annotation.clade <- rbind(annotation.clade, r)
    } else {
      r <- data.frame("gene"=i, "clade"=s, "annotation"="unannotated")
      annotation.clade <- rbind(annotation.clade, r)
    }
  }
}

# Calculate proportions
annotation.clade.count <- annotation.clade %>% dplyr::count(clade, annotation)
annotation.clade.count2 <- as.data.frame(annotation.clade.count %>% group_by(clade) %>% dplyr::mutate(freq = n / sum(n)))
annotation.clade.count2$annotation <- factor(annotation.clade.count2$annotation, 
                                             levels = c("annotated", "uncharacterized", "unannotated"))
annotation.clade.count2$clade <- factor(annotation.clade.count2$clade,levels = c("Berghia","Aeolidina","Nudibranchia","Heterobranchia","Gastropoda","Mollusca","Other"))
write.table(annotation.clade.count2, paste0("annotation-data-clades.txt"), quote=FALSE, row.names = FALSE, sep="\t")

# Plot figure
ac <- ggplot(annotation.clade.count2, aes(x = clade, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(values=rev(c("chartreuse4", "purple", "gray57"))) +
  labs(x="Tissue", y="Proportion of Genes") + 
  theme_bw() +
  theme(axis.ticks.x=element_blank(), legend.title=element_blank(), text = element_text(size = 25))
ac

png(file=paste0("annotations-clade-specific-genes.png"), width=1000, height=1000)
ac
dev.off()

# Plot full figure
c <- ggplot(counts, aes(x = Clades, y = novel_genes, fill = clade)) +
  geom_col() + 
  scale_fill_manual(values=cols) +
  labs(x="All Genes", y="Number of clade-specific genes") + 
  theme_bw() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 30), legend.position = "none")

d <- ggplot(cf_tidy_all, aes(x = Chromosome, y = num_genes, fill=Clade)) +
  geom_col() +
  scale_fill_manual(values=cols) +
  xlab("Chromosome") + 
  ylab("# of Genes") + 
  theme_bw() +
  theme(text = element_text(size = 30))

ac <- ggplot(annotation.clade.count2, aes(x = clade, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(values=rev(c("chartreuse4", "purple", "gray57"))) +
  labs(x=NULL, y="Proportion of Genes") + 
  theme_bw() +
  theme(axis.ticks.x=element_blank(), text = element_text(size = 30), 
        axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1), legend.title=element_blank())

nested <- ((c + theme(axis.title.x = element_text(margin = margin(t = -200))))|ac)/(d) + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1,2))
nested

png(file=paste0("genes-figure.png"), width = 1100, height = 1200)
nested 
dev.off()

####################################
# All genomes
####################################
# Pull in and parse data
spp_clusters = read.table("TAXON.attribute_metrics.txt", header=TRUE)
spp_clusters.filt = subset(spp_clusters,select=c(2,10,13))
colnames(spp_clusters.filt) = c("Species","Species-restricted","Shared")

# Tidy data
spp_clusters.filt.tidy = spp_clusters.filt %>% gather(Grouping, Proteins, "Species-restricted":"Shared")

# Plot data
g.plot <- ggplot(spp_clusters.filt.tidy, aes(x = Species, y = Proteins))+
  geom_col(aes(fill = Grouping), width = 0.7) + 
  coord_flip()
g.plot

# Bring in Tree
require("phytools")
tree = read.newick("species_tree_jessica_Sept2023.tre")
tree.r = rotate(tree, node=60)
tree.r = rotate(tree.r, node=61)
tree.r = rotate(tree.r, node=65)
tree.r = rotate(tree.r, node=67)
tree.r = rotate(tree.r, node=72)
tree.r = rotate(tree.r, node=73)
tree.r = rotate(tree.r, node=78)
tree.r = rotate(tree.r, node=81)
tree.r = rotate(tree.r, node=83)
tree.r = rotate(tree.r, node=88)
tree.r = rotate(tree.r, node=90)
tree.r = rotate(tree.r, node=92)
tree.r = rotate(tree.r, node=94)
tree.r = rotate(tree.r, node=95)

# Plot tree
phy <- plot.phylo(tree.r, show.tip.label=FALSE)

# Plot tree next to bar charts
is_tip <- tree.r$edge[,2] <= length(tree.r$tip.label)
ordered_tips <- tree.r$edge[is_tip, 2]
spp_clusters.filt.tidy$Species <- factor(spp_clusters.filt.tidy$Species, levels =  tree.r$tip.label[ordered_tips])

t.plot <- ggplot(spp_clusters.filt.tidy, aes(x = Species, y = Proteins))+
  geom_col(aes(fill = Grouping), width = 0.7) +
  scale_fill_manual(values=c("darkgrey","darkcyan")) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 10), axis.text.x = )
t.plot

# Also by percentages
spp_clusters.filt.tidy.p = spp_clusters.filt.tidy %>%
  dplyr::group_by(Species) %>%
  dplyr::mutate(proportion = Proteins/sum(Proteins)) %>% 

p.plot <- ggplot(spp_clusters.filt.tidy.p, aes(x = Species, y = proportion))+
  geom_col(aes(fill = Grouping), width = 0.7) +
  scale_fill_manual(values=c("darkgrey","#F8766D")) +
  coord_flip()
p.plot

# Proportions
spp_clusters.filt.tidy.prot = select(spp_clusters.filt.tidy.p, -(proportion))
spp_clusters.filt.tidy.prot.spread = spread(spp_clusters.filt.tidy.prot, Grouping, Proteins)
spp_clusters.filt.tidy.prot.spread.mat = data.matrix(spp_clusters.filt.tidy.prot.spread)[,2:3]
spp_clusters.filt.tidy.prot.spread.mat[ , c(1,2)] <- spp_clusters.filt.tidy.prot.spread.mat[ , c(2,1)]
colnames(spp_clusters.filt.tidy.prot.spread.mat)[c(1,2)] <- colnames(spp_clusters.filt.tidy.prot.spread.mat)[c(2,1)]
rownames(spp_clusters.filt.tidy.prot.spread.mat) = spp_clusters.filt.tidy.prot.spread$Species

rooted_cladogram = ape::compute.brlen(tree.r)

pdf("species-specific-genes.pdf", width=8, height=8)
plotTree.barplot(rooted_cladogram,spp_clusters.filt.tidy.prot.spread.mat,scal=0.5,
                 args.plotTree=list(fsize=0.7,lwd=2),
                 args.barplot=list(col=c("darkcyan","grey")))
dev.off()


# Proportions
spp_clusters.filt.tidy.prop = select(spp_clusters.filt.tidy.p, -(Proteins))
spp_clusters.filt.tidy.prop.spread = spread(spp_clusters.filt.tidy.prop, Grouping, proportion)
spp_clusters.filt.tidy.prop.spread.mat = data.matrix(spp_clusters.filt.tidy.prop.spread)[,2:3]
spp_clusters.filt.tidy.prop.spread.mat[ , c(1,2)] <- spp_clusters.filt.tidy.prop.spread.mat[ , c(2,1)]
colnames(spp_clusters.filt.tidy.prop.spread.mat)[c(1,2)] <- colnames(spp_clusters.filt.tidy.prop.spread.mat)[c(2,1)]
rownames(spp_clusters.filt.tidy.prop.spread.mat) = spp_clusters.filt.tidy.prop.spread$Species

pdf("species-specific-genes-proportion.pdf", width=8, height=8)
plotTree.barplot(rooted_cladogram,spp_clusters.filt.tidy.prop.spread.mat,scal=0.5,
                 args.plotTree=list(fsize=0.7,lwd=2),
                 args.barplot=list(col=c("darkcyan","grey"),legend.text=TRUE,args.legend=list(x=1,y=-16)))
dev.off()
  
