####################################
# size_barplot.R
# Written by: Jessica A. Goodheart
# Last Updated: 12 April 2023
# Purpose: To generate a barplot of Berghia scaffold sizes
# Inputs used: Genome fasta file to generate size data
####################################

####################################
# Initial setup
####################################

library("seqinr")
library("ggplot2")

setwd("[PATH_TO]/2_genome_stats/barplot")

# Pull out genome data
genome <- read.fasta(file = "../../../../April_2021_data_analysis/genome_files/data_purged_filtered/Berghia_Apr2021_hirise_purged.filtered.fasta")
g_stats <- as.vector(lapply(genome, function(x) summary(x)$length))

l <- vector()

for (i in 1:18) {
  l <- append(l, g_stats[[i]], i-1)
  }

seq_lengths <- data.frame(names(g_stats), l)
seq_lengths_sorted <-seq_lengths[order(seq_lengths$l, decreasing=TRUE),]

seq_lengths_sorted$Scaffold <- 1:18

write.csv(seq_lengths_sorted, "seq_lengths_sorted.csv", row.names=FALSE)

# Plot data
seq_lengths_sorted <- read.csv("seq_lengths_sorted.csv")
ggplot(seq_lengths_sorted, aes(x=Scaffold, y=l)) + 
  geom_bar(stat = "identity") +
  labs(y="Length (bp)") +
  labs(x="Scaffold Number") +
  theme_classic()

tiff("seq_lengths.tiff", units="in", width=5.5, height=5, res=300)
ggplot(seq_lengths_sorted, aes(x=Scaffold, y=l)) + 
  geom_bar(stat = "identity") +
  labs(y="Length (bp)") +
  labs(x="Scaffold Number") +
  theme_classic()  +
  theme(axis.title.x=element_text(size = 16),
        axis.text.x=element_text(size = 12),
        axis.title.y=element_text(size = 16),
        axis.text.y=element_text(size = 12),
        plot.margin = margin(1,1,1,1, "cm"))
dev.off()

png("seq_lengths.png", width = 550, height = 480)
ggplot(seq_lengths_sorted, aes(x=Scaffold, y=l)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.title.x=element_text(size = 15),
        axis.text.x=element_text(size = 11),
        axis.title.y=element_text(size = 15),
        axis.text.y=element_text(size = 11),
        plot.margin = margin(2,2,2,2, "cm")) +
  labs(x="Scaffold Number", y="Length (bp)")
dev.off()
