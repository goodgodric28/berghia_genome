####################################
# mollusc_genome_plot.R
# Written by: Jessica A. Goodheart
# Last Updated: 25 May 2023
# Purpose: To generate a plot comparing Berghia to other mollusc genomes
# Inputs used: Mollusca genome data from NCBI
####################################

####################################
# Initial setup
####################################

setwd("[PATH_TO]/2_genome_stats/genome_stats_chart/")

# Library
library(ggplot2)
library(ggrepel)

# Pull in dataset
genomes <- read.csv("Mollusc_genome_stats.csv")
genomes2 <- subset(genomes,genomes$Species!="Crepidula atrasolea")

genomes2$Class <- factor(genomes2$Class,levels = c("Berghia","other Heterobranchia","other Gastropoda","Bivalvia","Cephalopoda","Polyplacophora"))
cols <- c("#F8766D","#C49A00","#7CAE00", "#00C1A3", "#35A2FF", "#A58AFF")

########### Longest Scaffold ##############
# Standard plot
png("Mollusca_stats_chart.png")
ggplot(genomes2, aes(x=BUSCO...genome, y=Longest.scaffold, color=Class)) + 
  geom_point(aes(size=genome_size)) +
  theme_classic() +
  theme(axis.title.x=element_text(size = 15),
        axis.text.x=element_text(size = 10),
        axis.title.y=element_text(size = 15),
        axis.text.y=element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10)) +
  labs(x= "BUSCO Completeness (%)", y= "Longest Scaffold (Mb)", size="Genome Size", colour="Group")
dev.off()

# No labels, with legends
png("Mollusca_stats_chart_nolabelswithlegends.png")
ggplot(genomes2, aes(x=BUSCO...genome, y=Longest.scaffold, color=Class)) + 
  geom_point(aes(size=genome_size)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10),) +
  labs(size="Genome Size", colour="Group")
dev.off()

# No labels, no legends
png("Mollusca_stats_chart_nolabelsnolegends.png")
ggplot(genomes2, aes(x=BUSCO...genome, y=Longest.scaffold, color=Class)) + 
  geom_point(aes(size=genome_size), show.legend = FALSE) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),)
dev.off()

# No labels, no legends
png("Mollusca_stats_chart_labelsnolegends.png", width=1100, height=650)
ggplot(genomes2, aes(x=BUSCO...genome, y=Longest.scaffold, color=Class)) + 
  geom_point(aes(size=genome_size)) +
  labs(x= "BUSCO Completeness (%)", y= "Longest Scaffold (Mb)", size="Genome Size", colour="Group") +
  theme_classic() +
  theme(text = element_text(size = 20),legend.position = c(.14,.7)) +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  scale_size_continuous(range = c(5, 10))
dev.off()


ggplot(genomes2, aes(x=BUSCO...genome, y=Longest.scaffold, color=Class)) +
  geom_point(aes(size=genome_size)) +
  labs(x= "BUSCO Completeness (%)", y= "Longest Scaffold (Mb)", size="Genome Size", colour="Group") +
  theme_classic() +
  theme(text = element_text(size = 20),legend.position = c(.2,.75)) +
  scale_size_continuous(range = c(5, 10))

########### N50 ##############
png("Mollusca_stats_chart_N50.png")
ggplot(genomes2, aes(x=BUSCO...genome, y=N50..Mb., color=Class)) + 
  geom_point(aes(size=genome_size)) +
  theme_classic() +
  theme(axis.title.x=element_text(size = 15),
        axis.text.x=element_text(size = 10),
        axis.title.y=element_text(size = 15),
        axis.text.y=element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10)) +
  labs(x= "BUSCO Completeness (%)", y= "Scaffold N50 (Mb)", size="Genome Size (Mb)", colour="Group") +
  scale_size_continuous(range = c(2, 8))
dev.off()

png("Mollusca_stats_chart_logN50.png", width = 1000, height = 550)
ggplot(genomes2, aes(x=BUSCO...genome, y=Log.scaffold.N50...Mb, color=Class)) + 
  geom_point(aes(size=genome_size), alpha=0.6) +
  scale_color_manual(values=cols) +
  theme_classic() +
  theme(axis.title.x=element_text(size = 25),
        axis.text.x=element_text(size = 18),
        axis.title.y=element_text(size = 25),
        axis.text.y=element_text(size = 18),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 18),
        legend.position = c(.16,.57),
        plot.margin = margin(2,2,2,2, "cm")) +
  labs(x= "BUSCO Completeness (%)", y= "Log(Scaffold N50 (Mb)", size="Genome Size (Mb)", colour="Group") +
  scale_size_continuous(range = c(1, 12)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

tiff("Mollusca_stats_chart_logN50.tiff", units="in", width=10.5, height=5.5, res=300)
ggplot(genomes2, aes(x=BUSCO...genome, y=Log.scaffold.N50...Mb, color=Class)) + 
  geom_point(aes(size=genome_size), alpha=0.6) +
  scale_color_manual(values=cols) +
  theme_classic() +
  theme(axis.title.x=element_text(size = 20),
        axis.text.x=element_text(size = 16),
        axis.title.y=element_text(size = 20),
        axis.text.y=element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.position = c(.14,.56),
        plot.margin = margin(1,1,1,1, "cm")) +
  labs(x= "BUSCO Completeness (%)", y= "Log(Scaffold N50 (Mb)", size="Genome Size (Mb)", colour="Group") +
  scale_size_continuous(range = c(1, 10)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

