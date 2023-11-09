####################################
# mollusc_genome_plot.R
# Written by: Jessica A. Goodheart
# Last Updated: 31 October 2023
# Purpose: To generate a plot comparing Berghia to other mollusc genomes
# Inputs used: Mollusca genome data from NCBI
####################################

####################################
# Initial setup
####################################

setwd("[PATH_TO_DIR]")

# Library
library(ggplot2)
library(ggrepel)

# Pull in dataset
genomes <- read.table("mollusca_genomes_data.txt",sep = "\t",header = TRUE)

genomes$Group <- factor(genomes$Group,levels = c("Berghia","other Nudibranchia","other Heterobranchia","other Gastropoda","Scaphopoda","Bivalvia","Cephalopoda","Polyplacophora","Aplacophora"))
cols <- c("#F8766D","#E7861B","#BB9D00","#5BB300", "#00C1A3", "#00ABFD", "#9590FF","#DC71FA","#FF62BC") 

########### Longest Scaffold ##############
# Standard plot
png("Mollusca_stats_chart.png")
ggplot(genomes, aes(x=BUSCO.completeness, y=Longest.scaffold, color=Group)) + 
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


ggplot(genomes2, aes(x=BUSCO...genome, y=Longest.scaffold, color=Group)) +
  geom_point(aes(size=genome_size)) +
  labs(x= "BUSCO Completeness (%)", y= "Longest Scaffold (Mb)", size="Genome Size", colour="Group") +
  theme_classic() +
  theme(text = element_text(size = 20),legend.position = c(.2,.75)) +
  scale_size_continuous(range = c(5, 10))

########### N50 ##############
png("Mollusca_stats_chart_N50.png")
ggplot(genomes, aes(x=BUSCO.completeness, y=N50, color=Group)) + 
  geom_point(aes(size=Assembled.Genome.Size)) +
  theme_classic() +
  theme(axis.title.x=element_text(size = 15),
        axis.text.x=element_text(size = 10),
        axis.title.y=element_text(size = 15),
        axis.text.y=element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10)) +
  labs(x= "BUSCO Completeness", y= "Scaffold N50 (Mb)", size="Genome Size (Mb)", colour="Group") +
  scale_size_continuous(range = c(2, 8))
dev.off()

png("Mollusca_stats_chart_logN50.png", width = 1000, height = 550)
ggplot(genomes, aes(x=BUSCO.completeness, y=log.N50., color=Group)) + 
  geom_point(aes(size=Assembled.Genome.Size), alpha=0.6) +
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
  labs(x= "BUSCO Completeness", y= "Log(Scaffold N50 (Mb)", size="Genome Size (Mb)", colour="Group") +
  scale_size_continuous(range = c(1, 12)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

tiff("Mollusca_stats_chart_logN50.tiff", units="in", width=11, height=5.5, res=300)
ggplot(genomes, aes(x=BUSCO.completeness, y=log.N50., color=Group)) + 
  geom_point(aes(size=Assembled.Genome.Size), alpha=0.6) + 
  scale_color_manual(values=cols) +
  theme_classic() +
  theme(axis.title.x=element_text(size = 20),
        axis.text.x=element_text(size = 16),
        axis.title.y=element_text(size = 20),
        axis.text.y=element_text(size = 16),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        #legend.position = c(.2,.6),
        plot.margin = margin(1,1,1,1, "cm")) +
  labs(x= "BUSCO Completeness (%)", y= "Log(Scaffold N50)", size="Genome Size (bases)", colour="Group") +
  scale_size_continuous(range = c(1, 10)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

# Labeled
tiff("Mollusca_stats_chart_logN50_labeled.tiff", units="in", width=33, height=15.5, res=300)
ggplot(genomes, aes(x=BUSCO.completeness, y=log.N50., color=Group)) + 
  geom_point(aes(size=Assembled.Genome.Size), alpha=0.6) + 
  scale_color_manual(values=cols) +
  theme_classic() +
  theme(axis.title.x=element_text(size = 20),
        axis.text.x=element_text(size = 16),
        axis.title.y=element_text(size = 20),
        axis.text.y=element_text(size = 16),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        #legend.position = c(.2,.6),
        plot.margin = margin(1,1,1,1, "cm")) +
  labs(x= "BUSCO Completeness (%)", y= "Log(Scaffold N50 (Mb)", size="Genome Size (Mb)", colour="Group") +
  scale_size_continuous(range = c(1, 10)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  geom_label_repel(aes(label = Species),
                   box.padding   = 0.1, 
                   point.padding = 0.5,
                   segment.color = 'grey50')
dev.off()
