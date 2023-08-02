####################################
# Bs_genome_DESeq2_tissues.R
# Written by: Jessica A. Goodheart
# Last Updated: 30 May 2023
# Purpose: To analyze the clade distribution of upregulated genes in particular tissues in Berghia
# Inputs used: Expression counts from htseq-count and Orthofinder/kinfin orthogroup results
####################################

####################################
# Initial setup
####################################
# Packages
require(ggplot2)
require(tidyr)
require(dplyr)
require(stringr)
require(patchwork)
require(DESeq2)
require(RColorBrewer)
require(venn)
require(GGally)
require(ggbreak)

# Set output prefix
outputPrefix <- "Bs_tissues_DESeq2"

# Set the working directory
directory <- "[PATH TO]/tissue_diffexp/deseq2_tissues"
setwd(directory)


####################################
# Tissue-specific genes
####################################
# Pull in data and create data frame
brain <- read.table("brain/Bs_tissues_brain_DESeq2-lin-spec-genes-summary.txt", header=TRUE)
rhinophore <- read.table("rhinophore/Bs_tissues_rhinophore_DESeq2-lin-spec-genes-summary.txt", header=TRUE)
oraltentacle <- read.table("oraltentacle/Bs_tissues_oraltentacle_DESeq2-lin-spec-genes-summary.txt", header=TRUE)
distalceras <- read.table("distceras/Bs_tissues_distceras_DESeq2-lin-spec-genes-summary.txt", header=TRUE)
proximalceras <- read.table("proximalceras/Bs_tissues_proximalceras_DESeq2-lin-spec-genes-summary.txt", header=TRUE)
foot <- read.table("foot/Bs_tissues_foot_DESeq2-lin-spec-genes-summary.txt", header=TRUE)
tail <- read.table("tail/Bs_tissues_tail_DESeq2-lin-spec-genes-summary.txt", header=TRUE)

brain.2 <- data.frame(brain$Clade,brain$lin_spec_genes,brain$lin_spec_genes/sum(brain$lin_spec_genes),rep("Brain",6))
rhinophore.2 <- data.frame(rhinophore$Clade,rhinophore$lin_spec_genes,rhinophore$lin_spec_genes/sum(rhinophore$lin_spec_genes),rep("Rhinophore",6))
oraltentacle.2 <- data.frame(oraltentacle$Clade,oraltentacle$lin_spec_genes,oraltentacle$lin_spec_genes/sum(oraltentacle$lin_spec_genes),rep("Oral Tentacle",6))
distalceras.2 <- data.frame(distalceras$Clade,distalceras$lin_spec_genes,distalceras$lin_spec_genes/sum(distalceras$lin_spec_genes),rep("Distal Ceras",6))
proximalceras.2 <- data.frame(proximalceras$Clade,proximalceras$lin_spec_genes,proximalceras$lin_spec_genes/sum(proximalceras$lin_spec_genes),rep("Proximal Ceras",6))
foot.2 <- data.frame(foot$Clade,foot$lin_spec_genes,foot$lin_spec_genes/sum(foot$lin_spec_genes),rep("Foot",6))
tail.2 <- data.frame(tail$Clade,tail$lin_spec_genes,tail$lin_spec_genes/sum(tail$lin_spec_genes),rep("Tail",6))

colnames(brain.2) <- c("Clade","Genes","Props","Tissue")
colnames(rhinophore.2) <- c("Clade","Genes","Props","Tissue")
colnames(oraltentacle.2) <- c("Clade","Genes","Props","Tissue")
colnames(distalceras.2) <- c("Clade","Genes","Props","Tissue")
colnames(proximalceras.2) <- c("Clade","Genes","Props","Tissue")
colnames(foot.2) <- c("Clade","Genes","Props","Tissue")
colnames(tail.2) <- c("Clade","Genes","Props","Tissue")

# Tidy
#lineage_data <- rbind(brain.2,rhinophore.2,ot.2,distalceras.2)
#lineage_data$Clade <- factor(lineage_data$Clade,levels = c("all others","Mollusca","Gastropoda","Nudibranchia","Aeolidoidea","Berghia"))

# Tidy 
lineage_data <- rbind(brain.2,rhinophore.2,oraltentacle.2,distalceras.2,proximalceras.2,foot.2,tail.2)
lineage_data$Clade <- factor(lineage_data$Clade,levels = c("Other","Mollusca","Gastropoda","Nudibranchia","Aeolidoidea","Berghia"))

# Plot figure
c <- ggplot(lineage_data, aes(x = Tissue, y = Genes)) +
  geom_col() + 

  ylab("Number of tissue-specific genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.position = "none", text = element_text(size = 25))

png(file=paste0(outputPrefix, "-number-of-tissue-specific-genes.png"), width=1000, height=1000)
c 
dev.off()

##### Percentages ######

# Plot figure
p <- ggplot(lineage_data, aes(x = Tissue, y = Props, fill = Clade)) +
  geom_col() + 
  scale_fill_manual(values=rev(c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3","grey"))) +
  ylab("Proportion of tissue-specific genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), text = element_text(size = 25))

png(file=paste0(outputPrefix, "-clade-distribution-tissue-specific-genes.png"), width=1000, height=1000)
p
dev.off()

####################################
# Upregulated genes for each tissue
####################################
# Pull in annotation and clade/tissue specificity data
annot.summary <- read.delim("../deseq2_clades/gene-specificity-annotation-summary.txt", sep = "\t", header=TRUE)

# Pull in lists of genes upregulated in each tissue 
brain.up <- read.table("brain/Bs_tissues_brain_DESeq2-upregulated-stats-filtered.csv", header=TRUE, sep = ",")
rhinophore.up <- read.table("rhinophore/Bs_tissues_rhinophore_DESeq2-upregulated-stats-filtered.csv", header=TRUE, sep = ",")
oraltentacle.up <- read.table("oraltentacle/Bs_tissues_oraltentacle_DESeq2-upregulated-stats-filtered.csv", header=TRUE, sep = ",")
distalceras.up <- read.table("distceras/Bs_tissues_distceras_DESeq2-upregulated-stats-filtered.csv", header=TRUE, sep =",")
proximalceras.up <- read.table("proximalceras/Bs_tissues_proximalceras_DESeq2-upregulated-stats-filtered.csv", header=TRUE, sep =",")
foot.up <- read.table("foot/Bs_tissues_foot_DESeq2-upregulated-stats-filtered.csv", header=TRUE, sep = ",")
tail.up <- read.table("tail/Bs_tissues_tail_DESeq2-upregulated-stats-filtered.csv", header=TRUE, sep = ",")

brain.up <- brain.up[,1]
rhinophore.up <- rhinophore.up[,1]
oraltentacle.up <- oraltentacle.up[,1]
distalceras.up <- distalceras.up[,1]
proximalceras.up <- proximalceras.up[,1]
foot.up <- foot.up[,1]
tail.up <- tail.up[,1]

# Subset annotation summary to only upreg genes for each tissue
brain.up.clade <- annot.summary[annot.summary$genes %in% brain.up,c(1:2)]
rhinophore.up.clade <- annot.summary[annot.summary$genes %in% rhinophore.up,c(1:2)]
oraltentacle.up.clade <- annot.summary[annot.summary$genes %in% oraltentacle.up,c(1:2)]
distalceras.up.clade <- annot.summary[annot.summary$genes %in% distalceras.up,c(1:2)]
proximalceras.up.clade <- annot.summary[annot.summary$genes %in% proximalceras.up,c(1:2)]
foot.up.clade <- annot.summary[annot.summary$genes %in% foot.up,c(1:2)]
tail.up.clade <- annot.summary[annot.summary$genes %in% tail.up,c(1:2)]

# Summarize number of genes that are clade-specific for each tissue
brain.up.clade.sum <- summary(as.factor(brain.up.clade$clade))
rhinophore.up.clade.sum <- summary(as.factor(rhinophore.up.clade$clade))
oraltentacle.up.clade.sum <- summary(as.factor(oraltentacle.up.clade$clade))
distalceras.up.clade.sum <- summary(as.factor(distalceras.up.clade$clade))
proximalceras.up.clade.sum <- summary(as.factor(proximalceras.up.clade$clade))
foot.up.clade.sum <- summary(as.factor(foot.up.clade$clade))
tail.up.clade.sum <- summary(as.factor(tail.up.clade$clade))

# Create data frame for each tissue
brain.up.df <- data.frame(names(brain.up.clade.sum), brain.up.clade.sum, brain.up.clade.sum/sum(brain.up.clade.sum),rep("Brain",length(brain.up.clade.sum)))
rhinophore.up.df <- data.frame(names(rhinophore.up.clade.sum), rhinophore.up.clade.sum, rhinophore.up.clade.sum/sum(rhinophore.up.clade.sum),rep("Rhinophore",length(rhinophore.up.clade.sum)))
oraltentacle.up.df <- data.frame(names(oraltentacle.up.clade.sum), oraltentacle.up.clade.sum, oraltentacle.up.clade.sum/sum(oraltentacle.up.clade.sum),rep("Oral Tentacle",length(oraltentacle.up.clade.sum)))
distalceras.up.df <- data.frame(names(distalceras.up.clade.sum), distalceras.up.clade.sum, distalceras.up.clade.sum/sum(distalceras.up.clade.sum),rep("Distal Ceras",length(distalceras.up.clade.sum)))
proximalceras.up.df <- data.frame(names(proximalceras.up.clade.sum), proximalceras.up.clade.sum, proximalceras.up.clade.sum/sum(proximalceras.up.clade.sum),rep("Proximal Ceras",length(proximalceras.up.clade.sum)))
foot.up.df <- data.frame(names(foot.up.clade.sum), foot.up.clade.sum, foot.up.clade.sum/sum(foot.up.clade.sum),rep("Foot",length(foot.up.clade.sum)))
tail.up.df <- data.frame(names(tail.up.clade.sum), tail.up.clade.sum, tail.up.clade.sum/sum(tail.up.clade.sum),rep("Tail",length(tail.up.clade.sum)))

colnames(brain.up.df) <- colnames(rhinophore.up.df) <- colnames(oraltentacle.up.df) <- colnames(distalceras.up.df) <- colnames(proximalceras.up.df) <- colnames(foot.up.df) <- colnames(tail.up.df) <- c("Clade","Genes","Props","Tissue")

# Tidy (no brain)
up.lineage_data <- rbind(brain.up.df,rhinophore.up.df,oraltentacle.up.df,distalceras.up.df,proximalceras.up.df,foot.up.df,tail.up.df)
rownames(up.lineage_data) <- 1:nrow(up.lineage_data)
up.lineage_data$Clade <- factor(up.lineage_data$Clade,levels = c("Other","Mollusca","Gastropoda","Nudibranchia","Aeolidoidea","Berghia"))

# Tables of upregulated genes and props
up.genes <- pivot_wider(up.lineage_data, id_cols = Clade, names_from = Tissue, values_from = Genes)
up.genes[is.na(up.genes)] <- 0
up.genes <- data.frame(up.genes %>% bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total")))
up.genes$Total <- rowSums(up.genes[,2:8])

up.props <- data.frame(pivot_wider(up.lineage_data, id_cols = Clade, names_from = Tissue, values_from = Props))

write.table(up.genes, paste0(outputPrefix, "-clade-distribution-upregulated-gene-nums.txt"), quote=FALSE, row.names = FALSE, sep="\t")
write.table(up.props, paste0(outputPrefix, "-clade-distribution-upregulated-gene-props.txt"), quote=FALSE, row.names = FALSE, sep="\t")

# Plot figure
cu <- ggplot(up.lineage_data, aes(x = Tissue, y = Genes)) +
  geom_col() + 
  ylab("Number of upregulated genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.position = "none", text = element_text(size = 25))
cu + 
  scale_y_break(c(300, 500), scales=0.5, ticklabels=c(600)) +
  scale_y_break(c(700, 15000), scales=0.3, ticklabels=c(15500, 16000)) +
  ylim(c(0,16000))

png(file=paste0(outputPrefix, "-upregulated-genes.png"), width=1000, height=1000)
cu 
dev.off()

# Genes only upregulated in one tissue type
# Distribution of genes upregulated across each tissue type
up.list <- c(as.character(brain.up), as.character(rhinophore.up), as.character(oraltentacle.up), as.character(distalceras.up),
             as.character(proximalceras.up), as.character(foot.up), as.character(tail.up))
up.list2 <- unique(up.list)

up.dat <- do.call("data.frame", list(brain = as.integer(up.list2 %in% as.character(brain.up)),
                                     rhinophore = as.integer(up.list2 %in% as.character(rhinophore.up)),
                                     oraltentacle = as.integer(up.list2 %in% as.character(oraltentacle.up)),
                                     distalceras = as.integer(up.list2 %in% as.character(distalceras.up)),
                                     proximalceras = as.integer(up.list2 %in% as.character(proximalceras.up)),
                                     foot = as.integer(up.list2 %in% as.character(foot.up)),
                                     tail = as.integer(up.list2 %in% as.character(tail.up))))
rownames(up.dat) <- up.list2

# Plot venn diagram
cols = brewer.pal(7, "Set2")

png(file=paste0(outputPrefix, "-upregulated-venn-diagram.png"), width=1000, height=1000)
venn::venn(up.dat, zcolor=cols, ilcs = 1, sncs = 1.5, plotsize=100)
dev.off()

# Genes only upregulated in one tissue
# Pull out gene names for genes only expressed in certain tissues and save
nrow(annotation.up[!duplicated(annotation.up$gene),])
colnames(up.dat) <- c("BR",  "RH", "OT", "DC", "PC", "FO", "TA")

dc_only <- rownames(filter(up.dat, DC=="1" & PC=="0" & OT=="0" & RH=="0" & BR=="0" & TA=="0" & FO=="0"))
ot_only <- rownames(filter(up.dat, DC=="0" & PC=="0" & OT=="1" & RH=="0" & BR=="0" & TA=="0" & FO=="0"))
rh_only <- rownames(filter(up.dat, DC=="0" & PC=="0" & OT=="0" & RH=="1" & BR=="0" & TA=="0" & FO=="0"))
br_only <- rownames(filter(up.dat, DC=="0" & PC=="0" & OT=="0" & RH=="0" & BR=="1" & TA=="0" & FO=="0"))
ta_only <- rownames(filter(up.dat, DC=="0" & PC=="0" & OT=="0" & RH=="0" & BR=="0" & TA=="1" & FO=="0"))
fo_only <- rownames(filter(up.dat, DC=="0" & PC=="0" & OT=="0" & RH=="0" & BR=="0" & TA=="0" & FO=="1"))
pc_only <- rownames(filter(up.dat, DC=="0" & PC=="1" & OT=="0" & RH=="0" & BR=="0" & TA=="0" & FO=="0"))

all_only <- c(br_only,rh_only,ot_only,dc_only,pc_only,fo_only,ta_only)

# Subset upregulated dataframes to only those genes upregulated in one tissue
brain.up.clade.only <- brain.up.clade[brain.up.clade$genes %in% br_only,]
rhinophore.up.clade.only <- rhinophore.up.clade[rhinophore.up.clade$genes %in% rh_only,]
oraltentacle.up.clade.only <- oraltentacle.up.clade[oraltentacle.up.clade$genes %in% ot_only,]
distalceras.up.clade.only <- distalceras.up.clade[distalceras.up.clade$genes %in% dc_only,]
proximalceras.clade.up.only <- proximalceras.up.clade[proximalceras.up.clade$genes %in% pc_only,]
foot.up.clade.only <- foot.up.clade[foot.up.clade$genes %in% fo_only,]
tail.up.clade.only <- tail.up.clade[tail.up.clade$genes %in% ta_only,]

# Summarize number of genes that are clade-specific for each tissue
brain.up.clade.sum.only <- summary(as.factor(brain.up.clade.only$clade))
rhinophore.up.clade.sum.only <- summary(as.factor(rhinophore.up.clade.only$clade))
oraltentacle.up.clade.sum.only <- summary(as.factor(oraltentacle.up.clade.only$clade))
distalceras.up.clade.sum.only <- summary(as.factor(distalceras.up.clade.only$clade))
proximalceras.up.clade.sum.only <- summary(as.factor(proximalceras.clade.up.only$clade))
foot.up.clade.sum.only <- summary(as.factor(foot.up.clade.only$clade))
tail.up.clade.sum.only <- summary(as.factor(tail.up.clade.only$clade))

# Create data frame for each tissue
brain.up.only.df <- data.frame(names(brain.up.clade.sum.only), brain.up.clade.sum.only, brain.up.clade.sum.only/sum(brain.up.clade.sum.only),rep("Brain",length(brain.up.clade.sum.only)))
rhinophore.up.only.df <- data.frame(names(rhinophore.up.clade.sum.only), rhinophore.up.clade.sum.only, rhinophore.up.clade.sum.only/sum(rhinophore.up.clade.sum.only),rep("Rhinophore",length(rhinophore.up.clade.sum.only)))
oraltentacle.up.only.df <- data.frame(names(oraltentacle.up.clade.sum.only), oraltentacle.up.clade.sum.only, oraltentacle.up.clade.sum.only/sum(oraltentacle.up.clade.sum.only),rep("Oral Tentacle",length(oraltentacle.up.clade.sum.only)))
distalceras.up.only.df <- data.frame(names(distalceras.up.clade.sum.only), distalceras.up.clade.sum.only, distalceras.up.clade.sum.only/sum(distalceras.up.clade.sum.only),rep("Distal Ceras",length(distalceras.up.clade.sum.only)))
proximalceras.up.only.df <- data.frame(names(proximalceras.up.clade.sum.only), proximalceras.up.clade.sum.only, proximalceras.up.clade.sum.only/sum(proximalceras.up.clade.sum.only),rep("Proximal Ceras",length(proximalceras.up.clade.sum.only)))
foot.up.only.df <- data.frame(names(foot.up.clade.sum.only), foot.up.clade.sum.only, foot.up.clade.sum.only/sum(foot.up.clade.sum.only),rep("Foot",length(foot.up.clade.sum.only)))
tail.up.only.df <- data.frame(names(tail.up.clade.sum.only), tail.up.clade.sum.only, tail.up.clade.sum.only/sum(tail.up.clade.sum.only),rep("Tail",length(tail.up.clade.sum.only)))

colnames(brain.up.only.df) <- colnames(rhinophore.up.only.df) <- colnames(oraltentacle.up.only.df) <- colnames(distalceras.up.only.df) <- colnames(proximalceras.up.only.df) <- colnames(foot.up.only.df) <- colnames(tail.up.only.df) <- c("Clade","Genes","Props","Tissue")

# Tidy (no brain)
up.lineage_data.only <- rbind(brain.up.only.df,rhinophore.up.only.df,oraltentacle.up.only.df,distalceras.up.only.df,proximalceras.up.only.df,foot.up.only.df,tail.up.only.df)
rownames(up.lineage_data.only) <- 1:nrow(up.lineage_data.only)
up.lineage_data.only$Clade <- factor(up.lineage_data.only$Clade,levels = c("Other","Mollusca","Gastropoda","Nudibranchia","Aeolidoidea","Berghia"))

# Tables of upregulated genes and props
up.genes.only <- pivot_wider(up.lineage_data.only, id_cols = Clade, names_from = Tissue, values_from = Genes)
up.genes.only[is.na(up.genes.only)] <- 0
up.genes.only <- data.frame(up.genes.only %>% bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total")))

up.props.only <- data.frame(pivot_wider(up.lineage_data.only, id_cols = Clade, names_from = Tissue, values_from = Props))

write.table(up.genes.only, paste0(outputPrefix, "-clade-distribution-upregulated-onlyone-gene-nums.txt"), quote=FALSE, row.names = FALSE, sep="\t")
write.table(up.props.only, paste0(outputPrefix, "-clade-distribution-upregulated-onlyone-gene-props.txt"), quote=FALSE, row.names = FALSE, sep="\t")

# Plot figure
cou <- ggplot(up.lineage_data.only, aes(x = Tissue, y = Genes)) +
  geom_col() + 
  
  ylab("Number of upregulated genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.position = "none", text = element_text(size = 25))

png(file=paste0(outputPrefix, "-upregulated-onlyone-genes.png"), width=1000, height=1000)
cou 
dev.off()

##### Percentages ######

# Plot figure
pou <- ggplot(up.lineage_data.only, aes(x = Tissue, y = Props, fill = Clade)) +
  geom_col() + 
  scale_fill_manual(values=rev(c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3","grey"))) +
  ylab("Proportion of upregulated genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), text = element_text(size = 25))

png(file=paste0(outputPrefix, "-clade-distribution-upregulated-onlyone-genes.png"), width=1000, height=1000)
pou
dev.off()

####################################
# Downregulated genes for each tissue
####################################
# Pull in lists of genes upregulated in each tissue 
brain.down <- read.table("brain/Bs_tissues_brain_DESeq2-downregulated-stats-filtered.csv", header=TRUE, sep = ",")
rhinophore.down <- read.table("rhinophore/Bs_tissues_rhinophore_DESeq2-downregulated-stats-filtered.csv", header=TRUE, sep = ",")
oraltentacle.down <- read.table("oraltentacle/Bs_tissues_oraltentacle_DESeq2-downregulated-stats-filtered.csv", header=TRUE, sep = ",")
distalceras.down <- read.table("distceras/Bs_tissues_distceras_DESeq2-downregulated-stats-filtered.csv", header=TRUE, sep =",")
proximalceras.down <- read.table("proximalceras/Bs_tissues_proximalceras_DESeq2-downregulated-stats-filtered.csv", header=TRUE, sep =",")
foot.down <- read.table("foot/Bs_tissues_foot_DESeq2-downregulated-stats-filtered.csv", header=TRUE, sep = ",")
tail.down <- read.table("tail/Bs_tissues_tail_DESeq2-downregulated-stats-filtered.csv", header=TRUE, sep = ",")

brain.down <- brain.down[,1]
rhinophore.down <- rhinophore.down[,1]
oraltentacle.down <- oraltentacle.down[,1]
distalceras.down <- distalceras.down[,1]
proximalceras.down <- proximalceras.down[,1]
foot.down <- foot.down[,1]
tail.down <- tail.down[,1]

# Subset annotation summary to only downreg genes for each tissue
brain.down.clade <- annot.summary[annot.summary$genes %in% brain.down,c(1:2)]
rhinophore.down.clade <- annot.summary[annot.summary$genes %in% rhinophore.down,c(1:2)]
oraltentacle.down.clade <- annot.summary[annot.summary$genes %in% oraltentacle.down,c(1:2)]
distalceras.down.clade <- annot.summary[annot.summary$genes %in% distalceras.down,c(1:2)]
proximalceras.down.clade <- annot.summary[annot.summary$genes %in% proximalceras.down,c(1:2)]
foot.down.clade <- annot.summary[annot.summary$genes %in% foot.down,c(1:2)]
tail.down.clade <- annot.summary[annot.summary$genes %in% tail.down,c(1:2)]

# Summarize number of genes that are clade-specific for each tissue
brain.down.clade.sum <- summary(as.factor(brain.down.clade$clade))
rhinophore.down.clade.sum <- summary(as.factor(rhinophore.down.clade$clade))
oraltentacle.down.clade.sum <- summary(as.factor(oraltentacle.down.clade$clade))
distalceras.down.clade.sum <- summary(as.factor(distalceras.down.clade$clade))
proximalceras.down.clade.sum <- summary(as.factor(proximalceras.down.clade$clade))
foot.down.clade.sum <- summary(as.factor(foot.down.clade$clade))
tail.down.clade.sum <- summary(as.factor(tail.down.clade$clade))

# Create data frame for each tissue
brain.down.df <- data.frame(names(brain.down.clade.sum), brain.down.clade.sum, brain.down.clade.sum/sum(brain.down.clade.sum),rep("Brain",length(brain.down.clade.sum)))
rhinophore.down.df <- data.frame(names(rhinophore.down.clade.sum), rhinophore.down.clade.sum, rhinophore.down.clade.sum/sum(rhinophore.down.clade.sum),rep("Rhinophore",length(rhinophore.down.clade.sum)))
oraltentacle.down.df <- data.frame(names(oraltentacle.down.clade.sum), oraltentacle.down.clade.sum, oraltentacle.down.clade.sum/sum(oraltentacle.down.clade.sum),rep("Oral Tentacle",length(oraltentacle.down.clade.sum)))
distalceras.down.df <- data.frame(names(distalceras.down.clade.sum), distalceras.down.clade.sum, distalceras.down.clade.sum/sum(distalceras.down.clade.sum),rep("Distal Ceras",length(distalceras.down.clade.sum)))
proximalceras.down.df <- data.frame(names(proximalceras.down.clade.sum), proximalceras.down.clade.sum, proximalceras.down.clade.sum/sum(proximalceras.down.clade.sum),rep("Proximal Ceras",length(proximalceras.down.clade.sum)))
foot.down.df <- data.frame(names(foot.down.clade.sum), foot.down.clade.sum, foot.down.clade.sum/sum(foot.down.clade.sum),rep("Foot",length(foot.down.clade.sum)))
tail.down.df <- data.frame(names(tail.down.clade.sum), tail.down.clade.sum, tail.down.clade.sum/sum(tail.down.clade.sum),rep("Tail",length(tail.down.clade.sum)))

colnames(brain.down.df) <- colnames(rhinophore.down.df) <- colnames(oraltentacle.down.df) <- colnames(distalceras.down.df) <- colnames(proximalceras.down.df) <- colnames(foot.down.df) <- colnames(tail.down.df) <- c("Clade","Genes","Props","Tissue")

# Tidy (no brain)
down.lineage_data <- rbind(brain.down.df,rhinophore.down.df,oraltentacle.down.df,distalceras.down.df,proximalceras.down.df,foot.down.df,tail.down.df)
rownames(down.lineage_data) <- 1:nrow(down.lineage_data)
down.lineage_data$Clade <- factor(down.lineage_data$Clade,levels = c("Other","Mollusca","Gastropoda","Nudibranchia","Aeolidoidea","Berghia"))

down.genes <- pivot_wider(down.lineage_data, id_cols = Clade, names_from = Tissue, values_from = Genes)
down.genes[is.na(down.genes)] <- 0
down.genes <- data.frame(down.genes %>% bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total")))

down.props <- data.frame(pivot_wider(down.lineage_data, id_cols = Clade, names_from = Tissue, values_from = Props))

write.table(down.genes, paste0(outputPrefix, "-clade-distribution-downregulated-gene-nums.txt"), quote=FALSE, row.names = FALSE, sep="\t")
write.table(down.props, paste0(outputPrefix, "-clade-distribution-downregulated-gene-props.txt"), quote=FALSE, row.names = FALSE, sep="\t")


# Plot figure
cd <- ggplot(down.lineage_data, aes(x = Tissue, y = Genes)) +
  geom_col() + 
  
  ylab("Number of downregulated genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.position = "none", text = element_text(size = 25))

png(file=paste0(outputPrefix, "-downregulated-genes.png"), width=1000, height=1000)
cd
dev.off()

##### Percentages ######

# Plot figure
pd <- ggplot(down.lineage_data, aes(x = Tissue, y = Props, fill = Clade)) +
  geom_col() + 
  scale_fill_manual(values=rev(c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3","grey"))) +
  ylab("Proportion of downregulated genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), text = element_text(size = 25))

png(file=paste0(outputPrefix, "-clade-distribution-downregulated-genes.png"), width=1000, height=1000)
pd
dev.off()

####################################
# Annotation data
####################################
strings <- c("brain", "rhinophore", "oral tentacle", "distal ceras", "proximal ceras", "foot", "tail")

# Pull in upregulated annotation data and create data frame
annotation.up = data.frame(matrix(ncol = 3, nrow = 0))
colnames(annotation.up) = c("gene", "tissue", "annotation")
for (s in strings) {
  sub <- gsub(" ", "", s) 
    
  assign(paste0(sub, ".annot.up"), annot.summary %>% dplyr::filter(genes %in% eval(as.name(paste0(sub, ".up"))), !is.na(db_id), 
                                                            !grepl("putative|hypothetical|uncharacterized",Protein.names, ignore.case=TRUE), 
                                                            !is.na(Protein.names)))
  assign(paste0(sub, ".annot.up.g"), eval(as.name(paste0(sub, ".annot.up")))[,1])
  assign(paste0(sub, ".annot.up2"), annot.summary %>% dplyr::filter(genes %in% eval(as.name(paste0(sub, ".up"))), !is.na(db_id), 
                                                                grepl("putative|hypothetical|uncharacterized",Protein.names, ignore.case=TRUE), 
                                                                !is.na(Protein.names)))
  assign(paste0(sub, ".annot.up.u"), eval(as.name(paste0(sub, ".annot.up2")))[,1])
  
  for (i in eval(as.name(paste0(sub, ".up")))) {
    if (i %in% as.character(eval(as.name(paste0(sub, ".annot.up.g"))))) {
      r <- data.frame("gene"=i, "tissue"=s, "annotation"="annotated")
      annotation.up <- rbind(annotation.up, r)
    } else if (i %in% as.character(eval(as.name(paste0(sub, ".annot.up.u"))))) {
      r <- data.frame("gene"=i, "tissue"=s, "annotation"="uncharacterized")
      annotation.up <- rbind(annotation.up, r)
    } else {
      r <- data.frame("gene"=i, "tissue"=s, "annotation"="unannotated")
      annotation.up <- rbind(annotation.up, r)
    }
  }
}

# Calculate proportions
annotation.up.count <- annotation.up %>% dplyr::count(tissue, annotation)
annotation.up.count2 <- as.data.frame(annotation.up.count %>% group_by(tissue) %>% mutate(freq = n / sum(n)))
annotation.up.count2$annotation <- factor(annotation.up.count2$annotation, 
                                             levels = c("annotated", "uncharacterized", "unannotated"))
write.table(annotation.up.count2, paste0(outputPrefix, "-annotation-data-upregulated.txt"), quote=FALSE, row.names = FALSE, sep="\t")

# Plot figure
au <- ggplot(annotation.up.count2, aes(x = tissue, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(values=c("gray57",  "purple", "chartreuse4")) +
  labs(x="Tissue", y="Proportion of upregulated genes") + 
  scale_x_discrete(labels=c("Brain", "Rhinophore", "Oral Tentacle", "Distal Ceras", "Proximal Ceras", "Foot", "Tail")) +
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.title=element_blank(), text = element_text(size = 25))

png(file=paste0(outputPrefix, "-annotations-upregulated-genes.png"), width=1000, height=1000)
au
dev.off()

# Annotations for genes only upregulated in one tissue type
annotation.up.only <- annotation.up[annotation.up$gene %in% all_only,]

# Calculate proportions
annotation.up.only.count <- annotation.up.only %>% dplyr::count(tissue, annotation)
annotation.up.only.count2 <- as.data.frame(annotation.up.only.count %>% group_by(tissue) %>% mutate(freq = n / sum(n)))
annotation.up.only.count2$annotation <- factor(annotation.up.only.count2$annotation, 
                                          levels = c("annotated", "uncharacterized", "unannotated"))
write.table(annotation.up.only.count2, paste0(outputPrefix, "-annotation-data-upregulated-onlyone.txt"), quote=FALSE, row.names = FALSE, sep="\t")

# Plot figure
bu <- ggplot(annotation.up.only.count2, aes(x = tissue, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(values=c("gray57",  "purple", "chartreuse4")) +
  labs(x="Tissue", y="Proportion of upregulated genes") + 
  scale_x_discrete(labels=c("Brain", "Rhinophore", "Oral Tentacle", "Distal Ceras", "Proximal Ceras", "Foot", "Tail")) +
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.title=element_blank(), text = element_text(size = 25))

png(file=paste0(outputPrefix, "-annotations-upregulated-genes-onlyone.png"), width=1000, height=1000)
bu
dev.off()

# Pull in downregulated annotation data and create data frame
strings <- c("brain", "rhinophore", "oral tentacle", "distal ceras", "proximal ceras", "foot", "tail")
annotation.down = data.frame(matrix(ncol = 4, nrow = 0))
colnames(annotation.down) = c("gene", "tissue", "clade", "annotation")
for (s in strings) {
  sub <- gsub(" ", "", s) 
  
  assign(paste0(sub, ".annot.down"), annot.summary %>% dplyr::filter(genes %in% eval(as.name(paste0(sub, ".down"))), !is.na(db_id), 
                                                                   !grepl("putative|hypothetical|uncharacterized",Protein.names, ignore.case=TRUE), 
                                                                   !is.na(Protein.names)))
  assign(paste0(sub, ".annot.down.g"), eval(as.name(paste0(sub, ".annot.down")))[,1])
  assign(paste0(sub, ".annot.down2"), annot.summary %>% dplyr::filter(genes %in% eval(as.name(paste0(sub, ".down"))), !is.na(db_id), 
                                                                    grepl("putative|hypothetical|uncharacterized",Protein.names, ignore.case=TRUE), 
                                                                    !is.na(Protein.names)))
  assign(paste0(sub, ".annot.down.u"), eval(as.name(paste0(sub, ".annot.down2")))[,1])
  
  for (i in eval(as.name(paste0(sub, ".down")))) {
    if (i %in% as.character(eval(as.name(paste0(sub, ".annot.down.g"))))) {
      r <- data.frame("gene"=i, "tissue"=s, "annotation"="annotated")
      annotation.down <- rbind(annotation.down, r)
    } else if (i %in% as.character(eval(as.name(paste0(sub, ".annot.down.u"))))) {
      r <- data.frame("gene"=i, "tissue"=s, "annotation"="uncharacterized")
      annotation.down <- rbind(annotation.down, r)
    } else {
      r <- data.frame("gene"=i, "tissue"=s, "annotation"="unannotated")
      annotation.down <- rbind(annotation.down, r)
    }
  }
}

# Distribution of genes downregulated across each tissue type
down.list <- c(brain.down, rhinophore.down, oraltentacle.down, distalceras.down,
             proximalceras.down, foot.down, tail.down)
down.list <- unique(down.list)
down.dat <- do.call("data.frame", list(brain = as.integer(down.list %in% brain.down),
                                     rhinophore = as.integer(down.list %in% rhinophore.down),
                                     oraltentacle = as.integer(down.list %in% oraltentacle.down),
                                     distalceras = as.integer(down.list %in% distalceras.down),
                                     proximalceras = as.integer(down.list %in% proximalceras.down),
                                     foot = as.integer(down.list %in% foot.down),
                                     tail = as.integer(down.list %in% tail.down)))
rownames(down.dat) <- down.list

# Plot venn diagram
cols = brewer.pal(7, "Set2")

png(file=paste0(outputPrefix, "-downregulated-venn-diagram.png"), width=1000, height=1000)
venn::venn(down.dat, zcolor=cols, ilcs = 2, sncs = 2, plotsize=100)
dev.off()

# Calculate proportions
annotation.down.count <- annotation.down %>% dplyr::count(tissue, annotation)
annotation.down.count2 <- as.data.frame(annotation.down.count %>% group_by(tissue) %>% mutate(freq = n / sum(n)))
annotation.down.count2$annotation <- factor(annotation.down.count2$annotation, 
                                             levels = c("annotated", "uncharacterized", "unannotated"))
write.table(annotation.down.count2, paste0(outputPrefix, "-annotation-data-downregulated.txt"), quote=FALSE, row.names = FALSE, sep="\t")

# Plot figure
ad <- ggplot(annotation.down.count2, aes(x = tissue, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(values=c("gray57",  "purple", "chartreuse4")) +
  labs(x="Tissue", y="Proportion of downregulated genes") + 
  scale_x_discrete(labels=c("Brain", "Rhinophore", "Oral Tentacle", "Distal Ceras", "Proximal Ceras", "Foot", "Tail")) +
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.title=element_blank(), text = element_text(size = 25))

png(file=paste0(outputPrefix, "-annotations-downregulated-genes.png"), width=1000, height=1000)
ad
dev.off()


####################################
# Create full figure
####################################
# Upregulated
cu <- ggplot(up.lineage_data, aes(x = Tissue, y = Genes, fill = Clade)) +
  geom_col() + 
  scale_fill_manual(values=rev(c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3","grey"))) +
  ylab("Num. upregulated genes") + 
  scale_x_discrete(labels=c("BR", "RH", "OT", "DC", "PC", "FO", "TA")) +
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.position = "none", text = element_text(size = 20))
cu2 <- cu + 
  scale_y_break(c(300, 500), scales=0.5, ticklabels=c(500, 1000, 1500, 2000)) +
  scale_y_break(c(2000, 14900), scales=0.3, ticklabels=c(15000, 15200), space=0.3) +
  theme(axis.text.x.top = element_blank(), axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank())

pu <- ggplot(up.lineage_data, aes(x = Tissue, y = Props, fill = Clade)) +
  geom_col() + 
  scale_fill_manual(values=rev(c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3","grey"))) +
  ylab("Prop. upregulated genes") + 
  scale_x_discrete(labels=c("BR", "RH", "OT", "DC", "PC", "FO", "TA")) +
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.position="none", text = element_text(size = 20))

au <- ggplot(annotation.up.count2, aes(x = tissue, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(name="Annotation", values=c("gray57",  "purple", "chartreuse4")) +
  labs(x="Tissue", y="Prop. upregulated genes") + 
  scale_x_discrete(labels=c("BR", "RH", "OT", "DC", "PC", "FO", "TA")) +
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.position="none", text = element_text(size = 20))

colorBlind7  <- c("#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
upreg.venn <- venn::venn(up.dat, zcolor=colorBlind7, ilcs = 0.8, sncs = 1, 
                         plotsize=1000, ggplot=T, box=F) 
nested <- upreg.venn + pu + au +
  plot_annotation(tag_levels = 'A') +
  plot_layout(ncol=3, nrow=1, widths=c(2,1,1), heights=c(1,1,1), guides="collect") &
  theme(legend.text = element_text(size = 15), legend.position="bottom", 
        legend.spacing = unit(2, "cm")) &
  guides(fill=guide_legend(nrow=2))
nested 
dev.off()

ggsave(file=paste0(outputPrefix, "-upregulated-genes-figure.pdf"), 
       plot=nested, width = 15, height = 6)
dev.off()

png(file=paste0(outputPrefix, "-upregulated-genes-figure.png"), width = 1500, height = 650)
nested
dev.off()

# Upregulated in only one tissue type
cou <- ggplot(up.lineage_data.only, aes(x = Tissue, y = Genes, fill = Clade)) +
  geom_col() +
  scale_fill_manual(values=rev(c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3","grey"))) +
  ylab("Num. only upreg. genes") + 
  scale_x_discrete(labels=c("BR", "RH", "OT", "DC", "PC", "FO", "TA")) +
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.position = "none", text = element_text(size = 35))
cou2 <- cou + 
  scale_y_break(c(100, 200), scales=0.3, ticklabels=c(200, 1000, 2000), space=0.3) +
  scale_y_break(c(2000, 14000), scales=0.3, ticklabels=c(14200, 14500), space=0.3) +
  theme(axis.text.y.right = element_blank(), axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.x.top = element_blank(), axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank())

pou <- ggplot(up.lineage_data.only, aes(x = Tissue, y = Props, fill = Clade)) +
  geom_col() + 
  scale_fill_manual(values=rev(c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3","grey"))) +
  ylab("Prop. only upreg. genes") + 
  scale_x_discrete(labels=c("BR", "RH", "OT", "DC", "PC", "FO", "TA")) +
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.position="none", text = element_text(size = 20))

aou <- ggplot(annotation.up.only.count2, aes(x = tissue, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(name="Annotation", values=c("gray57",  "purple", "chartreuse4")) +
  labs(x="Tissue", y="Prop. only upreg. genes") + 
  scale_x_discrete(labels=c("BR", "RH", "OT", "DC", "PC", "FO", "TA")) +
  theme_classic() +
  theme(axis.ticks.x=element_blank(), legend.position="none", text = element_text(size = 20))

upreg.venn <- venn::venn(up.dat, zcolor=cols, ilcs = 2, sncs = 2, plotsize=500, ggplot=TRUE)

nested2 <- pou + aou 
nested2

ggsave(file=paste0(outputPrefix, "-upreg-onlyone-genes-figure.pdf"), 
       plot=nested2, width = 10, height = 6)
dev.off()

#####################################
# Annotations based on clades/tissues
#####################################

# Genes with clade-specificity and annotation inf0
clades.data <- annot.summary[,1:2]
clades.data.sub <- clades.data[match(annotation.up$gene, clades.data$genes),]

# Add clade data to annotation.up
annotation.up$clade <- clades.data.sub$clade

# Subset to only genes upregulated in non-brain tissues
annotation.up.table.nobrain <- annotation.up.table[annotation.up.table$tissue != 'brain',]
annotation.up.table.sub <- annotation.up.table[annotation.up.table$gene %in% annotation.up.table.nobrain$gene,]

# Organize data into factors
annotation.up.table.sub$clade <- factor(annotation.up.table.sub$clade,levels = c("Other","Mollusca","Gastropoda","Nudibranchia","Aeolidoidea","Berghia"))
annotation.up.table.sub$tissue <- factor(annotation.up.table.sub$tissue,levels = c("brain","rhinophore","oral tentacle",
                                                                                   "distal ceras","proximal ceras","foot","tail"))
annotation.up.table.sub$annotation <- factor(annotation.up.table.sub$annotation,levels = c("annotated","uncharacterized","unannotated"))
annotation.up.table.sub$annot <- rep("annot", times=nrow(annotation.up.table.sub))
annotation.up.table.sub <- annotation.up.table.sub[order(annotation.up.table.sub$clade, annotation.up.table.sub$annotation),]
annotation.up.table.sub$gene <- factor(annotation.up.table.sub$gene, levels=unique(annotation.up.table.sub$gene))


f <- ggplot(annotation.up.table.sub, aes(gene, tissue, fill = clade)) +
  geom_tile() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_fill_manual(values=rev(c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3","grey"))) +
  facet_grid2(.~clade, scales = "free_x", space = "free_x", switch = "x", 
              remove_labels = 'x', 
              strip = strip_themed(
                background_x = elem_list_rect(
                  fill = rev(c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3","grey"))))) +
  theme(axis.text.x = element_blank(), axis.line.x = element_blank(),
        axis.title.x = element_blank(), axis.line.x.bottom = element_blank(),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"))
g <- ggplot(annotation.up.table.sub, aes(gene, annot, fill = annotation)) +
  geom_tile() +
  scale_fill_manual(values=c("gray57",  "purple", "chartreuse4"))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme(axis.text.x = element_blank(), axis.line.x = element_blank(),
        axis.title.x = element_blank(), axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(), plot.margin = margin(0, 0, 0, 0, "pt"))
g/f +
  plot_layout(widths=c(1,1), heights = c(1,7))
dev.off()

# Calculate proportions, clade-specificity
annotation.up.count.clades <- annotation.up %>% dplyr::count(tissue, annotation, clade, .drop=FALSE)
annotation.up.count2.clades <- as.data.frame(annotation.up.count.clades %>% group_by(clade, tissue) %>% mutate(freq = n / sum(n)))
annotation.up.count2.clades$annotation <- factor(annotation.up.count2.clades$annotation, 
                                                 levels = c("annotated", "uncharacterized", "unannotated"))
annotation.up.count2.clades$clade <- factor(annotation.up.count2.clades$clade,levels = c("Other","Mollusca","Gastropoda","Nudibranchia","Aeolidoidea","Berghia"))
write.table(annotation.up.count2.clades, paste0(outputPrefix, "-annotation-data-upregulated-clades.txt"), quote=FALSE, row.names = FALSE, sep="\t")

# Table for each tissue
br.up.clade.annot <- subset(annotation.up.count2.clades, tissue == "brain")
rh.up.clade.annot <- subset(annotation.up.count2.clades, tissue == "rhinophore")
ot.up.clade.annot <- subset(annotation.up.count2.clades, tissue == "oral tentacle")
dc.up.clade.annot <- subset(annotation.up.count2.clades, tissue == "distal ceras")
pc.up.clade.annot <- subset(annotation.up.count2.clades, tissue == "proximal ceras")
fo.up.clade.annot <- subset(annotation.up.count2.clades, tissue == "foot")
ta.up.clade.annot <- subset(annotation.up.count2.clades, tissue == "tail")


# plots
au.br <- ggplot(br.up.clade.annot, aes(x = clade, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(name="Annotation", values=c("gray57",  "purple", "chartreuse4")) +
  labs(x="Tissue", y="Proportion of upregulated genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), text = element_text(size = 35), 
        axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1), legend.title=element_blank())

au.rh <- ggplot(rh.up.clade.annot, aes(x = clade, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(name="Annotation", values=c("gray57",  "purple", "chartreuse4")) +
  labs(x="Tissue", y="Proportion of upregulated genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), text = element_text(size = 35), 
        axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1), legend.title=element_blank())

au.ot <- ggplot(ot.up.clade.annot, aes(x = clade, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(name="Annotation", values=c("gray57",  "purple", "chartreuse4")) +
  labs(x="Tissue", y="Proportion of upregulated genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), text = element_text(size = 35), 
        axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1), legend.title=element_blank())

au.dc <- ggplot(dc.up.clade.annot, aes(x = clade, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(name="Annotation", values=c("gray57",  "purple", "chartreuse4")) +
  labs(x="Tissue", y="Proportion of upregulated genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), text = element_text(size = 35), 
        axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1), legend.title=element_blank())

au.pc <- ggplot(pc.up.clade.annot, aes(x = clade, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(name="Annotation", values=c("gray57",  "purple", "chartreuse4")) +
  labs(x="Tissue", y="Proportion of upregulated genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), text = element_text(size = 35), 
        axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1), legend.title=element_blank())

au.fo <- ggplot(fo.up.clade.annot, aes(x = clade, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(name="Annotation", values=c("gray57",  "purple", "chartreuse4")) +
  labs(x="Tissue", y="Proportion of upregulated genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), text = element_text(size = 35), 
        axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1), legend.title=element_blank())

au.ta <- ggplot(ta.up.clade.annot, aes(x = clade, y = freq, fill = annotation)) +
  geom_col() + 
  scale_fill_manual(name="Annotation", values=c("gray57",  "purple", "chartreuse4")) +
  labs(x="Tissue", y="Proportion of upregulated genes") + 
  theme_classic() +
  theme(axis.ticks.x=element_blank(), text = element_text(size = 35), 
        axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1), legend.title=element_blank())

# Combine plots
nested <- ((au.br|au.rh|au.ot)/(au.dc|au.pc|au.fo)) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides="collect") &
  theme(legend.text = element_text(size = 25), legend.position="bottom", 
        legend.spacing = unit(0.8, "cm")) &
  guides(fill=guide_legend(nrow=2))
nested 


####################################
# Summary Table
####################################

up.totals <- up.genes[7,]
down.totals <- down.genes[7,]
