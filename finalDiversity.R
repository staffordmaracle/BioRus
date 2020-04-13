# Read in FASTA file
library(seqinr)
library(dplyr)
library(plyr)
library(reshape2)
library(ggplot2)
library(ape)
library(phylotools)
library(ggtree)
library(vegan)
bkpn <- read.fasta("C:\\Users\\yukth\\Documents\\BIOL432\\Final Project\\BKPN_metazoa_pullout_align.fa")

# Collapse sequences to strings
bkpn.col <- lapply(bkpn, function(x){
  paste(toupper(as.character(x)),
        collapse="")})

# Only 2155 of the sequences between the two sites are unique!
bkpn.col.unique <- bkpn.col[!duplicated(bkpn.col)]


# Divide up into location
bkpn.df <- data.frame("ID"=names(bkpn.col.unique), 
                      "Seq"=unlist(bkpn.col.unique),
                      stringsAsFactors=F)
id.header <- gsub("(\\w)_\\d_(.*)", "\\1", bkpn.df$ID)
bkpn.df$Reservoir <- gsub("(..)(\\d)(..)(\\d).*", "\\1", id.header)
bkpn.df$Site <- gsub("(..)(\\d)(..)(\\d)", "\\1\\2", id.header)


# Use the "count" function to get a counts table depending on whatever
# variables you are interested in.

site.counts <- bkpn.df %>% 
  dcast(Seq ~ Site)

# Subsetting count sites dataframe to only include sequences which were found more than 5 times
site.counts$total = rowSums(site.counts[,2:15])
abundant.bkpn = site.counts[site.counts$total>5,1:15]

# Subset the bkpn.df to only include sequences that were found more than 5 times
merged = merge(x = abundant.bkpn,y = bkpn.df, "Seq", by.y = "Seq", all.x = FALSE, all.y = FALSE, sort = FALSE)
abundantMerged = select(merged,ID,Seq)

# Put new dataframe into format that can be used by dat2fasta functions
colnames(abundantMerged)[1] = "seq.name"
colnames(abundantMerged)[2] = "seq.text"

# Write fasta file using full data set
write.fasta(bkpn.col.unique, names=names(bkpn.col.unique), file.out="C:\\Users\\yukth\\Documents\\BIOL432\\Final Project\\bkpn_unique.fa")
bkpn.bin <- read.dna("C:\\Users\\yukth\\Documents\\BIOL432\\Final Project\\bkpn_unique.fa", format="fasta")

# Write fasta file using only abundant sequences
dat2fasta(abundantMerged,"bkpn_abundant.fa")
bkpn.abundant.bin <- read.dna("C:\\Users\\yukth\\Documents\\BIOL432\\Final Project\\bkpn_abundant.fa", format="fasta")

# Make distance matrices of both datasets
dist.mat <- as.matrix(dist.dna(bkpn.bin))
dist.mat2 <- as.matrix(dist.dna(bkpn.abundant.bin))

# Make neighbourjoined tree of both matrices
bkpn.tree <- nj(dist.mat)
bkpn.tree.abundant <- nj(dist.mat2)

# Create phylogenies of both data sets
ggtree(bkpn.tree)
ggtree(bkpn.tree.abundant)

# Calculate shannon diversity index and create bar graph for full dataset
site.counts.transposed = as.data.frame(t(site.counts[,2:15]))
diversity.full = diversity(site.counts.transposed)
graph1.df = as.data.frame(diversity.full)
graph1.df$Sites = row.names(graph1.df)
graph1.df$Colours = c("a","a","a","a","a","a","a","b","b","b","b","b","b","b")
colnames(graph1.df)[1] = "Shannon_Index"
ggplot(graph1.df, aes(Sites,Shannon_Index)) + geom_col(aes(fill = Colours)) +  theme(legend.position = "none") + scale_fill_manual(values = c("#6399F0", "#ED9962"))

# Calculate shannon diversity index and create bar graph for abundant dataset
abundant.counts.transposed = as.data.frame((t(abundant.bkpn[,2:15])))
diversity.abundant = diversity(abundant.counts.transposed)
graph2.df = as.data.frame(diversity.abundant)
graph2.df$Sites = row.names(graph2.df)
graph2.df$Colours = c("a","a","a","a","a","a","a","b","b","b","b","b","b","b")
colnames(graph2.df)[1] = "Shannon_Index"
ggplot(graph2.df, aes(Sites,Shannon_Index)) + geom_col(aes(fill = Colours)) +  theme(legend.position = "none") + scale_fill_manual(values = c("#70c41c", "#f02390"))

