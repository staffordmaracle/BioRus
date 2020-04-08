# Read in FASTA file
library(seqinr)
bkpn <- read.fasta("C:/Users/Shannon/Documents/BIOL 432/data/BKPN_metazoa_pullout_align.fa")

# Collapse sequences to strings
bkpn.col <- lapply(bkpn, function(x){
  paste(toupper(as.character(x)),
                        collapse="")})

# Only 2155 of the sequences between the two sites are unique!
bkpn.col.unique <- bkpn.col[!duplicated(bkpn.col)]

# Blast sequences
library(annotate)
res <- list()
for (i in (length(bkpn.col.unique)-14):1) {
  blast.results <- tryCatch(
                   blastSequences(bkpn.col.unique[[i]], 
                                  as='data.frame',
                                  timeout=1000,
                                  hitListSize=10),
                   error=function(e){
                     print(e)
                   }
  )
  res[[names(bkpn.col.unique)[i]]] <- blast.results
  save.image("C:/Users/Shannon/Documents/BIOL 432/data/blast.RData")
}

# Divide up into location
bkpn.df <- data.frame("ID"=names(bkpn.col), 
                      "Seq"=unlist(bkpn.col),
                      stringsAsFactors=F)
id.header <- gsub("(\\w)_\\d_(.*)", "\\1", bkpn.df$ID)
bkpn.df$Reservoir <- gsub("(..)(\\d)(..)(\\d).*", "\\1", id.header)
bkpn.df$Site <- gsub("(..)(\\d)(..)(\\d)", "\\1\\2", id.header)


# Use the "count" function to get a counts table depending on whatever
# variables you are interested in.
library(plyr)
site.counts <- bkpn.df %>% 
  dcast(Seq ~ Site)

set.seed(13)
library(vegan)

bray.dist <- vegdist(site.counts[,-1], method="bray", binary = F)
NMDSdat <- metaMDS(bray.dist, k=2, trymax=100)
PDat <- data.frame(NMDS1=NMDSdat$points[,1],
                   NMDS2=NMDSdat$points[,2],
                   Site=bkpn.df$Site)

library(ggplot2)
qplot(x=NMDS1,y=NMDS2,colour=Site,alpha=I(0.6),data=PDat)+theme_bw()


pca.dat <- prcomp(site.counts[,-1])


# Try BLAST+
cmd <- "blastn -db nr -query C:/Users/Shannon/Documents/BIOL\ 432/data/small.fa -remote -out C:/Users/Shannon/Documents/BIOL\ 432/small.out"
system(cmd)


# Couldn't figure out how to run enough BLASTs fast enough
# So instead I will align the sequences and see how that goes

library(ape)
write.fasta(bkpn.col.unique, names=names(bkpn.col.unique), file.out="C:/Users/Shannon/Documents/BIOL 432/data/bkpn_unique.fa")
bkpn.bin <- read.dna("C:/Users/Shannon/Documents/BIOL 432/data/bkpn_unique.fa", format="fasta")

checkAlignment(bkpn.bin)

dist.mat <- as.matrix(dist.dna(bkpn.bin))

library(reshape2)
library(ggplot2)
pdat <- melt(dist.mat)
ggplot(data=pdat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()

bkpn.tree <- nj(dist.mat)
library(ggtree)
ggtree(bkpn.tree)


muscle.align <- muscle(bkpn.bin, quiet=F)
dist.muscle <- dist.dna(muscle.align)

