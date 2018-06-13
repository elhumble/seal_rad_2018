# Generate circos plot for full seal : dog genome alignment
# For LAST alignment, see: alignment_notes

library(data.table)
library(dplyr)
options(scipen = 999)
library(RColorBrewer)
library(car)
source("scripts/RCircos.Long.Genome.Link.Plot.R")
library(RCircos)


#~~~~~~~~~~~~
# Read in coordinate files parsed from last output (split by species, order retained)

dog <- fread("data/raw/CanFam_coords.txt",
             col.names = c("Contig", "Start", "AlignLength", "Strand", "Length"))

seal <- fread("data/raw/PBJelly_v1.4_coords.txt", 
              col.names = c("Contig", "Start", "AlignLength", "Strand", "Length"))

summary(dog$AlignLength)
summary(seal$AlignLength)

# 0.1 - 52.8 kb
# mean 2.1 kb

#~~~~~~~~~~~~
# Read in seal scaffold length stats
# Get cumulative coordinate positions as though one chromsome
# Correct alignment positions

seal_lengths <- read.table("data/raw/PBJelly_v1.4_lengths", header = F) %>%  
  `colnames<-`(c("Contig", "Length")) %>%
  mutate(Length = as.numeric(Length)) %>%
  mutate(cum_end = cumsum(Length)) %>%
  mutate(cum_start = cum_end) %>%
  mutate(cum_start=lag(cum_start)) %>%
  mutate(cum_start = replace(cum_start, is.na(cum_start), 0))

seal_align <- seal_lengths %>%
  right_join(seal, by = "Contig") %>%
  .[c(1,5,6,4,3)] %>%
  `colnames<-`(c("Contig", "Start", "AlignLength", "ContigStartCoord", "ContigEndCoord")) %>%
  mutate(NewAlignStart = Start + (ContigStartCoord)) %>% # was contigstartcoord - 1
  mutate(NewAlignEnd = NewAlignStart + AlignLength)

#~~~~~~~~~~~~
# Get dog chromsome length stats

dog_align <- read.table("data/raw/CanFam3.1_lengths.txt", header = F) %>%
  `colnames<-`(c("Contig", "Length")) %>%
  mutate(Contig = gsub(">", "", Contig)) %>%
  mutate(Length = as.numeric(Length)) %>%
  right_join(dog, by = "Contig") %>%
  mutate(End = Start + AlignLength) %>%
  .[c(1,3,7)]


#~~ Coverage of dog

library(IRanges)
library(GenomicRanges) # bioconductor

data <- read.table("data/raw/CanFam3.1_lengths.txt", header = F) %>%
  `colnames<-`(c("Contig", "Length")) %>%
  mutate(Contig = gsub(">", "", Contig)) %>%
  mutate(Length = as.numeric(Length)) %>%
  right_join(dog, by = "Contig") %>%
  mutate(End = Start + AlignLength) %>%
  filter(grepl("^NC", Contig))
  
gr <- GRanges(seqnames = data$Contig, strand = "+",
              ranges = IRanges(start = data$Start, width = data$AlignLength))

# merge overlapping ranges
gr <- reduce(gr)
data <- data.frame(Contig = seqnames(gr),
                   Start = start(gr)-1,
                   End = end(gr)) %>%
  mutate(width = width(gr))

sum(data$width)

dog_lengths <- read.table("data/raw/CanFam3.1_lengths.txt", header = F) %>%
  mutate(V1 = gsub(">", "", V1)) %>%
  filter(grepl("^NC", V1)) %>%
  mutate(V2 = as.numeric(V2))

sum(dog_lengths$V2)

sum(data$width) / sum(dog_lengths$V2) * 100

#~~~~~~~~~~~~
# Explore alignment lengths
# Join both dataframes
# Filter by alignment length

hist(seal_align$AlignLength)
seal_filt <- filter(seal_align, AlignLength >= 20000)
hist(seal_filt$AlignLength)

alignlength <- 5000

align <- cbind(seal_align, dog_align) %>%
  .[c(1,3,6,7,8,9,10)] %>%
  `colnames<-`(c("Chromosome", "AlignLength", "chromStart", "chromEnd", 
                 "Chromosome.1", "chromStart.1", "chromEnd.1")) %>%
  mutate(Chromosome = "chr1") %>%
  filter(grepl("^NC", Chromosome.1)) %>%
  filter(AlignLength >= alignlength) %>%
  .[c(1,3,4,5,6,7)]

# Error checking
link.lengths <- align$chromEnd.1 - align$chromStart.1
which(link.lengths==0)

# Recode dog chromosome names
align$Chromosome.1 <- recode(align$Chromosome.1, '"NC_006583.3" = "chr1"; "NC_006584.3" = "chr2"; 
                             "NC_006585.3" = "chr3"; "NC_006586.3" = "chr4"; 
                             "NC_006587.3" = "chr5"; "NC_006588.3" = "chr6"; 
                             "NC_006589.3" = "chr7"; "NC_006590.3" = "chr8"; 
                             "NC_006591.3" = "chr9"; "NC_006592.3" = "chr10"; 
                             "NC_006593.3" = "chr11"; "NC_006594.3" = "chr12"; 
                             "NC_006595.3" = "chr13"; "NC_006596.3" = "chr14"; 
                             "NC_006597.3" = "chr15"; "NC_006598.3" = "chr16"; 
                             "NC_006599.3" = "chr17"; "NC_006600.3" = "chr18"; 
                             "NC_006601.3" = "chr19"; "NC_006602.3" = "chr20"; 
                             "NC_006603.3" = "chr21"; "NC_006604.3" = "chr22"; 
                             "NC_006605.3" = "chr23"; "NC_006606.3" = "chr24"; 
                             "NC_006607.3" = "chr25"; "NC_006608.3" = "chr26";
                             "NC_006609.3" = "chr27"; "NC_006610.3" = "chr28"; 
                             "NC_006611.3" = "chr29"; "NC_006612.3" = "chr30"; 
                             "NC_006613.3" = "chr31"; "NC_006614.3" = "chr32"; 
                             "NC_006615.3" = "chr33"; "NC_006616.3" = "chr34"; 
                             "NC_006617.3" = "chr35"; "NC_006618.3" = "chr36"; 
                             "NC_006619.3" = "chr37"; "NC_006620.3" = "chr38"; 
                             "NC_006621.3" = "chrX"; "NC_002008.4" = "chrMT"')


# Get link colours
color <- colorRampPalette(brewer.pal(11,"Paired"))(40)
# color <- sample(color) # get a good sample
pie(rep(1, length(color)), col = color , main="") 

# Or read in good color swatch
# color <- read.table("colors_1")
# color <- my_color$x

# Generate link.data
linkcolours <- matrix(nrow = 40,
                      ncol = 2)
linkcolours[,1] <- color
linkcolours[,2] <- c(paste0("chr", seq(1:38)), "chrX", "chrMT")
colnames(linkcolours) <- c("PlotColor", "Chromosome.1")
linkcolours <- as.data.frame(linkcolours)

link.data <- left_join(align, linkcolours, by = "Chromosome.1")

#~~~~~~~~~~~~
# Make Ideograms

# SEAL
data(UCSC.Baylor.3.4.Rat.cytoBandIdeogram);
seal_ideogram <- UCSC.Baylor.3.4.Rat.cytoBandIdeogram[1,]
seal_ideogram[1,3] <- seal_lengths[c(length(seal_lengths[,1])),3]

# DOG
dog_ideogram <- read.table("data/raw/CanFam3.1_lengths.txt", header = F) %>%
  `colnames<-`(c("Chromosome", "ChromEnd")) %>%
  mutate(Chromosome = gsub(">", "", Chromosome)) %>%
  mutate(ChromStart = 0) %>%
  .[c(1,3,2)] %>%
  filter(grepl("^NC", Chromosome)) %>%
  mutate(Band = "qA1", Stain = rep(c("gneg", "gpos"),20))

dog_chr <- c(paste0("chr", seq(1:38)), "chrX", "chrMT")
dog_ideogram$Chromosome <- dog_chr

# Correct species names
species.list <- c("S", "D")
cyto.list <- list(seal_ideogram, 
                  dog_ideogram)
RCircos.Multiple.Species.Core.Components(cyto.list, 
                                         species.list, NULL, 1, 0)

link.data[,1] <- paste(species.list[1], link.data[,1], sep="")
link.data[,4] <- paste(species.list[2], link.data[,4], sep="")

#~~~~~~~~~~~~
# Make Circos Plot

# Set plotting parameters
params <- RCircos.Get.Plot.Parameters();
params$base.per.unit <- 300000;
params$highlight.width <- 0.7;
params$line.color <- "gray";
RCircos.Reset.Plot.Parameters(params);

# Initialize graphic device (GUI or image file)
pdf(file="figs/RCircos.Seal.and.Dog.5000.pdf", height=8, width=8)

png(file=paste0("figs/circos_seal_dog_full_",alignlength,".png"), units = "in", res = 300, height=8, width=8)
par(cex=0.75)
RCircos.Set.Plot.Area()
# Plot chromosome ideogram and link lines
RCircos.Chromosome.Ideogram.Plot()
RCircos.Long.Genome.Link.Plot(link.data, track.num = 1)
dev.off()



#~~~~~~~~~~~~
# Large Blocks
# From maffilter merge

# blocks in dog
# blocks in seal

# maffilter input.file=PBJelly_v1.4_CanFam3.1.maf input.file.compression=none output.log=PBJelly_v1.4_CanFam3.1.maffilter.merge.log maf.filter=“Merge(species=(CanFam), dist_max=10000),Output(file=merge.maf, compression=none)”
