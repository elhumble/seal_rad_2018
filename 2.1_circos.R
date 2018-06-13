# Generate circos plot for seal : dog genome alignment
# Longest 40 seal contigs only
# For LAST alignment, see: alignment_notes

library(data.table)
library(dplyr)
options(scipen = 999)
library(RColorBrewer)
library(car)
source("scripts/RCircos.Long.Genome.Link.Plot.R")
library(RCircos)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Longest 40 scaffolds         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~
# Read in coordinate files parsed from last output (split by species, order retained)

dog <- fread("data/raw/CanFam_coords.txt",
             col.names = c("Contig", "Start", "AlignLength", "Strand", "Length"))

seal <- fread("data/raw/PBJelly_v1.4_coords.txt", 
              col.names = c("Contig", "Start", "AlignLength", "Strand", "Length")) %>%
  mutate(Contig = gsub("_quiver_pilon", "", Contig))

#~~~~~~~~~~~~
# Get dog chromsome length stats

seal_align <- read.table("data/raw/PBJelly_v1.4_lengths", header = F) %>%  
  mutate(V1 = gsub("_quiver_pilon", "", V1)) %>%
  `colnames<-`(c("Contig", "Length")) %>%
  mutate(Length = as.numeric(Length)) %>%
  right_join(seal, by = "Contig") %>%
  mutate(End = Start + AlignLength) 
#.[c(1,3,7)]


dog_align <- read.table("data/raw/CanFam3.1_lengths.txt", header = F) %>%
  `colnames<-`(c("Contig", "Length")) %>%
  mutate(Contig = gsub(">", "", Contig)) %>%
  mutate(Length = as.numeric(Length)) %>%
  right_join(dog, by = "Contig") %>%
  mutate(End = Start + AlignLength)
#.[c(1,3,7)]

#~~~~~~~~~~~~
# Join both dataframes
# Filter by alignment length


align <- cbind(seal_align, dog_align) %>%
  .[c(1,2,3,4,7,8,9,10,14)] %>%
  `colnames<-`(c("Chromosome", "SealContigLength", "chromStart", "AlignLength", "chromEnd", 
                 "Chromosome.1", "DogContigLength", "chromStart.1", "chromEnd.1"))

# get largest 40 seal contigs

top_chr <- dplyr::select(align, Chromosome, SealContigLength) %>%
  distinct(Chromosome, SealContigLength) %>%
  dplyr::arrange(desc(SealContigLength)) %>%
  top_n(40, SealContigLength)

summary(top_chr$SealContigLength)

align <- top_chr %>%
  left_join(align, by = "Chromosome") %>%
  .[-c(2,3,8)] %>%
  mutate(Chromosome = gsub("Contig", "chr", Chromosome)) %>%
  filter(grepl("^NC", Chromosome.1)) %>%
  filter(AlignLength >= 5000) %>% # filter for alignment length
  dplyr::select(-AlignLength)



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
seal_ideogram <- read.table("data/raw/PBJelly_v1.4_lengths", header = F) %>%
  mutate(V1 = gsub("_quiver_pilon", "", V1)) %>%
  `colnames<-`(c("Chromosome", "ChromEnd")) %>%
  mutate(ChromStart = 0) %>%
  dplyr::arrange(desc(ChromEnd)) %>%
  top_n(40, ChromEnd) %>%       
  .[c(1,3,2)] %>%
  mutate(Band = "qA1", Stain = rep(c("gneg", "gpos"),20)) %>%
  mutate(Chromosome = gsub("Contig", "chr", Chromosome))


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


#~~~~~~~~~~~~

# scale seal scaffolds three times larger to that both span around half of plot

link.data <- link.data %>%
  mutate(chromStart = chromStart * 3, chromEnd = chromEnd * 3)

seal_ideogram <- seal_ideogram %>%
  mutate(ChromEnd = ChromEnd * 3)

#~~~~~~~~~~~~



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
params <- RCircos.Get.Plot.Parameters()
params$base.per.unit <- 300000
params$highlight.width <- 0.7
params$line.color <- "gray"
RCircos.Reset.Plot.Parameters(params)

# Initialize graphic device (GUI or image file)
png(file="figs/circos_top40.png", units = "in", res = 300, height=8, width=8)

par(cex=0.75)
RCircos.Set.Plot.Area()
# Plot chromosome ideogram and link lines
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(link.data, track.num, FALSE)
#RCircos.Long.Genome.Link.Plot(link.data, track.num = 1)

dev.off()


jpeg(file="figs/circos_top40.jpg", units = "in", res = 300, height=8, width=8)

par(cex=0.75)
RCircos.Set.Plot.Area()
# Plot chromosome ideogram and link lines
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(link.data, track.num, FALSE)
#RCircos.Long.Genome.Link.Plot(link.data, track.num = 1)

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Alignment Stats         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

detach(package:plyr)

# Using align df which includes NC dog chrosomes and longest 40 seal scaffs
# Contigs mapping to just one unique Dog chromsome

nchr <- align %>%
  dplyr::group_by(Chromosome) %>%
  summarise(length(unique(Chromosome.1)))

nrow(filter(nchr, `length(unique(Chromosome.1))` == 1)) # 26
nrow(filter(nchr, `length(unique(Chromosome.1))` == 1)) / nrow(nchr)  # 65 %

# Contigs mapping mostly to one Dog chromosome

# number of hits to each dog chr
all_hits <- align %>%
  dplyr::group_by(Chromosome, Chromosome.1) %>%
  summarise(n = length(Chromosome.1))

# number of total hits

total_hits <- all_hits %>%
  group_by(Chromosome) %>%
  summarise(sum = sum(n))

# proportion of total hits for each dog chr

# largest number of hits

biggest_hit <- all_hits %>% 
  group_by(Chromosome) %>%
  filter(n == max(n))

df <- left_join(total_hits, biggest_hit, by = "Chromosome") %>%
  mutate(prop = n/sum)

nrow(filter(df, prop > 0.99))
nrow(filter(df, prop > 0.90))
nrow(filter(df, prop > 0.80))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Alignment divergence         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

genome_divergence <- fread("~/alignments/LAST/PBJelly_v1.4_CanFam3.1/maffilter/PBJelly_v1.4_CanFam3.1.statistics.csv")
colnames(genome_divergence) <- c("Chr", "Start", "Stop", "Div")
genome_divergence <- filter(genome_divergence, Div != Inf)
mean(genome_divergence$Div) # 13.8


