# structure analysis on pups not adults

# Population structure analysis
# PCA using adegenet

library(adegenet)
library(ggplot)
library(ggthemr)
library(dplyr)
library(data.table)
library(tidyr)
library(gridExtra)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Remove adults and related pups & MAF / geno filter      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Read in SSB id's (not in fam format)

adults <- fread("data/raw/pops/SSB.txt", header = F) %>%
  filter(!grepl("AGP", V1))

# Identify full and half sibs

full_sibs <- fread("data/raw/pops/SSB.txt", header = F) %>%
  filter(grepl("AGP", V1)) %>%
  filter(grepl("AGP03049|AGP98046|AGP01087|AGP01036|AGP97025|
               AGP02078|AGP02094|AGP03013|AGP01109|AGP96056|AGP96155|
               AGP99059|AGP01088|AGP97086|AGP00091", V1))

remove <- rbind(adults, full_sibs)

write.table(remove[,2], "data/processed/adults_fullSibs.txt",
            col.names = F, row.names = F, quote = F)

# structure_pups

system("mkdir data/structure")
system("vcftools --vcf data/ArcGaz_pilon_biallelic_sub.recode.vcf --remove data/processed/adults_fullSibs.txt --maf 0.05 --max-missing 0.99 --recode --recode-INFO-all --out data/structure/maf05_g99")


# hardy weinberg equilibrium

system("vcftools --vcf data/structure/maf05_g99.recode.vcf --hwe 0.01 --recode --recode-INFO-all --out data/structure/maf05_g99_hwe")

system("vcftools --vcf data/structure/maf05_g99.recode.vcf --keep data/raw/pops/bouvetoya.txt --hardy --out data/structure/maf05_g99_hweBouv")
system("vcftools --vcf data/structure/maf05_g99.recode.vcf --keep data/raw/pops/cape_shirreff.txt --hardy --out data/structure/maf05_g99_hweCS")
system("vcftools --vcf data/structure/maf05_g99.recode.vcf --keep data/raw/pops/heard_island.txt --hardy --out data/structure/maf05_g99_hweHI")
system("vcftools --vcf data/structure/maf05_g99.recode.vcf --keep data/raw/pops/iles_kerguelen.txt --hardy --out data/structure/maf05_g99_hweIK")
system("vcftools --vcf data/structure/maf05_g99.recode.vcf --keep data/raw/pops/macquarie_island.txt --hardy --out data/structure/maf05_g99_hweMI")
system("vcftools --vcf data/structure/maf05_g99.recode.vcf --keep data/raw/pops/SSB.txt --hardy --out data/structure/maf05_g99_hweSG")


hw_b <- fread("data/structure/maf05_g99_hweBouv.hwe")
hw_cs <- fread("data/structure/maf05_g99_hweCS.hwe")
hw_hi <- fread("data/structure/maf05_g99_hweHI.hwe")
hw_ik <- fread("data/structure/maf05_g99_hweIK.hwe")
hw_mi <- fread("data/structure/maf05_g99_hweMI.hwe")
hw_sg <- fread("data/structure/maf05_g99_hweSG.hwe")

nrow(filter(hw_b, P_HWE > 0.01))
nrow(filter(hw_cs, P_HWE > 0.01))
nrow(filter(hw_hi, P_HWE > 0.01))
nrow(filter(hw_ik, P_HWE > 0.01))
nrow(filter(hw_mi, P_HWE > 0.01))
nrow(filter(hw_sg, P_HWE > 0.01))

# file of snps out of hwe

write.table(filter(hw_sg, P_HWE < 0.01) %>%
              select(CHR, POS) %>%
              unite(X, c(CHR, POS), sep = "\t"),
            "data/structure/hwe_filter.txt", 
            quote = F, row.names = F, col.names = F)

# filter for SNPs out of HWE

system("vcftools --vcf data/structure/maf05_g99.recode.vcf --exclude-positions data/structure/hwe_filter.txt --recode --recode-INFO-all --out data/structure/maf05_g99_hwe")


# Make plink

system("vcftools --vcf data/structure/maf05_g99_hwe.recode.vcf --plink --out data/structure/maf05_g99_hwe")
system("plink --file data/structure/maf05_g99_hwe --make-bed --recode --out data/structure/maf05_g99_hwe")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Linkage pruning          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# recode map file for LD pruning
system("scripts/recode_full_map.sh data/structure/maf05_g99_hwe.map data/structure/maf05_g99_hwe.map")

# Make bed
system("plink --file data/structure/maf05_g99_hwe --make-bed --out data/structure/maf05_g99_hwe --allow-extra-chr --debug")

# Run PLINK --indep
system("plink --bfile data/structure/maf05_g99_hwe --indep 50 5 2 --nonfounders --out data/structure/maf05_g99_hwe --allow-extra-chr --debug")
system("wc -l data/structure/maf05_g99_hwe.prune.in")

# LD filter
# Make raw file for adegenet
system("plink --bfile data/structure/maf05_g99_hwe --extract data/structure/maf05_g99_hwe.prune.in --out data/structure/maf05_g99_hwe_ld --recodeA --allow-extra-chr --debug")

# Made ped file for structure
system("plink --bfile data/structure/maf05_g99_hwe --extract data/structure/maf05_g99_hwe.prune.in --out data/structure/maf05_g99_hwe_ld --recode --allow-extra-chr --debug")


#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       adegenet PCA       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

gl <- read.PLINK("data/structure/maf05_g99_hwe_ld.raw")

# assign population colours

library(RColorBrewer)

cbPalette <- c(brewer.pal(6, "Dark2"))
col <- ifelse(grepl("^H", gl@pop), cbPalette[2],
              ifelse(grepl("^K", gl@pop), cbPalette[3],
                     ifelse(grepl("^M", gl@pop), cbPalette[4],
                            ifelse(grepl("^B", gl@pop), cbPalette[5],
                                   ifelse(grepl("^I", gl@pop) | grepl("^i", gl@pop), cbPalette[6], cbPalette[1])))))

unique(col)

# PCA

pca1 <- glPca(gl, returnDotProd=T) # naxes 10
4


# ggplot

ggthemr(palette = "pale", layout = "clean", 
        line_weight = 0.7, text_size = 20, type = "outer")

pc1 <- pca1$scores[,1]
pc2 <- pca1$scores[,2]
pc3 <- pca1$scores[,3]
pc4 <- pca1$scores[,4]
ind_names <- gl@ind.names

ggplot_pca <- as.data.frame(cbind(pc1, pc2, pc3, pc4)) %>%
  mutate(ind_names = ind_names) %>%
  mutate(pop = col) 

# eig

eig <- data.frame(pca1$eig)
eig$percentage = (eig[, 1]/sum(eig$pca1.eig))*100
sum(eig$percentage)
sum(eig$percentage[1:2])

eig$percentage <- round(eig$percentage, digits = 1)
eig$percentage[1]
eig$percentage[2]
eig$percentage[3]


rad_pca_1_2 <- ggplot(ggplot_pca, aes(pc1, pc2)) + 
  geom_point(aes(colour = factor(pop)), size = 3) +
  scale_color_manual(values = cbPalette, 
                     name="Population",
                     breaks=c("#1B9E77", "#E6AB02", "#66A61E", "#7570B3", "#D95F02", "#E7298A"),
                     labels=c("South Georgia", "South Shetlands", "Bouvetoya", "Kerguelen", "Heard Island", "Macquarie")) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        legend.position = "none") +
  ylab(paste0("PC2 (",eig$percentage[2],"%)")) +
  xlab(paste0("PC1 (",eig$percentage[1],"%)")) +
  ggtitle('(B) 27,592 SNPs') + 
  theme(plot.title=element_text(hjust=0, size = 20, face = "plain")) +
  scale_x_reverse()



# PCA 1 & 3

rad_pca_1_3 <- ggplot(ggplot_pca, aes(pc1, pc3)) + 
  geom_point(aes(colour = factor(pop)), size = 3) +
  scale_color_manual(values = cbPalette, 
                     name="Population",
                     breaks=c("#1B9E77", "#E6AB02", "#66A61E", "#7570B3", "#D95F02", "#E7298A"),
                     labels=c("South Georgia", "South Shetlands", "Bouvetoya", "Kerguelen", "Heard Island", "Macquarie")) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        legend.position = "none") +
  ylab(paste0("PC3 (",eig$percentage[3],"%)")) +
  xlab(paste0("PC1 (",eig$percentage[1],"%)")) +
  ggtitle('(D) 27,592 SNPs') + 
  theme(plot.title=element_text(hjust=0, size = 20, face = "plain")) +
  scale_x_reverse()



# Plot to get legend
png("figs/pca_legend.png", units = "in", res = 300, width = 8, height = 6)

ggplot(ggplot_pca, aes(pc1, pc2)) + 
  geom_point(aes(colour = factor(pop)), size = 3) +
  scale_color_manual(values = cbPalette, 
                     name="Population",
                     breaks=c("#1B9E77", "#E6AB02", "#66A61E", "#7570B3", "#D95F02", "#E7298A"),
                     labels=c("South Georgia", "South Shetlands", "Bouvetoya", "Kerguelen", "Heard Island", "Macquarie")) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        legend.title = element_text(face = "plain")) +
  ylab("PC2") +
  xlab("PC1")

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Make structure file   #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

# Handmade structure file

# get genotypes
system("mkdir data/structure/temp")
system("cut -d ' ' -f 7- data/structure/maf05_g99_hwe_ld.ped > data/structure/temp/A")

# convert values
system("sed 's/A/1/g' data/structure/temp/A > data/structure/temp/B")
system("sed 's/T/2/g' data/structure/temp/B > data/structure/temp/C")
system("sed 's/G/3/g' data/structure/temp/C > data/structure/temp/D")
system("sed 's/C/4/g' data/structure/temp/D > data/structure/temp/E")
system("sed 's/0/-9/g' data/structure/temp/E > data/structure/temp/F")

# get ids
system("awk '{print $1}' data/structure/maf05_g99_hwe_ld.ped > data/structure/temp/ids")


# make structure input file
gen <- fread("data/structure/temp/F", header = F)
ids <- fread("data/structure/temp/ids", header = F)
pop <- rep(1, nrow(ids))
stru_in <- cbind(ids, pop)
stru_in <- stru_in %>% mutate(pop = ifelse(grepl("^H", V1), 5,
                                           ifelse(grepl("^K", V1), 4,
                                                  ifelse(grepl("^M", V1), 6,
                                                         ifelse(grepl("^B", V1), 3,
                                                                ifelse(grepl("^I", V1) | grepl("^i", V1), 2, 1))))))

write.table(stru_in, "data/structure/temp/ids_pop.txt", col.names = F, row.names = F, quote = F)

system("paste -d ' ' data/structure/temp/ids_pop.txt data/structure/temp/F > data/structure/handmade_rad.stru")


# sort structure files by population

system("sort -k 2 data/structure/handmade_rad.stru > data/structure/handmade_rad_sort.stru")

system("rm -r data/structure/temp")





#~~~~~~~~~~~~~~~~~~~~~~~~~#
#     Microsatellites     #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Get microsatellite summary stats for Supp

ids <- read.table("data/raw/individuals_to_keep.txt") %>%
  dplyr::select(V1)

ms_plate <- fread("data/raw/Rack61_microsatellites_Jul2017.csv", header = T) %>%
  right_join(ids, by = c("ID" = "V1"))

pop_ids <- ms_plate %>%
  dplyr::select(ID) %>%
  mutate(pop = ifelse(grepl("^H", ID), 5,
                      ifelse(grepl("^K", ID), 4,
                             ifelse(grepl("^M", ID), 6,
                                    ifelse(grepl("^B", ID), 3,
                                           ifelse(grepl("^I", ID) | grepl("^i", ID), 2, 1))))))

ms_plate <- cbind(pop_ids, ms_plate) %>%
  .[-3]  %>% # rm third col
  arrange(pop)

write.table(ms_plate, "data/processed/ms.stru", col.names = F, row.names = F, quote = F, sep = " ")

file <- "data/processed/ms.stru"
gen <- fread("data/processed/ms.stru", header = F)

struc <- read.structure(file = file,
                        n.ind = nrow(gen),
                        n.loc = (ncol(gen)-2)/2,
                        onerowperind = T,
                        col.lab = 0,
                        col.pop = 2,
                        col.others = 1, 
                        row.marknames = 0,
                        ask = F)

#~~ ms summary
div <- summary(struc)
ms_df <- data.frame(div$Hobs, div$Hexp, div$loc.n.all)

#~~ allele size range:
ms_sr <- ms_plate %>%
  dplyr::select(-c(ID, pop))

ms_sr[ms_sr == -9] <- NA

# stack alleles
library(stringr)
# remove .1 and .2
names(ms_sr) <- substring(names(ms_sr),1,str_length(names(ms_sr))-2)
head(ms_sr)

# melt same alleles
ms_sr <- sapply(unique(names(ms_sr)), function(x) unname(unlist(ms_sr[,names(ms_sr)==x])))

# col min and max
min <- as.data.frame(apply(ms_sr,2,min, na.rm = T))
max <- as.data.frame(apply(ms_sr,2,max, na.rm = T))
locus_names <- rownames(max)
ms_sr <- cbind(min, max)
colnames(ms_sr) <- c("min", "max")

ms_df <- cbind(ms_df, ms_sr)
ms_df <- mutate(ms_df, locus = locus_names)
write.csv(ms_df, "data/processed/microsatellite_summary_stats.csv", row.names = F, quote = F)



#~~ Get microsatellite structure input file

ms_plate <- fread("data/raw/Rack61_microsatellites_Jul2017.csv", header = T) %>%
  right_join(ids, by = c("ID" = "V1"))

pop_ids <- ms_plate %>%
  select(ID) %>%
  mutate(pop = ifelse(grepl("^H", ID), 5,
                      ifelse(grepl("^K", ID), 4,
                             ifelse(grepl("^M", ID), 6,
                                    ifelse(grepl("^B", ID), 3,
                                           ifelse(grepl("^I", ID) | grepl("^i", ID), 2, 1))))))

ms_plate <- cbind(pop_ids, ms_plate) %>%
  .[-3]  %>% # rm third col
  arrange(pop)

system("mkdir data/structure/ms")
write.table(ms_plate, "data/structure/ms/stru_in.stru", col.names = F, row.names = F, quote = F, sep = " ")


#~~~~~~~~~~~~~~~~~~~#
#        PCA        #
#~~~~~~~~~~~~~~~~~~~#

file <- "data/structure/ms/stru_in.stru"
gen <- fread("data/structure/ms/stru_in.stru", header = F)


struc <- read.structure(file = file,
                        n.ind = nrow(gen),
                        n.loc = (ncol(gen)-2)/2,
                        onerowperind = T,
                        col.lab = 0,
                        col.pop = 2,
                        col.others = 1, 
                        row.marknames = 0,
                        ask = F)



sum(is.na(struc$tab))
X <- scaleGen(struc, NA.method="mean")

# dudi.pca object
pca2 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)

#png(file = paste0(out, "_PCA.png"), units = "in", res = 300, height = 8, width = 11)
# s.class(pca2$li, pop(struc), xax = 1, yax = 2, col = cbPalette, grid = T, clabel = 0.7, addaxes = T)

xaxis1 <- pca2$li$Axis1
yaxis_2 <- pca2$li$Axis2
yaxis_3 <- pca2$li$Axis3

ggplot_pca_ms <- as.data.frame(cbind(xaxis1, yaxis_2, yaxis_3)) %>%
  mutate(pop = pop(struc)) %>%
  mutate(col = ifelse(grepl(5, pop), "#D95F02",
                      ifelse(grepl(4, pop), "#7570B3",
                             ifelse(grepl(6, pop), "#E7298A",
                                    ifelse(grepl(3, pop), "#66A61E",
                                           ifelse(grepl(2, pop), "#E6AB02", "#1B9E77")))))) %>%
  mutate(as.factor(col))


# eig

ms_eig <- data.frame(pca2$eig)
ms_eig$percentage = (ms_eig[, 1]/sum(ms_eig$pca2.eig))*100
sum(ms_eig$percentage)
sum(ms_eig$percentage[1:2])

ms_eig$percentage <- round(ms_eig$percentage, digits = 1)
ms_eig$percentage[1]
ms_eig$percentage[2]
ms_eig$percentage[3]


#~~~~~~~~~~~~~~~~#
#     Plot       #
#~~~~~~~~~~~~~~~~#

cbpal2 <- cbPalette[c(1,5,2,3,4,6)]

# pca 1 & 2

#png("figs/ms_pca1_2.png", units = "in", res = 300, width = 8, height = 6)

ms_pca_1_2 <- ggplot(ggplot_pca_ms, aes(xaxis1, yaxis_2)) + 
  geom_point(aes(colour = ggplot_pca_ms$pop), size = 3) +
  scale_color_manual(values = cbpal2,
                     name="Population",
                     breaks=c("1", "2", "3", "4", "5", "6"),
                     labels=c("South Georgia", "South Shetlands", "Bouvetoya", "Kerguelen", "Heard Island", "Macquarie")) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        legend.position="none") +
  ylab(paste0("PC2 (",ms_eig$percentage[2],"%)")) +
  xlab(paste0("PC1 (",ms_eig$percentage[1],"%)")) +
  ggtitle('(A) 27 microsatellites') + 
  theme(plot.title=element_text(hjust=0, size = 20, face = "plain"))




# pca 1 & 3


ms_pca_1_3 <- ggplot(ggplot_pca_ms, aes(xaxis1, yaxis_3)) + 
  geom_point(aes(colour = ggplot_pca_ms$pop), size = 3) +
  scale_color_manual(values = cbpal2,
                     name="Population",
                     breaks=c("1", "2", "3", "4", "5", "6"),
                     labels=c("South Georgia", "South Shetlands", "Bouvetoya", "Kerguelen", "Heard Island", "Macquarie")) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        legend.position="none") +
  ylab(paste0("PC3 (",ms_eig$percentage[3],"%)")) +
  xlab(paste0("PC1 (",ms_eig$percentage[1],"%)")) +
  ggtitle('(C) 27 microsatellites') + theme(plot.title=element_text(hjust=0, size = 20, face = "plain"))


#~~ Plot for manuscript


png("figs/pca_1_2.png", units = "in", res = 300, width = 11, height = 6)
grid.arrange(ms_pca_1_2, rad_pca_1_2, ncol=2, nrow=1)
dev.off()


png("figs/pca_1_3.png", units = "in", res = 300, width = 11, height = 6)
grid.arrange(ms_pca_1_3, rad_pca_1_3, ncol=2, nrow=1)
dev.off()


png("figs/pca_1_2_3.png", units = "in", res = 300, width = 11, height = 11)
grid.arrange(ms_pca_1_2, rad_pca_1_2,
             ms_pca_1_3,rad_pca_1_3, ncol=2, nrow=2)
dev.off()


jpeg("figs/pca_1_2_3.jpg", units = "in", res = 300, width = 11, height = 11)
grid.arrange(ms_pca_1_2, rad_pca_1_2,
             ms_pca_1_3,rad_pca_1_3, ncol=2, nrow=2)
dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       adegenet DAPC       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(RColorBrewer)
cbPalette <- c(brewer.pal(6, "Dark2"))

gl <- read.PLINK("data/structure/maf05_g99_hwe_ld.raw")

# specify colours for plotting later
col <- ifelse(grepl("^H", gl@pop), cbPalette[2],
              ifelse(grepl("^K", gl@pop), cbPalette[3],
                     ifelse(grepl("^M", gl@pop), cbPalette[4],
                            ifelse(grepl("^B", gl@pop), cbPalette[5],
                                   ifelse(grepl("^I", gl@pop) | grepl("^i", gl@pop), cbPalette[6], cbPalette[1])))))

# assign geographic locations to gl object
gl$pop <- as.factor(ifelse(grepl("^H", gl@pop), "HI",
                           ifelse(grepl("^K", gl@pop), "K",
                                  ifelse(grepl("^M", gl@pop), "Mac",
                                         ifelse(grepl("^B", gl@pop), "B",
                                                ifelse(grepl("^I", gl@pop) | grepl("^i", gl@pop), "CS", "SG"))))))



#~~ Find clusters

# Best BIC is 0
grp <- find.clusters(gl, max.n.clust = 15)
60 # Number PCs
10 # Number clusters

# Run with pop info, npcs 19 (n/3)
dapc1 <- dapc(gl, gl$pop, n.da=60, n.pca=19)

# find optimal n pcas to keep
temp <- optim.a.score(dapc1)
temp$best
dapc2 <- dapc(gl, gl$pop, n.pca = temp$best, n.da = 10) # naxes 10


#~~ Plot results

ggthemr(palette = "pale", layout = "clean", 
        line_weight = 0.7, text_size = 20, type = "outer")

# Get DAPC coords
pc1 <- dapc2$ind.coord[,1]
pc2 <- dapc2$ind.coord[,2]
pc3 <- dapc2$ind.coord[,3]
pc4 <- dapc2$ind.coord[,4]

ggplot_dapc <- as.data.frame(cbind(pc1, pc2, pc3, pc4)) %>%
  mutate(pop = col) 


png("figs/rad_dapc1_2.png", units = "in", res = 300, width = 10, height = 8)
ggplot(ggplot_dapc, aes(pc1, pc2)) + 
  geom_point(aes(colour = factor(pop)), size = 3) +
  scale_color_manual(values = cbPalette, 
                     name="Population",
                     breaks=c("#1B9E77", "#E6AB02", "#66A61E", "#7570B3", "#D95F02", "#E7298A"),
                     labels=c("South Georgia", "South Shetlands", "Bouvetoya", "Kerguelen", "Heard Island", "Macquarie")) +
  theme(axis.text.x = element_text(face = "plain", colour = "white"),
        axis.text.y = element_text(face = "plain", colour = "white")) +
  ylab("PC2") +
  xlab("PC1")
dev.off()


png("figs/rad_dapc1_3.png", units = "in", res = 300, width = 10, height = 8)
ggplot(ggplot_dapc, aes(pc1, pc3)) + 
  geom_point(aes(colour = factor(pop)), size = 3) +
  scale_color_manual(values = cbPalette, 
                     name="Population",
                     breaks=c("#1B9E77", "#E6AB02", "#66A61E", "#7570B3", "#D95F02", "#E7298A"),
                     labels=c("South Georgia", "South Shetlands", "Bouvetoya", "Kerguelen", "Heard Island", "Macquarie")) +
  theme(axis.text.x = element_text(face = "plain", colour = "white"),
        axis.text.y = element_text(face = "plain", colour = "white")) +
  ylab("PC2") +
  xlab("PC1")
dev.off()
