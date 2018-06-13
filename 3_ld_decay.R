# Estimating LD decay

library(data.table)
library(dplyr)
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Prepare files for PLINK subsetting          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Prepare random list of seals to remove for even pop (including all pups)

set.seed(1)
ran <- sample(1:57, 51)

pups <- fread("data/raw/pops/SSB.txt", header = F) %>%
  filter(grepl("AGP", V1))

set.seed(1)
random <- fread("data/raw/pops/SSB.txt", header = F) %>%
  filter(!grepl("AGP", V1)) %>%
  sample_n(27, replace = F) # random 27 individuals to remove

remove <- rbind(random, pups)

system("mkdir data/processed")
write.table(remove[,2], "data/processed/random51_IDs.txt",
            col.names = F, row.names = F, quote = F)


#~~ Unrelated individuals

# Write list of pups to remove in fam format

pups <- read.table("data/raw/pops/SSB.txt", header = F) %>%
  filter(grepl("^AGP", V2))


write.table(pups, "data/processed/pups.txt", col.names = F,
            row.names = F, quote = F)

# Write list of adults to keep in fam format

adults <- read.table("data/raw/pops/SSB.txt", header = F) %>%
  filter(!grepl("^AGP", V2))


write.table(adults, "data/processed/adults.txt", col.names = F,
            row.names = F, quote = F)


#~~ Even population size (used for structure analysis)

# Get random seals and add second column
random51 <- read.table("data/processed/random51_IDs.txt", header = F, col.names = "V1") %>%
  mutate(V2 = V1)

write.table(random51, "data/processed/random51.txt", col.names = F,
            row.names = F, quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Prepare raw vcf files           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Filter for high depth of coverage
system("mkdir data/ld")

system("vcftools --vcf data/ArcGaz_pilon_biallelic_sub_maf.recode.vcf --out data/ld/GQ_DP --minGQ 5 --minDP 8 --maxDP 30 --recode --recode-INFO-all")
system("vcftools --vcf data/ld/GQ_DP.recode.vcf --plink --out data/ld/GQ_DP")

# recode and reorder map
system("scripts/recode_full_map.sh data/ld/GQ_DP.map data/ld/GQ_DP.map")
system("plink --file data/ld/GQ_DP --make-bed --out data/ld/GQ_DP --allow-extra-chr")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Select SNPs in top 100 contigs    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

lengths <- fread("data/raw/PBJelly_v1.4_lengths") %>%
  arrange(desc(V2)) %>%
  .[c(1:100),]

sum(lengths$V2)
summary(lengths$V2)

snps <- fread("data/ld/GQ_DP.bim") %>%
  select(V1, V4) %>%
  right_join(lengths, by = "V1") %>%
  unite(ID, V1, V4, sep = ":")

write.table(snps[,1], "data/ld/snps_top40.txt",
            row.names = F, col.names = F, quote = F)


# Extract SNPs in longest contigs from PLINK files
system("plink --bfile data/ld/GQ_DP --extract data/ld/snps_top100.txt --nonfounders --recode --make-bed --out data/ld/top100snps --allow-extra-chr")

# Remove SNPs with Mendel Errors
system("plink --bfile data/ld/top100snps --exclude data/mendel/MendelErrorStrict_SNPs.txt --nonfounders --recode --make-bed --out data/ld/top100snpsME --allow-extra-chr")




#~~ Subset files


# Paths to pops
pop_files <- paste("data/raw/pops/", list.files(path = "data/raw/pops", pattern="*.txt"), sep = "")
pop_files

out_names <- gsub(".txt", "", list.files(path = "data/raw/pops"))


# Filter for 50% genotyping rate and 10% MAF
# geno 0.01 = 99%
# geno 0.05 = 95%
# geno 0.1 = 90%


for (i in 1:length(pop_files)){
  system(paste0("plink --bfile data/ld/top100snpsME --keep ", pop_files[i], " --nonfounders --geno 0.5 --maf 0.1 --hwe 0.001 --make-bed --out data/ld/", out_names[i], " --allow-extra-chr"))
}


# Filter subsetted SSB files in the same way

system("plink --bfile data/ld/SSB --remove data/processed/pups.txt --nonfounders --geno 0.5 --maf 0.1 --hwe 0.001 --make-bed --out data/ld/SSB_unrelated --allow-extra-chr")
system("plink --bfile data/ld/SSB --remove data/processed/random51.txt --nonfounders --geno 0.5 --maf 0.1 --hwe 0.001 --make-bed --out data/ld/SSB_even --allow-extra-chr")

system("plink --bfile data/ld/SSB --nonfounders --geno 0.35 --make-bed --out data/ld/SSB_strict --allow-extra-chr")

# adults

system("plink --bfile data/ld/top100snpsME --keep data/processed/adults.txt --nonfounders --geno 0.5 --maf 0.1 --hwe 0.001 --make-bed --out data/ld/SSB_adults --allow-extra-chr")


# Unrelated pups

unrelated <- read.table("data/processed/adults_fullSibs.txt") %>%
  mutate(V2 = V1)

write.table(unrelated, "data/processed/adults_fullSibs_plink.txt", 
            col.names = F, row.names = F, quote = F)

system("plink --bfile data/ld/SSB --remove data/processed/adults_fullSibs_plink.txt --nonfounders --geno 0.5 --maf 0.1 --hwe 0.001 --make-bed --out data/ld/SSB_unpups --allow-extra-chr")


#~~~~~~~~~~~~~~~~~~~~~~~#
#     Run PLINK LD      #
#~~~~~~~~~~~~~~~~~~~~~~~#

out_names <- gsub(".txt", "", list.files(path = "data/raw/pops"))

# --ld-window-kb 50 *
# --ld-window 50
# --ld-window 10 = huisman pnas
# --ld-window-kb 1000 * polar bear

for (i in 1:length(pop_files)){
  system(paste0("plink --bfile data/ld/",out_names[i], " --r2 --ld-window-kb 500 --allow-extra-chr --debug --out data/ld/", out_names[i], " --ld-window-r2 0"))
}

system("plink --bfile data/ld/SSB_unrelated --r2 --ld-window-kb 500 --allow-extra-chr --debug --out data/ld/SSB_unrelated --ld-window-r2 0")
system("plink --bfile data/ld/SSB_even --r2 --ld-window-kb 500 --allow-extra-chr --nonfounders --debug --out data/ld/SSB_even --ld-window-r2 0")
system("plink --bfile data/ld/SSB_strict --r2 --ld-window-kb 500 --allow-extra-chr --nonfounders --debug --out data/ld/SSB_strict --ld-window-r2 0")
system("plink --bfile data/ld/SSB_unpups --r2 --ld-window-kb 500 --allow-extra-chr --nonfounders --debug --out data/ld/SSB_unpups --ld-window-r2 0")
system("plink --bfile data/ld/SSB_adults --r2 --ld-window-kb 500 --allow-extra-chr --nonfounders --debug --out data/ld/SSB_adults --ld-window-r2 0")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Visualise LD Decay      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# IMPORT the data

import_files <- function(file){
  file <- fread(file) %>%
    mutate(distance = abs(BP_B - BP_A)) %>%
    arrange(distance)
}


# paths to datasets

ld <- paste("data/ld/", list.files(path = "data/ld", pattern="*.ld"), sep = "")

ld_out <- lapply(ld, import_files)

dataset_names <-list.files(path = "data/ld", pattern="*.ld")
dataset_names <-  gsub(".ld", "", dataset_names)
names(ld_out) <- dataset_names
names(ld_out) 

library(plyr)
ld_out_df <- ldply(ld_out)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          LD curve            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Nonlinear logarithmic regression curve of r2 on genetic distance
#~~ Source: https://fabiomarroni.wordpress.com/2011/08/09/estimate-decay-of-linkage-disequilibrium-with-distance/

nonlinear_ld_curve <- function(r2){
  r2 <- r2 %>%
    mutate(distanceKb = distance/1000)
  
  distance <- r2$distanceKb
  LD.data <- r2$R2
  n <- 33*2

  HW.st<-c(C=0.1)
  HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*
                      (1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),
                    start=HW.st,control=nls.control(maxiter=100))
  
  tt<-summary(HW.nonlinear)
  new.rho<-tt$parameters[1]
  fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*
    (1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
  
  ld.df<-data.frame(distance,fpoints)
  ld.df<-ld.df[order(ld.df$distance),]
  ld.df <- ld.df
}


ld_curves <- lapply(ld_out, nonlinear_ld_curve)
ld_curves <- ldply(ld_curves)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Plot results          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggthemr)
ggthemr(palette = "pale", layout = "clean", 
        line_weight = 0.5, text_size = 16, type = "outer", spacing = 2)
cols <- c("#66A61E", "#E6AB02", "#D95F02", "#7570B3", "#E7298A", "#1B9E77", "#A6761D")

library(RColorBrewer)
cols <- c(brewer.pal(8, "Dark2"))


#~~ SSB adults for manuscript

png(file="figs/LD_decay_SSB_adults.png", units = "in", res = 300, height=6, width=7)

ggplot(filter(ld_out_df, .id == "SSB_adults"), aes(x=distance, y = R2)) +
  geom_point(aes(x = distance/1000, y = R2), size = 0.1, alpha = 1/5, col = "darkgrey") + # 1/1.5
  geom_line(data = filter(ld_curves, .id == "SSB_adults"), aes(x = distance, y = fpoints), size = 1, col = "grey30") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  labs(x = "Distance between SNPs (kb)", y = expression(paste("Linkage disequilibrium ", r^2))) +
  #ylim(c(0,0.5)) + xlim(c(0,1000)) # 1000 kb?
  ylim(c(0,0.5)) + xlim(c(0,500))

dev.off()

#~~ Other subsets:

#~~ All pops

png(file="figs/ld_decay_allpops.png", units = "in", res = 300, height=9, width=10)

ggplot(ld_out_df, aes(x=distance, y = R2)) +
  #geom_point(aes(x = distance/1000, y = R2), size = 0.1, alpha = 1/1.5, col = "lightgrey") +
  geom_line(data = ld_curves, aes(x = distance, y = fpoints, colour = factor(.id)), size = 1) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  labs(x = "Distance (kb)", y = expression(bold(paste("Linkage disequilibrium ", r^2)))) +
  scale_color_manual(values = cols, 
                     name="Population",
                     labels = c("Bouvetoya", "South Shetlands", "Heard Island", "Iles Kerguelen", "Macquarie", "SSB", "SSB even", "SSB unrelated"))

dev.off()

#~~ SSB unrelated

png(file="figs/LD_decay_SSB_unrelated.png", units = "in", res = 300, height=6, width=7)

ggplot(filter(ld_out_df, .id == "SSB_unrelated"), aes(x=distance, y = R2)) +
  geom_point(aes(x = distance/1000, y = R2), size = 0.1, alpha = 1/5, col = "darkgrey") + # 1/1.5
  geom_line(data = filter(ld_curves, .id == "SSB_unrelated"), aes(x = distance, y = fpoints), size = 1, col = "grey30") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  labs(x = "Distance between SNPs (kb)", y = expression(paste("Linkage disequilibrium ", r^2))) +
  ylim(c(0,0.5)) + xlim(c(0,500))

dev.off()



#~~ SSB

png(file="figs/LD_decay_SSB.png", units = "in", res = 300, height=6, width=7)

ggplot(filter(ld_out_df, .id == "SSB"), aes(x=distance, y = R2)) +
  geom_point(aes(x = distance/1000, y = R2), size = 0.1, alpha = 1/5, col = "darkgrey") + # 1/1.5
  geom_line(data = filter(ld_curves, .id == "SSB"), aes(x = distance, y = fpoints), size = 1, col = "grey30") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  labs(x = "Distance between SNPs (kb)", y = expression(paste("Linkage disequilibrium ", r^2))) +
  ylim(c(0,0.5)) + xlim(c(0,500))

dev.off()


#~~ SSB strict genotyping rate 75%

png(file="figs/LD_decay_SSB_strict.png", units = "in", res = 300, height=6, width=7)

ggplot(filter(ld_out_df, .id == "SSB_strict"), aes(x=distance, y = R2)) +
  geom_point(aes(x = distance/1000, y = R2), size = 0.1, alpha = 1/5, col = "darkgrey") + # 1/1.5
  geom_line(data = filter(ld_curves, .id == "SSB"), aes(x = distance, y = fpoints), size = 1, col = "grey30") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  labs(x = "Distance between SNPs (kb)", y = expression(paste("Linkage disequilibrium ", r^2))) +
  ylim(c(0,0.5)) + xlim(c(0,500))

dev.off()


#~~~~~~~~~~~~~~~~~~~~#
#       Stats        #
#~~~~~~~~~~~~~~~~~~~~#

ssb_ld <- filter(ld_out_df, .id == "SSB_adults")
mean(ssb_ld$R2)

head(filter(ld_curves, .id == "SSB_adults") %>%
  arrange(desc(fpoints)) %>%
  filter(fpoints <= 0.2))

head(filter(ld_curves, .id == "SSB_adults") %>%
       arrange(desc(fpoints)) %>%
       filter(fpoints <= 0.5))

head(filter(ld_curves, .id == "SSB") %>%
       arrange(desc(fpoints)) %>%
       filter(fpoints <= 0.2))

