# Quantifying inbreeding
# sMLH, IBCs

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
options(scipen=999)
library(reshape2)
library(inbreedR)

system("mkdir data/inbreeding")

#~~~~~~~~~~~~~~~~~~~~~~#
#  maf 0.05, 60% geno  #
#~~~~~~~~~~~~~~~~~~~~~~# 

system("vcftools --vcf data/ArcGaz_pilon_biallelic_sub_maf.recode.vcf --out data/inbreeding/maf05_g60 --maf 0.05 --max-missing 0.60 --recode --recode-INFO-all")
system("vcftools --vcf data/inbreeding/maf05_g60.recode.vcf --plink --out data/inbreeding/maf05_g60")
system("plink --file data/inbreeding/maf05_g60 --make-bed --recodeAD --out data/inbreeding/maf05_g60")


#~~ Filter raw vcf file

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  maf 0.05, 60% geno, GQ 5, DP 8  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

system("vcftools --vcf data/ArcGaz_pilon_biallelic_sub_maf.recode.vcf --geno-depth --out data/inbreeding/raw")

depth <- fread("data/inbreeding/raw.gdepth", colClasses = "numeric") %>%
  unite(SNP, CHROM, POS)
depth <- select(depth, -SNP)

depth <- depth %>%
  mutate(mean = rowMeans(.), 
         sum = rowSums(.), 
         logmeandepth = log10(mean))

CI <- 0.90
CI_meanDP <- stats::quantile(depth$mean, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_sumDP <- stats::quantile(depth$sum, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_meanDP
CI_sumDP

hist(depth$logmeandepth)

# Filter for DP
system("vcftools --vcf data/ArcGaz_pilon_biallelic_sub_maf.recode.vcf --out data/inbreeding/GQ_DP --minGQ 5 --minDP 8 --maxDP 30 --recode --recode-INFO-all")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#            Remove individuals with high amounts of missing data          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Identify individuals with lots of missing data

system("vcftools --vcf data/inbreeding/maf05_g60.recode.vcf --missing-indv --out data/inbreeding/maf05_g60")
system("vcftools --vcf data/inbreeding/GQ_DP.recode.vcf --missing-indv --out data/inbreeding/GQ_DP")


# make histogram

import_files <- function(file){
  fread(file)
}

IM_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.imiss"), sep = "")
IM <- lapply(IM_files, import_files)
names(IM) <- IM_files

hist_fun <- function(x){
  hist(x$F_MISS)
}

par(mfrow = c(1, 2))
histograms <- lapply(IM, hist_fun)
par(mfrow = c(1, 1))

# 40% maf05_g60
# 70% GQ_DP_maf05_g30

system("mawk '$5 > 0.4' data/inbreeding/maf05_g60.imiss | cut -f1 > data/inbreeding/maf05_g60_lowDP.indv") # 0.4
system("mawk '$5 > 0.9' data/inbreeding/GQ_DP.imiss | cut -f1 > data/inbreeding/GQ_DP_lowDP.indv")

# feed to vcftools (remove individuals with initial low coverage)
system("vcftools --vcf data/inbreeding/maf05_g60.recode.vcf --remove data/inbreeding/maf05_g60_lowDP.indv --recode --recode-INFO-all --out data/inbreeding/maf05_g60_ldi")
system("vcftools --vcf data/inbreeding/GQ_DP.recode.vcf --remove data/inbreeding/maf05_g60_lowDP.indv --recode --recode-INFO-all --out data/inbreeding/GQ_DP_ldi")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Restrict Dataset        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Now that poor coverage individuals have been removed:
# Remove mendel errors
# Keep variants called in a high percentage of individuals


#~~ Mendel errors

system("vcftools --vcf data/inbreeding/GQ_DP_ldi.recode.vcf --exclude-positions data/mendel/MendelErrorStrict_SNPs_tab.txt --recode --recode-INFO-all --out data/inbreeding/GQ_DP_ldi_ME")


#~~ Filter for SG individuals here

system("vcftools --vcf data/inbreeding/GQ_DP_ldi_ME.recode.vcf --keep data/raw/SSB_IDs.txt --recode --recode-INFO-all --out data/inbreeding/GQ_DP_ldi_ME_SSB")
system("bcftools query -l data/inbreeding/GQ_DP_ldi_ME_SSB.recode.vcf | wc -l ")



#~~ Genotyping rate

system("vcftools --vcf data/inbreeding/GQ_DP_ldi_ME_SSB.recode.vcf --max-missing 0.7 --maf 0.05 --recode --recode-INFO-all --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB")
system("bcftools query -l data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB.recode.vcf | wc -l ")

# Convert to PLINK
system("vcftools --vcf data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB.recode.vcf --plink --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB")
system("wc -l data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB.map")
system("plink --file data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB --recodeAD --make-bed --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB")
system("wc -l data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB.fam")


#~~ HWE

system("plink --file data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB --hwe 0.0001 --nonfounders --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe --recode --make-bed --allow-extra-chr --debug")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Linkage pruning          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# lapply these functions

# recode map file for LD pruning
system("scripts/recode_full_map.sh data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe.map data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe.map")

# Make bed
system("plink --file data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe --make-bed --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe --allow-extra-chr --debug")

# Run PLINK --indep
system("plink --bfile data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe --indep 50 5 2 --nonfounders --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe --allow-extra-chr --debug")
system("wc -l data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe.prune.in")

# LD filter

# Make plink raw file
system("plink --bfile data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe --extract data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe.prune.in --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD --recodeAD --allow-extra-chr --debug")
# Make standard plink files
system("plink --bfile data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe --extract data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe.prune.in --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD --make-bed --recode --allow-extra-chr --debug")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Individual missingness       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

system("plink --bfile data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD --missing --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD_missing --allow-extra-chr --debug")

miss <- fread("data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD_missing.imiss")

png("data/inbreeding/miss.png")
hist(miss$F_MISS)
dev.off()

# filter for individuals with missing data

system("plink --bfile data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD --mind 0.8 --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD_missing --recodeAD --allow-extra-chr --debug")
system("plink --bfile data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD --mind 0.8 --out data/inbreeding/GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD_missing --make-bed --recode --allow-extra-chr --debug")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Inbreeding coefficients     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Get sMLH values

get_sMLH_from_plinkraw <- function(file) {
  
  x <- fread(file, colClasses = "character")
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  row.names(x) <- ids
  NAs <- apply(x, 1, function(x) sum(is.na(x)))
  
  sMLH <- as.data.frame(sMLH(x))
  sMLH$ANIMAL <- ids
  sMLH$NAS <- NAs
  sMLH <- dplyr::filter(sMLH, grepl("^AG", ANIMAL))
  colnames(sMLH) <- c("sMLH", "ANIMAL", "NAs")
  sMLH <- sMLH
  
}


raw_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.raw"), sep = "")
sMLH <- lapply(raw_files[2], get_sMLH_from_plinkraw)
names(sMLH) <- lapply(raw_files[2], function(x) gsub("data/inbreeding/|\\.raw", "", x))

sMLH <- rbindlist(sMLH, idcol = "Run") %>%
  dplyr::rename(IID = ANIMAL)

#~~ Load microsatellite data

ms <- fread("data/raw/Rack61_microsatellites_Jul2017.csv", header = T) %>%
  right_join(filter(sMLH, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD"), by = c("ID" = "IID"))
IDs <- ms$ID
rownames(ms) <- ms$ID

ms <- dplyr::select(ms, -Run, -ID, -sMLH, -NAs)

raw_ms <- convert_raw(ms)
ms_sMLH <- data.frame(sMLH(raw_ms))
ms_sMLH$IID <- IDs
colnames(ms_sMLH) <- c("ms_sMLH", "IID")


#~~ g2 msats

g2_ms <- g2_microsats(raw_ms, nperm = 1000, nboot = 1000, CI = 0.95)
g2_ms

rm(ms)
rm(raw_ms)


#~~ g2 (boot over loci)

get_g2_from_plinkraw <- function(file) {
  
  x <- fread(file, colClasses = "numeric") %>%
    filter(grepl("^AG", IID))
  
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  row.names(x) <- ids
  
  library(inbreedR)
  g2 <- g2_snps(x, nperm = 1000, nboot = 1000, CI = 0.95)
  g2
  g2 <- g2
  
}

g2 <- lapply(raw_files[2], get_g2_from_plinkraw)
g2 <- g2[[1]]
plot(g2, col = "grey")
g2
g2$g2


#~~ Get Fhats

# recode bim files for GCTA

recode_bim_1chr <- function(file){
  file <- fread(file) %>%
    mutate(V1 = 1) %>%
    fwrite(file, quote = F, row.names = F,
           col.names = F, sep = " ")
}

bim_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.bim"), sep = "")
lapply(bim_files, recode_bim_1chr)

# get fhats using gcta

plink_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.ped"), sep = "")
plink_files <- lapply(plink_files, function(x) gsub(".ped", "", x))

for (i in 1:length(plink_files)){
  system(paste0("~/programs/gcta64 --bfile ", plink_files[i]," --autosome --ibc --out ", plink_files[i]," --thread-num 10"))
}

# load fhats

load_fhats <- function(file) {
  
  fhats <- fread(file, header = T) %>%
    mutate(Animal = ifelse(grepl("^AGP", IID), "Pup", "Adult"))
  
}


ibc_files <- paste("data/inbreeding/", list.files(path = "data/inbreeding", pattern="*.ibc"), sep = "")
fhats <- lapply(ibc_files[2], load_fhats)
names(fhats) <- lapply(ibc_files[2], function(x) gsub("data/inbreeding/|\\.ibc", "", x))

ibcs <- rbindlist(fhats, idcol = "Run") %>%
  left_join(sMLH, by = c("IID", "Run")) %>%
  left_join(ms_sMLH, by  = c("IID"))


library(ggthemr)

ggthemr(palette = "pale", layout = "clean", 
        line_weight = 0.7, text_size = 20, type = "outer")
swatch()
to_swap <- swatch()[3:4]


# ~~ Plot g2

g2$g2 / var(ibcs$sMLH, na.rm = T)
g2$g2 / var(ibcs$Fhat1)
g2$g2 / var(ibcs$Fhat2)
g2$g2 / var(ibcs$Fhat3)


g2_plot <- data.frame(g2$g2_boot)
lcl <- g2$CI_boot[1]
ucl <- g2$CI_boot[2]
g2_boot_summary <- data.frame(lcl, ucl)

ibcs_vars_summary <- data.frame(c(g2$g2,
                                  var(ibcs$sMLH, na.rm = T),
                                  var(ibcs$Fhat1),
                                  var(ibcs$Fhat2),
                                  var(ibcs$Fhat3)),
                                c("g2", "sMLH", "Fhat1", "Fhat2", "Fhat3"))

colnames(ibcs_vars_summary) <- c("val", "var")

plot(g2)


# g2 bootstrapping distribution showing empirical g2 with CIs

require(gridExtra)
library(sitools)
cbPalette <- c( "#1B9E77", "#66A61E", "#E6AB02", "black", "#7570B3", "#D95F02", "#E7298A")

png("data/inbreeding/g2_boot.png", units = "in", res = 300, width = 8, height = 7)

#g2_CI_plot <- 
  ggplot(g2_plot, aes(g2$g2_boot)) + 
  geom_histogram(colour = "grey45", fill = "grey45") +
  geom_errorbarh(aes(xmin = g2_boot_summary$lcl , xmax = g2_boot_summary$ucl , y = 90),
                 size = 0.8, color = "black", linetype = "solid", height = 0) +
  geom_vline(data = ibcs_vars_summary, aes(xintercept = val, colour = var), size = 0.8, 
             linetype = c("dashed", "solid", "solid", "solid", "solid"), show.legend = T) +
  scale_colour_manual(values = cbPalette, name = "",
                      breaks = c("g2","Fhat1", "Fhat2", "Fhat3", "sMLH"),
                      labels = c(expression(italic(g[2])), 
                                 expression("var"(italic(hat(F)["I"]))), 
                                 expression("var"(italic(hat(F)["II"]))), 
                                 expression("var"(italic(hat(F)["III"]))), 
                                 expression("var"("sMLH")))) +
  labs(y = "Counts", x = expression(italic(g[2]))) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  ggtitle('(B)') + theme(plot.title=element_text(hjust=0, size = 18, face = "plain"))


dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#               Distribution               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


distribution <- ggplot(filter(ibcs, Run.x == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD"), aes(x=Fhat3, color=.id, fill=.id), col = "grey33") +
  geom_histogram(aes(y=..density..),  position="identity", alpha=0.9, col = "grey33", fill = "grey33") +
  geom_density(alpha=0.6, col = "grey33", fill = "grey33") + #3262AB
  labs(y = "Density", x = expression(italic(hat(F)["III"]))) +
  #xlim(c(-0.5,0.5)) +
  xlim(c(-0.2,0.2)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  ggtitle('(A)') + theme(plot.title=element_text(hjust=0, size = 18, face = "plain"))

library(gridExtra)

g2_dist_panel <- grid.arrange(distribution, g2_CI_plot, ncol = 2, nrow = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Adult Pup Distribution          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Adult pup distribution

#library(RColorBrewer)
#cbPalette <- c(brewer.pal(6, "Dark2"))

png("figs/ibc_pup_adult.png", units = "in", res = 300, width = 11, height = 8)

adult_pup_fhat3 <- ggplot(filter(ibcs, Run.x == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD"), aes(x=Fhat3, fill=Animal, color = Animal)) +
  #geom_histogram(aes(y=..density..),  position="identity", alpha=0.9, col = "#3262AB", fill = "#3262AB") +
  geom_density(alpha=0.8) +
  labs(y = "Density", x = expression(italic(hat(F)["III"]))) +
  xlim(c(-0.2,0.2)) +
  #xlim(c(0.75, 1.3)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  scale_color_manual(values=c("black", "black")) +
  # scale_color_manual(values=c("grey33", "#EB9050")) +
  scale_fill_manual(values=c("grey33", "#EB9050"))

dev.off()



png("figs/ibc_and_g2.png", units = "in", res = 300, width = 13, height = 6)
grid.arrange(adult_pup_fhat3, g2_CI_plot,
             ncol=2, nrow=1)
dev.off()


#~~ Permutation test

pups <- filter(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_LD"), Animal == "Pup")$Fhat3
adults <- filter(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_LD"), Animal == "Adult")$Fhat3
obsdiff <- mean(pups)-mean(adults)

N <- length(pups)

avgdiff <- replicate(10000, {
  all <- sample(c(pups,adults)) # scrambles data
  newpups <- all[1:N] # first half of data
  newadults <- all[(N+1):(2*N)] # second half of data
  return(mean(newpups) - mean(newadults)) # new mean
})

hist(avgdiff)
abline(v=obsdiff, col="red", lwd=2)
(sum(abs(avgdiff) > abs(obsdiff)) + 1) / (length(avgdiff) + 1)

median(pups)
median(adults)



#~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Correlation plots     #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#lm_eqn = function(m) {
#  eq <- substitute(italic(r)^2~"="~r2, 
#                   list(r2 = format(summary(m)$r.squared, digits = 3)))
#  as.character(as.expression(eq));   
#}

corr_eqn <- function(x,y, digits = 2) {
  corr_coef <- round(cor(x, y), digits = digits)
  paste("italic(r) == ", corr_coef)
}


#m1 <- lm(Fhat1 ~ Fhat3, filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD"))
labels1 <- data.frame(x = -0.04, y = 0.34, label = corr_eqn(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD")$Fhat1, 
                                                           filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD")$Fhat3))

plot7 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD"), aes(x=Fhat3, y = Fhat1)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") + #3262AB
  labs(x = expression(italic(hat(F)["III"])), y = expression(italic(hat(F)["I"]))) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  # geom_text(aes(x = -0.03, y = 0.15, label = lm_eqn(m1)), parse = TRUE,
  #           size = 5, fontface = "plain") +
  geom_text(data = labels1, aes(x = x, y = y, label = label), parse = TRUE,
            size = 5, fontface = "plain") +
  ggtitle('(C)') + theme(plot.title=element_text(hjust=0, size = 18, face = "plain"))


labels2 <- data.frame(x = -0.04, y = 0.15, label = corr_eqn(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD")$Fhat2, 
                                                           filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD")$Fhat3))

plot8 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD"), aes(x=Fhat3, y = Fhat2)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") +
  labs(x = expression(italic(hat(F)["III"])), y = expression(italic(hat(F)["II"]))) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain")) +
  geom_text(data = labels2, aes(x = x, y = y, label = label), parse = TRUE,
            size = 5, fontface = "plain") +
  ggtitle('(D)') + theme(plot.title=element_text(hjust=0, size = 18, face = "plain"))


scaleFUN <- function(x) sprintf("%.1f", x)
labels3 <- data.frame(x = 0.12, y = 1.15, label = corr_eqn(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD")$sMLH, 
                                                          filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD")$Fhat3))
plot9 <- ggplot(filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD"), aes(x=Fhat3, y = sMLH)) +
  geom_point(size = 2, alpha = 1/1.5, col = "grey33") +
  labs(x = expression(italic(hat(F)["III"])), y = "sMLH") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  geom_text(data = labels3, aes(x = x, y = y, label = label), parse = TRUE,
            size = 5, fontface = "plain") +
  ggtitle('(E)') + theme(plot.title=element_text(hjust=0, size = 18, face = "plain")) +
  scale_y_continuous(labels=scaleFUN)


library(gridExtra)


png("figs/ibc_correlations.png", units = "in", res = 300, width = 10, height = 9) # 22, 6
grid.arrange(plot4, plot7,
             plot8, plot9, ncol=2, nrow=2)
dev.off()

corr_plots <- grid.arrange(plot7, plot8, plot9, ncol = 3, nrow = 1)



#~~ Plot for manuscript:

png("figs/inbreeding_fig.png", units = "in", res = 300, width = 14, height = 10) # 22, 6
grid.arrange(g2_dist_panel, corr_plots, nrow = 2)
dev.off()

jpeg("figs/inbreeding_fig.jpg", units = "in", res = 300, width = 14, height = 10) # 22, 6
grid.arrange(g2_dist_panel, corr_plots, nrow = 2)
dev.off()

#~~ ms RAD correlation

png("data/inbreeding/ms_vs_rad.png")
ggplot(ibcs, aes(x = ms_sMLH, y = sMLH)) + 
  geom_point(colour = "grey45") +
  labs(y = "SNP sMLH", x = "Microsatellite sMLH" ) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        axis.ticks = element_blank()) 
dev.off()

#~~ Heritability of heterozygosity

ped <- read.table("data/raw/new_pedigree.txt", header = T) %>%
  na.omit()
head(ped)

head(ibcs)

h2_ID <- left_join(ped, filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD"), by = c("ANIMAL" = "IID"))
h2_MOTHER <- left_join(ped, filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD"), by = c("MOTHER" = "IID"))
h2_FATHER <- left_join(ped, filter(ibcs, Run == "GQ_DP_maf05_G_ldi_ME_SSB_hwe_LD"), by = c("FATHER" = "IID"))
h2_midParent <- (h2_MOTHER$sMLH + h2_FATHER$sMLH) / 2

plot(h2_ID$sMLH ~ h2_MOTHER$sMLH) +
  abline(lm(h2_ID$sMLH ~ h2_MOTHER$sMLH))

summary(lm(h2_ID$sMLH ~ h2_MOTHER$sMLH))

plot(h2_ID$sMLH ~ h2_FATHER$sMLH) +
  abline(lm(h2_ID$sMLH ~ h2_FATHER$sMLH))

summary(lm(h2_ID$sMLH ~ h2_FATHER$sMLH))


plot(h2_ID$sMLH ~ h2_midParent) +
  abline(lm(h2_ID$sMLH ~ h2_midParent))

lm.het <-  lm(h2_ID$sMLH ~ h2_midParent)
library(ggfortify)
autoplot(lm.het)
lm.het
summary(lm.het)
summary.lm.het <- summary(lm.het)$coefficients 

