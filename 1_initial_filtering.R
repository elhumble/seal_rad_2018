# Initial filtering of RAD genotypes

library(data.table)
library(dplyr)
library(seqinr)
library(tidyr)
library(plyr)
options(scipen=999)
library(ggplot2)
library(readxl)
library(stringr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Load raw vcf file                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

getwd()

# Count number of SNPs in vcf file

system("grep -v '#' /external/data/Emily/Seal_RADseq/cebitec_2016/AFS/SAMPLES_PILON/GATK/ArcGaz_pilon_genotype_gvcf.vcf | wc -l")

# Update vcf to include only biallelic snps

system("bcftools view -m2 -M2 -v snps /external/data/Emily/Seal_RADseq/cebitec_2016/AFS/SAMPLES_PILON/GATK/ArcGaz_pilon_genotype_gvcf.vcf > data/ArcGaz_pilon_biallelic.vcf")
system("grep -v '#' data/ArcGaz_pilon_biallelic.vcf | wc -l") # 869941

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Recode sample names                                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# current sample names

system("bcftools query -l data/ArcGaz_pilon_biallelic.vcf")

# use list of new names in correct order to recode

system("bcftools reheader -s data/raw/sample_names.txt -o data/ArcGaz_pilon_biallelic_rename.vcf data/ArcGaz_pilon_biallelic.vcf")

# new sample names

system("bcftools query -l data/ArcGaz_pilon_biallelic_rename.vcf")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. VCFtools filtering                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# filter for individuals in correct triads

system("vcftools --vcf data/ArcGaz_pilon_biallelic_rename.vcf --out data/ArcGaz_pilon_biallelic_sub --keep data/raw/individuals_to_keep.txt --recode --recode-INFO-all")
system("bcftools query -l data/ArcGaz_pilon_biallelic_sub.recode.vcf")

# remove SNPs with 0 or > 0.5 maf

system("vcftools --vcf data/ArcGaz_pilon_biallelic_sub.recode.vcf --out data/ArcGaz_pilon_biallelic_sub_maf --maf 0.0000001 --max-maf 0.5 --recode --recode-INFO-all")
system("grep -v '#' data/ArcGaz_pilon_biallelic_sub_maf.recode.vcf | wc -l") # 677607

# get maf stats

system("vcftools --vcf data/ArcGaz_pilon_biallelic_sub_maf.recode.vcf --freq2 --out data/ArcGaz_pilon_biallelic_sub_maf")


# get depth stats
system("vcftools --vcf data/ArcGaz_pilon_biallelic_sub_maf.recode.vcf --geno-depth --out data/ArcGaz_pilon_biallelic_sub_maf")

depth <- fread("data/ArcGaz_pilon_biallelic_sub_maf.gdepth", colClasses = "numeric") %>%
  unite(SNP, CHROM, POS)

depth <- dplyr::select(depth, -SNP)

depth <- depth %>%
  dplyr::mutate(mean = rowMeans(.), 
                sum = rowSums(.), 
                logmeandepth = log10(mean))

CI <- 0.90
CI_meanDP <- stats::quantile(depth$mean, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_sumDP <- stats::quantile(depth$sum, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_meanDP
CI_sumDP

hist(depth$logmeandepth)
hist(log10(depth$sum))
hist(depth$sum)

good_depth <- filter(depth, sum > 400)


#~~ Average DOC per SNP

mean(depth$sum)
summary(depth$sum)

#~~ Mean depth per individual stats

system("vcftools --vcf data/ArcGaz_pilon_biallelic_sub_maf.recode.vcf --depth --out data/ArcGaz_pilon_biallelic_sub_maf")

idepth <- fread("data/ArcGaz_pilon_biallelic_sub_maf.idepth", colClasses = "numeric")
mean(idepth$MEAN_DEPTH)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Get Mendel error SNPs        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# make plink file

system("vcftools --vcf data/ArcGaz_pilon_biallelic_sub_maf.recode.vcf --plink --out data/ArcGaz_pilon_biallelic_sub_maf")

# recode map file to include contig in chr column

system("scripts/recode_full_map.sh data/ArcGaz_pilon_biallelic_sub_maf.map data/ArcGaz_pilon_biallelic_sub_maf.map")

# Make bed

system("plink --file data/ArcGaz_pilon_biallelic_sub_maf --make-bed --out data/ArcGaz_pilon_biallelic_sub_maf --allow-extra-chr --debug")

# Determine Working Pedigree for PLINK
# Pedigree code from Susan Johnston: https://goo.gl/iRw3Ro 

pedigree <- read.table("data/raw/new_pedigree.txt", header = T, 
                       colClasses = c("character", "character", "character"))
pedigree$MOTHER <- as.character(pedigree$MOTHER)
pedigree$FATHER <- as.character(pedigree$FATHER)

famped <- NULL

for(i in pedigree$ANIMAL){
  ped1 <- pedigree[which(pedigree$ANIMAL == i),]
  {
    ped2 <- ped1[,1:3]
    ped2 <- rbind(data.frame(ANIMAL = c(unlist(ped2[1, 2:3])), FATHER = 0, MOTHER = 0), ped2)
    famped <- rbind(famped, ped2)
    famped <- famped[which(famped[,1] != 0),]
    rm(ped2)
  }
}

fatherID <- famped$FATHER
fatherID <- fatherID[fatherID != 0]
fatherID <- na.omit(fatherID)
fatherID <- rep(fatherID, each = 3)
popID <- as.character(famped[c(73:98),1])

famIDs <- c(fatherID, popID)
famped$Family <- famIDs

famped <- famped %>%
  dplyr::mutate(Phenotype = 0,
                Sex = c(rep(c(2,1,"unknown"),24), rep("unknown", 26))) %>%
  dplyr::mutate(FATHER = ifelse(is.na(FATHER),0,FATHER),
                MOTHER = ifelse(is.na(MOTHER),0,MOTHER))

famped <- famped %>%
  dplyr::distinct(ANIMAL, .keep_all = T)

# Write over fam file
# Arrange rows to match plink
# COLS: FAMID, IID, FATHER, MOTHER, SEX, PHENOTYPE

head(famped)
fam <- fread("data/ArcGaz_pilon_biallelic_sub_maf.fam")


new_famped <- fam %>%
  left_join(famped, by = c("V2" = "ANIMAL"))
head(new_famped)

write.table(new_famped[c(9,1,7,8,11,10)], "data/ArcGaz_pilon_biallelic_sub_maf.fam", col.names = F,
            row.names = F, quote = F)




#~~ Run PLINK --mendel

system("mkdir data/mendel")
system("plink --bfile data/ArcGaz_pilon_biallelic_sub_maf --mendel --out data/mendel/ArcGaz_pilon_biallelic_sub_maf --allow-extra-chr --debug")

# Process Mendelian inconsistencies

# Load mendel output and NA errors due to missing parental genotypes

menderr_plink <- fread("data/mendel/ArcGaz_pilon_biallelic_sub_maf.mendel") %>%
  mutate(V6 = gsub("\\*", NA, V6), V8 = gsub("\\*", NA, V8)) %>%
  mutate(V6 = gsub("NA\\/NA", NA, V6), V8 = gsub("NA\\/NA", NA, V8)) %>%
  na.omit()

# Determine per SNP error rate

error_rate <- menderr_plink %>%
  group_by(V4) %>%
  dplyr::summarise(length(V4)) %>%
  `colnames<-`(c("V4", "N")) %>%
  filter(N > 0) %>%
  mutate(rate = (N/24)*100) # 24 N of triads

hist(error_rate$rate)

# write out list of Mendel Error SNPs

# strict
mendel_error_strict <- error_rate %>%
  filter(N > 1) %>% # wrong in one or more triads (4%)
  dplyr::select(V4) %>%
  separate(V4, c("Contig", "Position"), sep = ":") 

write.table(mendel_error_strict, "data/mendel/MendelErrorStrict_SNPs.txt", col.names = F,
            row.names = F, quote = F, sep = ":")
write.table(mendel_error_strict, "data/mendel/MendelErrorStrict_SNPs_tab.txt", col.names = F,
            row.names = F, quote = F, sep = "\t")

mendel_error_mod <- error_rate %>%
  filter(N > 4) %>% # wrong in more than 4 (16%)
  dplyr::select(V4) %>%
  separate(V4, c("Contig", "Position"), sep = ":") 

write.table(mendel_error_mod, "data/mendel/MendelErrorMod_SNPs.txt", col.names = F,
            row.names = F, quote = F, sep = ":")
write.table(mendel_error_mod, "data/mendel/MendelErrorMod_SNPs_tab.txt", col.names = F,
            row.names = F, quote = F, sep = "/t")

system("wc -l data/mendel/MendelErrorMod_SNPs.txt") # 617
system("wc -l data/mendel/MendelErrorStrict_SNPs.txt") # 26202


