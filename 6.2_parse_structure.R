# Parse and plot STRUCTURE output using mainly pophelper()

library(dplyr)
library(tidyr)
library(stringr)
#install.packages('devtools',dependencies=T)
#library(devtools)
#install_github("royfrancis/pophelper", force = T)
library(pophelper)
library(data.table)
options(scipen=999)
library(ggplot2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Collect output files                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Full

system("mkdir data/structure/results/run_files")
system("mv data/structure/results/*_f data/structure/results/run_files/")

system("mkdir data/structure/ms/results/run_files")
system("mv data/structure/ms/results/*_f data/structure/ms/results/run_files/")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Load files and collect clumpp output            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load_and_clumpp <- function(path_to_structure_out){
  
  all_files <- list.files(path_to_structure_out, pattern = "^results")
  
  # creates character vector containing path to relevant files
  struc_out_paths <- paste0(path_to_structure_out, all_files)
  
  slist <- readQ(files=struc_out_paths, filetype = "structure")
  
  # export for CLUMPP later on
  clumppExport(qlist=slist, useexe=T)
  # collect CLUMPP output
  collectClumppOutput(filetype="merged") # aligned
  system("rm -r pop_K*")
  
  # move clump files to correct directory
  system(paste0("mv pop-* ", path_to_structure_out))
  
}


load_and_clumpp("data/structure/results/run_files/")
load_and_clumpp("data/structure/ms/results/run_files/")



# Function to get K summary stats from run files

load_and_K <- function(path_to_structure_out){
  
  all_files <- list.files(path_to_structure_out, pattern = "^results")
  
  # creates character vector containing path to relevant files
  struc_out_paths <- paste0(path_to_structure_out, all_files)
  
  slist <- readQ(files=struc_out_paths, filetype = "structure")
  
  # Get summary stats
  
  em <- evannoMethodStructure(summariseQ(tabulateQ(slist)))
  plot(em$deltaK)
  
  evannoMethodStructure(data=em, exportplot=T, writetable = T)
  system(paste0("mv evannoMethodStructure* ", path_to_structure_out))
  
  em <- em
  
}


rad_ks <- load_and_K("data/structure/results/run_files/")
ms_ks <- load_and_K("data/structure/ms/results/run_files/")


#~~ optimal K 

# if mean estimated ln probability of data is highest at K=1, K=1, 
# if not, use evanno-method (delta K) to choose k

optimal_k <- function(x) {
  elpd <- which.max(x$elpdmean)
  if (elpd == 1) {
    return(1)
  } else {
    delta_k <- which.max(x$deltaK)
    delta_k
  }
}


optimal_k(rad_ks)
optimal_k(ms_ks)



# output best k's

both_k <- function(x) {
  elpd <- which.max(x$elpdmean)
  deltaK <- which.max(x$deltaK)
  ks <- c(lnk = x$k[elpd], deltak = x$k[deltaK])
  # write outfile
  #write.table(ks, paste0(path_to_structure_out, "Ks.txt"))
  ks <- ks
}


rad_bothKs <- both_k(rad_ks)
ms_bothKs <- both_k(ms_ks)


#~~ Make deltak and elpd plot for Supp

require(gridExtra)
library(ggthemr)
ggthemr(palette = "pale", layout = "clean", 
        line_weight = 0.5, text_size = 16, type = "outer", spacing = 2)

#~~ RAD Ks

rad_dk <- ggplot(rad_ks, aes(x=k, y = deltaK)) +
  geom_point(size = 1, col = "grey30") +
  geom_line(size = 1, col = "grey30") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  labs(x = "K", y = expression(paste(Delta,italic("K")))) +
  scale_x_continuous(breaks=c(1:6), labels=c(1:6),limits=c(1,6)) +
  ggtitle('(D) 27,592 SNPs') + theme(plot.title=element_text(hjust=0, size = 16, face = "plain"))

rad_elpd <- ggplot(rad_ks, aes(x=k, y = elpdmean)) +
  geom_point(size = 1, col = "grey30") + # 1/1.5
  geom_line(size = 1, col = "grey30") +
  geom_errorbar(aes(ymin = elpdmean - elpdsd, ymax= elpdmean + elpdsd), colour="grey30", width=0) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  labs(x = "K", y = expression(paste("Ln Pr(",italic("X"),"|",italic("K"),")"))) +
  scale_x_continuous(breaks=c(1:6), labels=c(1:6),limits=c(1,6)) +
  ggtitle('(B) 27,592 SNPs') + theme(plot.title=element_text(hjust=0, size = 16, face = "plain"))



#~~ Microsatellite Ks

ms_dk <- ggplot(ms_ks, aes(x=k, y = deltaK)) +
  geom_point(size = 1, col = "grey30") +
  geom_line(size = 1, col = "grey30") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  labs(x = "K", y = expression(paste(Delta,italic("K")))) +
  scale_x_continuous(breaks=c(1:6), labels=c(1:6),limits=c(1,6)) +
  ggtitle('(C) 27 microsatellites') + theme(plot.title=element_text(hjust=0, size = 16, face = "plain"))

ms_elpd <- ggplot(ms_ks, aes(x=k, y = elpdmean)) +
  geom_point(size = 1, col = "grey30") + # 1/1.5
  geom_line(size = 1, col = "grey30") +
  geom_errorbar(aes(ymin = elpdmean - elpdsd, ymax= elpdmean + elpdsd), colour="grey30", width=0) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  labs(x = "K", y = expression(paste("Ln Pr(",italic("X"),"|",italic("K"),")"))) +
  scale_x_continuous(breaks=c(1:6), labels=c(1:6),limits=c(1,6)) +
  ggtitle('(A) 27 microsatellites') + theme(plot.title=element_text(hjust=0, size = 16, face = "plain"))


png("figs/deltaK_elpd.png", units = "in", res = 300, width = 10, height = 9)
grid.arrange(ms_elpd, rad_elpd, ms_dk, rad_dk, ncol=2, nrow=2)
dev.off()

jpeg("figs/deltaK_elpd.jpg", units = "in", res = 300, width = 10, height = 9)
grid.arrange(ms_elpd, rad_elpd, ms_dk, rad_dk, ncol=2, nrow=2)
dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Structure Plots               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Structure runs not merged:

make_structure_plots <- function(path_to_structure_out, path_to_struc_file){
  
  # Get slist
  all_files <- list.files(path_to_structure_out, pattern = "^results")
  
  # creates character vector containing path to relevant files
  struc_out_paths <- paste0(path_to_structure_out, all_files)
  
  slist <- readQ(files=struc_out_paths, filetype = "structure")
  
  
  #~~ Customize plots
  
  # strip panel label showing k only
  fn1 <- function(x) attr(x,"k")
  spnames <- paste0("K=",sapply(slist,fn1))
  
  # custom colours
  library(ggthemr)
  ggthemr(palette = "solarized", layout = "clean",
          line_weight = 0.7, text_size = 20, type = "outer")
  swatch()
  #[1] "#073642" "#268bd2" "#dc322f" "#2aa198" "#b58900" "#6c71c4" "#d33682"
  
  # add group labels
  struc_file <- path_to_struc_file
  
  pops <- fread(struc_file) %>%
    select(V2) %>%
    mutate(V2 = ifelse(V2 == 1, "SSB",
                       ifelse(V2 == 2, "CS",
                              ifelse(V2 == 3, "B",
                                     ifelse(V2 == 4, "K", 
                                            ifelse(V2 == 5, "HI", "Mac"))))))
  
  #~~ Write out plots
  
  #correctly align plot panel and label panel
  
  plotQ(qlist=(slist)[1:5], splab=spnames[1:5], imgoutput = "join", grplab=pops,
        clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"K1_reps"))
  
  plotQ(qlist=(slist)[6:10], splab=spnames[6:10], imgoutput = "join", grplab=pops,
        clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"K2_reps"))
  
  plotQ(qlist=(slist)[11:15], splab=spnames[11:15], imgoutput = "join", grplab=pops,
        clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"K3_reps"))
  
  plotQ(qlist=(slist)[16:20], splab=spnames[16:20], imgoutput = "join", grplab=pops,
        clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"K4_reps"))
  
  plotQ(qlist=(slist)[21:25], splab=spnames[21:25], imgoutput = "join", grplab=pops,
        clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"K5_reps"))
  
  plotQ(qlist=(slist)[26:30], splab=spnames[26:30], imgoutput = "join", grplab=pops,
        clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"K6_reps"))
  
}



make_structure_plots("data/structure/results/run_files/", "data/structure/handmade_rad_sort.stru")
make_structure_plots("data/structure/ms/results/run_files/", "data/structure/ms/stru_in.stru")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. CLUMPP Structure Plots        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Plot with merged files

make_structure_clumpp_plots <- function(path_to_structure_out, path_to_struc_file){
  
  # Get slist
  all_files <- list.files(path_to_structure_out, pattern = "^results")
  
  # creates character vector containing path to relevant files
  struc_out_paths <- paste0(path_to_structure_out, all_files)
  
  slist <- readQ(files=struc_out_paths, filetype = "structure")
  # File paths
  clumpp_files <- list.files(paste0(path_to_structure_out,"pop-merged/"))
  clumpp_out_paths <- paste0(path_to_structure_out,"pop-merged/", clumpp_files)
  
  # Read in files
  clist <- readQ(files=clumpp_out_paths, filetype = "clumpp")
  
  # Plot aligned structure plots
  
  # strip panel label showing k only
  fn1 <- function(x) attr(x,"k")
  spnames <- paste0("K=",sapply(slist,fn1))
  
  struc_file <- path_to_struc_file
  
  pops <- fread(struc_file) %>%
    select(V2) %>%
    mutate(V2 = ifelse(V2 == 1, "SSB",
                       ifelse(V2 == 2, "CS",
                              ifelse(V2 == 3, "B",
                                     ifelse(V2 == 4, "K", 
                                            ifelse(V2 == 5, "HI", "Mac"))))))
  
  plotQ(qlist=(clist)[1:5], splab=spnames[c(6,11,16,22,26)], imgoutput = "join", grplab=pops,
        clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"K1-K6_merged"))
  
  # plot with individuals sorted by population assignment
  # plotQ(clist, sortind="all")
  
}


make_structure_clumpp_plots("data/structure/results/run_files/", "data/structure/handmade_rad_sort.stru")
make_structure_clumpp_plots("data/structure/ms/results/run_files/", "data/structure/ms/stru_in.stru")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  8. Membership Probabilites      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rad_k2_prob <- read.table("data/structure/results/run_files/pop-merged/pop_K2-combined-merged.txt") %>%
  mutate(top = ifelse(V2 > V3, V2, V3))
nrow(filter(rad_k2_prob, top >= 0.99)) / 37

ms_k2_prob <- read.table("data/structure/ms/results/run_files/pop-merged/pop_K2-combined-merged.txt")  %>%
  mutate(top = ifelse(V2 > V3, V2, V3))
nrow(filter(ms_k2_prob, top >= 0.99)) / 37

reb_k2_prob <- read.table("data/rebuttal/results/run_files/pop-merged/pop_K4-combined-merged.txt") %>%
  mutate(top = ifelse(V2 > V3, V2, V3))
nrow(filter(reb_k2_prob, top >= 0.99)) / 37


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  9. DISTRUCT Plots               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

cols <- c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E") # colours wrong

# select first run of K=2
distructExport(slist[6], grplabbottom=pops$V2, useexe = T, indwidth=1.2, printcolorbrewer = T, fontsize = 3, xscale = 5, yscale = 3,
               clustercol = cols)


system(paste0("mkdir ", path_to_structure_out, "distruct_plots"))
system(paste0("mv results_job** ", path_to_structure_out, "distruct_plots"))


