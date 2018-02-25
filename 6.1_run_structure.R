# Run this Rscript to run structure in parallel

library(ParallelStructure)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run ParallelStructure                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Usage: nohup Rscript 7_run_structure.R &

options(scipen=999)

# Structure files

gen <- fread("data/structure/handmade_rad_sort.stru", header = F)
#gen <- fread("data/structure/ms/stru_in.stru", header = F)


#~~ Specify in and out files for structure

infile <- "data/structure/handmade_rad_sort.stru"
outpath <- "data/structure/results/"

#infile <- "data/structure/ms/stru_in.stru"
#outpath <- "data/structure/ms/results/"

#~~ construct job matrix and write to job file

nrep <- 5
#burnin <- 50000
#niter <- 150000
up_to_k <- 6

niter <- 1000000
burnin <- 100000

#~~ define variables for job matrix

k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

#~~ make the job matrix
pop <- "1,2,3,4,5,6"

seal_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(seal_jobs), ncol = length(seal_jobs[1,]), file = "seal_jobs.txt")

#~~ file path to structure

STR_path='/usr/local/bin/'


#~~ Run Parallel Structure

system("mkdir data/structure/results")
system("mkdir data/structure/ms/results")


#~~ Run structure

setwd("~/rad_analysis_pilon_Dec2017")
ParallelStructure::parallel_structure(structure_path=STR_path, joblist="seal_jobs.txt", n_cpu=20, 
                                      infile=infile, 
                                      outpath= outpath, 
                                      numinds = nrow(gen),
                                      numloci=(ncol(gen)-2)/2, 
                                      #popdata = 1, 
                                      #popflag = 1, 
                                      #usepopinfo = 1, 
                                      printqhat=1, 
                                      plot_output=0, 
                                      onerowperind=1)

