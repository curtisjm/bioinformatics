#!/usr/bin/env Rscript

#SBATCH --job-name=whole_chromosome

#SBATCH --account=fc_williamslab

#SBATCH --partition=savio2_bigmem

#SBATCH --time=40:00:00


library(ggplot2)
library(dplyr)
library(reshape2)


#make sure to setwd()

#path to folder with all of the files - NO DUPLICATES
path <- "/global/scratch/users/nksmoot/emseq/whole_xsome"


#name of sample needs to be the same as the file name
#IN ORDER OF AGE
samples <- c('4-2_R1', '4-4_R1')

#assign colors to each sample in this format
#i think you need to keep it in the same order
colors = c("4-2_R1" = "#840729", 
           "4-4_R1" = "#b5152c")

chrom <- paste0("Chr", 1:5)

#window sizes
win_size <- 250000
small_win_size = 1000

#pericentromeric region endpoints
#NEED TO NARROW THIS DOWN MORE(?)
centromeres <- list('Chr1' = list('min' = 15000000, 'max' = 17000000),
                    'Chr2' = list('min' = 4000000, 'max' = 6000000),
                    'Chr3' = list('min' = 14000000, 'max' = 16000000),
                    'Chr4' = list('min' = 4000000, 'max' = 6000000),
                    'Chr5' = list('min' = 12000000, 'max' = 14000000))

#cutoff for filtering out low coverage windows
cutoff=50

samples_sort <- sort(samples)

#import all CpG data for each chromosome in each sample
#file names need to be in the format "samplename_all_CpG_Chr#.bed"
import_data <- function(path, tissue, chromosome) {
  bed_paths <- list.files(path, recursive = T, 
                          pattern = paste0("^(", paste(tissue, collapse = "|"), 
                                           ")", "_all_CpG_", chromosome, ".bed"), 
                          full.names = T)
  data <- lapply(bed_paths, read.table, header = F, 
                 col.names = c('chr', 'start', 'stop', 'unmeth', 'meth', 'level'))
  names(data) <- tissue
  return(data)
}


all_data <- lapply(chrom, import_data, path = path, tissue = samples_sort)

names(all_data) <- chrom



#####MAX STUFF MAY BE UNNECESSARY
#find max stop site for a chromosome
max_stop <- function(chrom, win_size) {
  #helper, max stop site for a sample at a specific chrom
  max_helper <- function(tissue) {
    return(max(tissue$stop))
  }
  max_by_tissue <- lapply(chrom, max_helper)
  max_by_tissue <- max(unlist(max_by_tissue))
  #round up to nearest win_size multiple
  max_by_tissue <- ceiling(max_by_tissue/win_size) * win_size 
  return(max_by_tissue)
}

#list of chromosome endpoints
max_list <- lapply(all_data, max_stop, win_size = win_size)

#mean for each window
split_windows <- function(chrom_data, max, win_size) {
  split_helper <- function(tissue, max, win_size) {
    id <- cut(tissue$start, seq(0, max, win_size))
    win_data <- as.data.frame(tapply(tissue$level, id, mean))
    names(win_data) <- "mean"
    win_data$totals <- tapply(tissue$unmeth + tissue$meth, id, sum)
    win_data$bins <- rownames(win_data)
    win_data$index <- 1:nrow(win_data)
    return(win_data)
  }
  chrom_wins <- lapply(chrom_data, split_helper, 
                       max=max, win_size=win_size)
  return(chrom_wins)
}

all_wins <- mapply(split_windows, chrom_data=all_data, max=max_list,
                   MoreArgs = list(win_size=win_size),
                   SIMPLIFY = F)
print("all_wins:")
head(all_wins)


#plotting

make_lines <- function(chrom_data, samples, colors, chrom, win_size) {
  line_helper <- function(tissue, sample) {
    return(geom_line(data=tissue, aes(x=index, y=mean, color=sample)))
  }
  
  geoms <- mapply(line_helper, tissue=chrom_data, sample=samples, SIMPLIFY=F)
  
  plt <- ggplot() + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25), expand = 
expansion(mult = c(0, 0)))+ 
    scale_x_continuous(breaks = seq(from = 0, to = 120, by = 40),expand = expansion(mult = c(0, 0))) + 
    
    ylab("Methylation Level") + xlab(paste0("Chromosome Location (", 
                                            as.character(win_size/1000000), " Mbp)")) +
    ggtitle(paste0(chrom, " Methylation Level")) +
    theme_classic() 


  
  for (i in names(colors)) {
    plt <- plt + geoms[[i]]
  }
  return(plt+ scale_color_manual(values = colors))
}
all_gg <- mapply(make_lines, chrom_data=all_wins, chrom=chrom,
                 MoreArgs=list(samples=samples_sort, colors=colors,
                               win_size=win_size),
                 SIMPLIFY = F)


all_gg

