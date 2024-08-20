.libPaths( c('/global/scratch/users/kenchen/R', .libPaths()) )

library(ggplot2)
library(dplyr)
library(reshape2)

#make sure to setwd()

#path to folder with all of the files - NO DUPLICATES
path <- "C:/Dawei_UCB/1_UCB_Research/Aging/bed_files/grep/CHG"


#name of sample needs to be the same as the file name
#IN ORDER OF AGE
samples <- c('sperm', 'embryo',
             'D9_Col0', 'D13_Col0', 'D18_Col0', 'D23_Col0', 'D30_Col0', 'D45_Col0', 'D60_Col0')

#assign colors to each sample in this format
#i think you need to keep it in the same order
colors = c("sperm" = "black",
           "embryo" = "#45102d", "D9_Col0" = "#840729", 
           "D13_Col0" = "#b5152c", "D18_Col0" = "#e33834",
           "D23_Col0" = "#fe7361", "D30_Col0" = "#f99c7c", "D45_Col0" = "#fcc7a9",
           "D60_Col0" = "#f8e3d3")

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
                                           ")", "_all_CHG_", chromosome, ".bed"), 
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

#plotting

make_lines <- function(chrom_data, samples, colors, chrom, win_size) {
  line_helper <- function(tissue, sample) {
    return(geom_line(data=tissue, aes(x=index, y=mean, color=sample)))
  }
  
  geoms <- mapply(line_helper, tissue=chrom_data, sample=samples, SIMPLIFY=F)
  
  plt <- ggplot() + scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 25), expand = expansion(mult = c(0, 0)))+ 
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
























#extract centromere region data
get_region <- function(chrom_data, centromere, win_size) {
  region_helper <- function(tissue, centromere, win_size) {
    tissue <- subset(tissue, index*win_size >= centromere$min &
                       index*win_size <= centromere$max)
    return(tissue)
  }
  chrom_data <- lapply(chrom_data, region_helper, centromere, win_size)
}

filt_wins <- mapply(get_region, chrom_data = all_wins, 
                    centromere = centromeres, 
                    win_size = win_size, SIMPLIFY = F)
filt_gg <- mapply(make_lines, chrom_data=filt_wins, chrom=chrom,
                  MoreArgs=list(samples=samples_sort, colors=colors,
                                win_size=win_size),
                  SIMPLIFY = F)


###############
#IF NOT INTERESTED IN LOOKING AT SMALLER WINDOWS, SKIP TO INSPECTING PLOTS AND DATA SECTION
#############

#####
#by small windows
#####

small_max_list <- lapply(all_data, max_stop, win_size = small_win_size)

#optimize this section by getting region of interest before splitting
small_wins <- mapply(split_windows, chrom_data=all_data, max=small_max_list,
                     MoreArgs = list(win_size=small_win_size),
                     SIMPLIFY = F)
filt_small_wins <- mapply(get_region, chrom_data = small_wins, 
                          centromere = centromeres, 
                          win_size = small_win_size, SIMPLIFY = F)
small_gg <- mapply(make_lines, chrom_data=filt_small_wins, chrom=chrom,
                   MoreArgs=list(samples=samples_sort, colors=colors, 
                                 win_size = small_win_size),
                   SIMPLIFY = F)

############
#merge + filter data from each chromosome for formatting for ggplot

merge_samples <- function(chrom_data, cutoff) {
  cutoff_filt <- function(chrom_data, cutoff) {
    for (col in seq(2, length(chrom_data), 2)) {
      for (row in 1:dim(chrom_data)[1]) {
        if (!is.na(chrom_data[[row, col]]) & chrom_data[[row,col]] < cutoff) {
          chrom_data[[row, col-1]] = NA
        }
      }
    }
    return(chrom_data)
  }
  data <- lapply(chrom_data, cutoff_filt, cutoff=cutoff)
  data <- lapply(data, select, "mean")
  data <- do.call(cbind, data)
  names(data) <- c(names(chrom_data))
  return(data)
}

small_merged <- lapply(filt_small_wins, merge_samples, cutoff=cutoff)
###########

plot_windows <- function(merged_data, win_size, samples) {
  plot_helper <- function(merged_chrom, chrom) {
    df <- as.data.frame(t(merged_chrom))
    df$tissue <- rownames(df)
    df$tissue <- factor(df$tissue, 
                        levels=samples)
    df <- melt(df, id.vars=c("tissue"))
    df <- na.omit(df)
    plot <- ggplot(df, aes(x=tissue, y=value, group=variable)) + 
      geom_line(alpha = 0.03)+#win_size/10^5) + 
      ggtitle(paste0(chrom, " Window Methylation (", win_size, " bp)")) +
      ylab("Methylation Level") +
      xlab("Sample") +
      ylim(40, 100)
    plot
    return(plot)
  }
  plots <- mapply(plot_helper, merged_data, names(merged_data), 
                  SIMPLIFY = F)
  return(plots)
}

small_win_plots <- plot_windows(small_merged, small_win_size, samples)


################
#Inspecting plots and data
############

#view plots with large windows
all_gg$Chr1 #change to Chr2 to view Chr2 and so on

#view where pericentromeric region cutoffs occur
#change Chr in all 3 places to view different ones
all_gg$Chr1 + geom_vline(xintercept = centromeres$Chr1$min/win_size) + 
  geom_vline(xintercept = centromeres$Chr1$max/win_size)


#view large windows, but only pericentromeric regions
filt_gg$Chr1


#view small window plots
#this will probably look very ugly. that's expected
small_gg$Chr1

#view how methylation in each small window changes over time
small_win_plots$Chr1



