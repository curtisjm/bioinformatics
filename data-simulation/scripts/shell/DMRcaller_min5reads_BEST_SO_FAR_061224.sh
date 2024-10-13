#!/usr/bin/env Rscript

#SBATCH --job-name=DMRcaller

#SBATCH --account=fc_williamslab

#SBATCH --partition=savio3

#SBATCH --time=01:00:00

#SBATCH --output=/global/scratch/users/nksmoot/emseq/dmr_caller_updated/rug/slurm-%j.out

library("DMRcaller")
library("GenomicRanges")
library("rtracklayer")

#Code we need to accept user-specified parameters 
#(longer and more clunky to make it easier for us
#to input variables with -1, -2, -n, -o)

args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

# Initialize variables
file1 <- ""
file2 <- ""
outputName <- ""
outputDir <- ""

# Function to process arguments
processArgs <- function(args) {
  for (i in seq_along(args)) {
    if (args[i] == '-1') {
      file1 <<- args[i + 1]
    } else if (args[i] == '-2') {
      file2 <<- args[i + 1]
    } else if (args[i] == '-n') {
      outputName <<- args[i + 1]
    } else if (args[i] == '-o') {
      outputDir <<- args[i + 1]
    }
  }
}

# Call the function to process arguments
processArgs(args)

# Check if all required arguments are provided
if (file1 == "" || file2 == "" || outputName == "" || outputDir == "") {
  stop("Error: Missing required arguments. Please provide -1, -2, -n, and -o options.")
}



#------------------------------------------------------------------------------------------------------------------------
#CALCULATING CONVERSION RATE

#Load the cytosine reports:
CX_report_file1 <- DMRcaller::readBismark(file1)
CX_report_file2 <- DMRcaller::readBismark(file2)

#Extract the chloroplast methylation data:
PtDNA_file1 <- CX_report_file1[seqnames(CX_report_file1) == "ChrC"]
PtDNA_file2 <- CX_report_file2[seqnames(CX_report_file2) == "ChrC"]

#Calculate conversion:
conversion_file1 <- 1 - (sum(mcols(PtDNA_file1)$readsM) / sum(mcols(PtDNA_file1)$readsN))
conversion_file2 <- 1 - (sum(mcols(PtDNA_file2)$readsM) / sum(mcols(PtDNA_file2)$readsN))



#------------------------------------------------------------------------------------------------------------------------

#CORRECTING FOR CONVERSION RATE
#First, correcting number of methylated reads based on conversion rate
CX_report_file1_adjusted <- CX_report_file1
CX_report_file1_adjusted$readsM <- round(CX_report_file1$readsM - CX_report_file1$readsN * (1-conversion_file1))
CX_report_file1_adjusted$readsM[CX_report_file1_adjusted$readsM < 0 ] <- 0

CX_report_file2_adjusted <- CX_report_file2
CX_report_file2_adjusted$readsM <- round(CX_report_file2$readsM - CX_report_file2$readsN * (1-conversion_file2))
CX_report_file2_adjusted$readsM[CX_report_file2_adjusted$readsM < 0 ] <- 0


#Now, correcting total number of reads based on conversion rate
CX_report_file1_adjusted$readsN <- round(CX_report_file1$readsN * conversion_file1)
CX_report_file2_adjusted$readsN <- round(CX_report_file2$readsN * conversion_file2)

#Then, a new CX report can be generated using DMRcaller:
DMRcaller::saveBismark(CX_report_file1_adjusted,"CX_report_file1_adjusted.txt") 
DMRcaller::saveBismark(CX_report_file2_adjusted,"CX_report_file2_adjusted.txt")

#From now on we're going to work with the adjusted CX reports.

#------------------------------------------------------------------------------------------------------------------------
#CALLING DMRs IN THE 3 SEQUENCE CONTEXTS USING BIN METHOD

#Defining a smaller region of 30kb to test out this script:
#chr_local <- GRanges(seqnames = Rle("Chr1"), ranges = IRanges(1,30000))

#ben cutoffs = .35, .2, .15
ben <- list(.35, .2, .15)
#jacobsen cutoffs = .4, .2, .15
jacobsen <- list(.4, .2, .15)
cutoff <- ben


{DMRsBinsCG <- DMRcaller::computeDMRs(
	CX_report_file1_adjusted,
	CX_report_file2_adjusted,
context = "CG", 
method = "bins", 
binSize = 300,
pValueThreshold = 0.01, 
minCytosinesCount = 10,
minProportionDifference = cutoff[[1]], 
minGap = 500,
minReadsPerCytosine = 5, 
#regions = chr_local, #testing on small region
cores = 1)}
#They say that their adjusted P-value of 0.01 using the Benjamini and Hochberg's method for controlling false discovery rate
#Is minCytosineCount of 4 good, or choose another value? That's just the min. number of cytosines we need
#to call a region a DMR.

{DMRsBinsCHG <- DMRcaller::computeDMRs(
	CX_report_file1_adjusted, 
	CX_report_file2_adjusted,
context = "CHG", 
method = "bins", 
binSize = 300,
pValueThreshold = 0.01, 
minCytosinesCount = 5,
minProportionDifference = cutoff[[2]], 
minGap = 500,
minReadsPerCytosine = 5, 
#regions = chr_local, 
cores = 1)}

{DMRsBinsCHH <- DMRcaller::computeDMRs(
	CX_report_file1_adjusted, 
	CX_report_file2_adjusted,
context = "CHH", 
method = "bins", 
binSize = 300,
pValueThreshold = 0.01, 
minCytosinesCount = 5,
minProportionDifference = cutoff[[3]], 
minGap = 500,
minReadsPerCytosine = 5, 
#regions = chr_local, 
cores = 1)}


#---------------------------------------------------------------------------------------------------------

#EXPORTING THE RESULTING DMRs AS BED FILES

write.table(as.data.frame(DMRsBinsCG), file = paste0(outputDir, outputName, "_DMRsBinsCG.txt"), sep = "\t", quote = FALSE)
write.table(as.data.frame(DMRsBinsCG), file = paste0(outputDir, outputName, "_DMRsBinsCG.bed"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



write.table(as.data.frame(DMRsBinsCHG), file = paste0(outputDir, outputName, "_DMRsBinsCHG.txt"), sep = "\t", quote = FALSE)
write.table(as.data.frame(DMRsBinsCHG), file = paste0(outputDir, outputName, "_DMRsBinsCHG.bed"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


write.table(as.data.frame(DMRsBinsCHH), file = paste0(outputDir, outputName, "_DMRsBinsCHH.txt"), sep = "\t", quote = FALSE)
write.table(as.data.frame(DMRsBinsCHH), file = paste0(outputDir, outputName, "_DMRsBinsCHH.bed"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



