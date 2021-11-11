#!/usr/bin/env Rscript

### Author: Juan Hurtado

### This program filters ftm6 BLAST tables
###
### Run this script from Terminal with, e.g.:
### 
###   $ Rscript --vanilla filter_blast.R \
###       -e 1e-5 -b 80 -n yes \
###       -i blast.outfmt6 \
###       -o filtered_blast.outftm6
###


### LIBRARIES
library("optparse")
library("stringr")

### ARGUMENTS
option_list <- list(
  make_option(c("--input", "-i"), default=NULL,
              help = "blast table file name"),
  make_option(c("--qlen", "-q"), default=NULL,
              help = "minimum query length"),
  make_option(c("--slen", "-s"), default=NULL,
              help = "minimum subject length"),
  make_option(c("--pident", "-p"), default=NULL,
              help = "minimum pident: percentage of identical matches"),
  make_option(c("--length", "-l"), default=NULL,
              help = "minimum alignment length"),
  make_option(c("--evalue", "-e"), default=NULL,
              help = "maximum e-value"),
  make_option(c("--bitscore", "-b"), default=40,
              help = "minimum bitscore"),
  make_option(c("--qcovs", "-Q"), default=NULL,
              help = "minimum query coverage per subject"),
  make_option(c("--names", "-n"), default="no",
              help = "Want headers in output table? 'yes' or 'no'"),
  make_option(c("--output", "-o"), default="filtered_blast.outftm6",
              help = "output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))



# Validating arguments
if (is.null(opt$qlen) && is.null(opt$slen) && is.null(opt$pident) && is.null(opt$length) && is.null(opt$evalue) && is.null(opt$bitscore) && is.null(opt$qcovs)) {
  stop("No filter criteria provided")
}
if (is.null(opt$input) || !file.exists(opt$input)) {
  stop(paste0("Input file not provided or does not exist. --input ", opt$input))
}
if (opt$names != "yes" && opt$names != "no") {
	stop(paste0("Must choose 'yes' or 'no' for --names ", opt$names))
}


### CONSTANTS
# Parameters
qlen <- as.integer(opt$qlen)
slen <- as.integer(opt$slen)
pident <- as.numeric(opt$pident)
length <- as.integer(opt$length)
evalue <- as.numeric(opt$evalue)
bitscore <- as.numeric(opt$bitscore)
qcovs <- as.numeric(opt$qcovs)
# Read Files
input_file <- opt$input
# Write Files
output_file <- opt$output



### MAIN PROGRAM
data <- read.table(input_file, header=FALSE)
names(data) <- c("qseqid", "qlen", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs")

if (!is.null(opt$qlen)) {
	data <- data[data$qlen >= qlen,]
}
if (!is.null(opt$slen)) {
	data <- data[data$slen >= slen,]
}
if (!is.null(opt$pident)) {
	data <- data[data$pident >= pident,]
}
if (!is.null(opt$length)) {
	data <- data[data$length >= length,]
}
if (!is.null(opt$evalue)) {
	data <- data[data$evalue <= evalue,]
}
if (!is.null(opt$bitscore)) {
	data <- data[data$bitscore >= bitscore,]
}
if (!is.null(opt$qcovs)) {
	data <- data[data$qcovs >= qcovs,]
}

if (opt$names == "yes") {
	write.table(data, file=output_file, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
} else {
	write.table(data, file=output_file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
}

