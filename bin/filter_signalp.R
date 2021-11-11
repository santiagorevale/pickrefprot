#!/usr/bin/env Rscript

### Author: Juan Hurtado

### This program filters protein sequences that may have
### signal peptide
###
### Only proteins with signal peptide according to signalP
### v5.0 or with incomplete N-terminus (5prime_partial) will
### pass the filter.
### Protein type (complete, 5prime_partial, 3prime_partial,
### or internal) is obtained from the table of
### gene/transcript/protein-ids correspodences.
###
### Run this script from Terminal with, e.g.:
###
###   $ Rscript --vanilla filter_signalp.R \
###       -f selected_translations.fasta \
###       -s summary.signalp5 \
###       -c gn_tr_pp.annotated.txt \
###       -o SigPeps.fasta
###


### LIBRARIES
library("optparse")
library("stringr")
library("seqinr")

### ARGUMENTS
option_list <- list(
  make_option(c("--translations", "-f"), default="selected_translations.fasta",
              help = "translations file name"),
  make_option(c("--signalp", "-s"), default="summary.signalp5",
              help = "signalP summary output file"),
  make_option(c("--correspondences", "-c"), default="gn_tr_pp.annotated.txt",
              help = "type-annotated proteins table"),
  make_option(c("--output", "-o"), default="SigPeps.fasta",
              help = "output file")
)
opt <- parse_args(OptionParser(option_list=option_list))



# Validating arguments
if (is.null(opt$translations) || !file.exists(opt$translations)) {
  stop(paste0("Translations file not provided or does not exist. --translations ", opt$translations))
}
if (is.null(opt$signalp) || !file.exists(opt$signalp)) {
  stop(paste0("SignalP Summary file not provided or does not exist. --signalp ", opt$signalp))
}
if (is.null(opt$correspondences) || !file.exists(opt$correspondences)) {
  stop(paste0("Correspondences Table file not provided or does not exist. --correspondences ", opt$correspondences))
}


### CONSTANTS
input1_file <- opt$translations
input2_file <- opt$signalp
input3_file <- opt$correspondences
output_file <- opt$output


### MAIN PROGRAM
translations <- read.fasta(file=input1_file, as.string=TRUE)

signalp <- read.table(input2_file, header=FALSE, fill=TRUE)
signalp <- signalp[signalp[,2] != "OTHER",]

correspondences <- read.table(input3_file, header=TRUE)
correspondences <- correspondences[correspondences$ppType == "5prime_partial",]

filter <- unique(c(as.character(signalp[,1]), as.character(correspondences$trID)))

translations <- translations[names(translations) %in% filter]
write.fasta(translations, names=names(translations), as.string=TRUE, file.out=output_file)
