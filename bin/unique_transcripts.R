#!/usr/bin/env Rscript

### Author: Jun Hurtado

### This program selects one transcript per gene (after
### having selected one protein per transcript): the one
### that codes the longest protein.
###
### Starting from the proteins and cds sequences (fasta
### files) selected from Transdecoder output, a table of
### gene/transcript/proteins-ids correspondences is
### generated: one protein and transcript per parental
### gene. Protein length is added to this table.
###
### Run this script from Terminal with, e.g.:
###
###   $ Rscript --vanilla unique_transcripts.R \
###       -a filtered_transdecoder.pep.fasta \
###       -b filtered_transdecoder.cds.fasta \
###       -c target_transcripts.fasta \
###       -o selected
###


### LIBRARIES
library("optparse")
library("stringr")
library("seqinr")

### ARGUMENTS
option_list <- list(
  make_option(c("--input1", "-a"), default="filtered_transdecoder.pep.fasta",
              help = "name of the file with the selected translations from transdecoder"),
  make_option(c("--input2", "-b"), default="filtered_transdecoder.cds.fasta",
              help = "name of the file with the selected cds from transdecoder"),
  make_option(c("--input3", "-c"), default="target_transcripts.fasta",
              help = "transcripts fasta file name"),
  make_option(c("--outprefix", "-o"), default="NULL",
              help = "prefix for the output fasta files")
)
opt <- parse_args(OptionParser(option_list=option_list))


# Validating arguments
if (is.null(opt$input1) || !file.exists(opt$input1)) {
  stop(paste0("First input file not provided or does not exist. --input1 ", opt$input1))
}
if (is.null(opt$input2) || !file.exists(opt$input2)) {
  stop(paste0("Second input file not provided or does not exist. --input2 ", opt$input2))
}
if (is.null(opt$input3) || !file.exists(opt$input3)) {
  stop(paste0("Third input file not provided or does not exist. --input3 ", opt$input3))
}
if (is.null(opt$outprefix)) {
  stop(paste0("Prefix for the output fasta files must be provided. --outprefix ", opt$outprefix))
}


### CONSTANTS
# Read Files
input1_file <- paste0(opt$input1)
input2_file <- paste0(opt$input2)
input3_file <- paste0(opt$input3)
# Write Files
output1_file <- paste0(opt$outprefix, "transcripts.fasta")
output2_file <- paste0(opt$outprefix, "CDSs.fasta")
output3_file <- paste0(opt$outprefix, "translations.fasta")
output4_file <- "gn_tr_pp.txt"
output5_file <- paste0(opt$outprefix, "gn_tr_pp.txt")


### MAIN PROGRAM
Trs <- read.fasta(file=input3_file, as.string=TRUE)
Pps <- read.fasta(file=input1_file, as.string=TRUE)
N <- length(Pps)
Pps_Trs <- character(N)
for (i in 1:N) {
	Pps_Trs[i] <- sub("\\..*", "", names(Pps)[i])
	print(paste0("First loop progress ", round(i*100/N, 4), "%"))
}

N <- length(Trs)
Corr <- data.frame(gnID=character(N), trID=character(N), ppID=character(N))
levels(Corr$trID) <- names(Trs)
levels(Corr$ppID) <- c(names(Pps), "remove")
for (i in 1:N) {
	levels(Corr$gnID) <- c(levels(Corr$gnID), sub(".*geneID=", "", attr(Trs[[i]], "Annot")))
	Corr$gnID[i] <- sub(".*geneID=", "", attr(Trs[[i]], "Annot"))
	Corr$trID[i] <- names(Trs)[i]
	if(Corr$trID[i] %in% Pps_Trs) {
		Corr$ppID[i] <- names(Pps)[which(Pps_Trs == Corr$trID[i])]
	} else {
		Corr$ppID[i] <- "remove"
	}
	print(paste0("Second loop progress ", round(i*100/N, 4), "%"))
}
Corr <- Corr[Corr$ppID != "remove",]
Trs <- Trs[names(Trs) %in% Corr$trID]
names(Pps) <- Pps_Trs

Cds <- read.fasta(file=input2_file, as.string=TRUE)
Cds <- Cds[names(Cds) %in% Corr$ppID]
N <- length(Cds)
for (i in 1:N) {
	names(Cds)[i] <- as.character(Corr$trID[Corr$ppID == names(Cds)[i]])
	print(paste0("Third loop progress ", round(i*100/N, 4), "%"))
}

write.table(Corr, file=output4_file, row.names=FALSE, sep="\t", quote=FALSE)

N <- length(Corr[,1])
Corr$ppL <- integer(N)
for (i in 1:N) {
	Corr$ppL[i] <- nchar(Pps[which(names(Pps) == Corr$trID[i])])
	print(paste0("Fourth loop progress ", round(i*100/N, 4), "%"))
}

filtro <- logical(N)
niveles <- unique(Corr$gnID)
N <- length(niveles)
for (i in 1:N) {
	pos <- which(Corr$gnID == niveles[i])
	if (length(pos) == 1) {
		filtro[pos] <- T
	} else {
		maximo <- max(Corr$ppL[pos])
		pos2 <- which(Corr$ppL[pos] == maximo)[1]
		filtro[pos[pos2]] <- T
	}
	print(paste0("Fifth loop progress ", round(i*100/N, 4), "%"))
}
Corr2 <- Corr[filtro,]
write.table(Corr2, file=output5_file, row.names=FALSE, sep="\t", quote=FALSE)

Trs <- Trs[names(Trs) %in% Corr2$trID]
write.fasta(Trs, names=names(Trs), as.string=TRUE, file.out=output1_file)
Cds <- Cds[names(Cds) %in% Corr2$trID]
write.fasta(Cds, names=names(Cds), as.string=TRUE, file.out=output2_file)
Pps <- Pps[names(Pps) %in% Corr2$trID]
write.fasta(Pps, names=names(Pps), as.string=TRUE, file.out=output3_file)
