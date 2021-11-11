#!/usr/bin/env Rscript

### Author: Juan Hurtado

### This program filters Transdecoder predictions
###
### Portions of CDS (if provided) and proteins with unknown nucleotides or aa 
### (Ns o Xx) will be eliminated, only the longest remaining fragment will be 
### kept, and the protein type (complete, 5prime_partial, 3prime_partial, or
### internal) will be reannotated. 
### Spurious prediction (i.e., incomplete proteins of less than 'l' aa) will be 
### eliminated. 
### Up to two proteins per transcript wil be preselected as potentially correct:
### the longest without hits against protein database and the longest with hits.
### If the longest of these preselected proteins is the one with hits, the other
### The shorter one will be eliminated unless 1) it is the one with hits and 2)
### it is less than 'r' times shorter than the longer one. In such a case, the
### longer one will be eliminated. 
### Fasta files with the filtered sequences are generated.
###
### Run this script from Terminal with, e.g.:
###
###   $ Rscript --vanilla filter_transdecoder.R \
###       -l 100 -r 3 \
###       -p target_transcripts.fasta.transdecoder.pep \
###       -o filtered_transdecoder.pep.fasta
###


### LIBRARIES
library("optparse")
library("stringr")
library("seqinr")

### ARGUMENTS
option_list <- list(
  make_option(c("--length", "-l"), default=100,
              help = "minimum translation length of incomplete proteins for retaining"),
  make_option(c("--ratio", "-r"), default=3,
              help = "minimum ratio of protein length between the longer protein and the one with hits for retaining the one without hits"),
  make_option(c("--pep", "-p"), default="target_transcripts.fasta.transdecoder.pep",
              help = "name of the input proteins fasta file"),
  make_option(c("--cds", "-c"), default="target_transcripts.fasta.transdecoder.cds",
              help = "name of the input cds fasta file"),
  make_option(c("--output", "-o"), default="filtered_transdecoder.pep.fasta",
              help = "output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))


# Validating arguments 
if (is.null(opt$pep) || !file.exists(opt$pep)) {
  stop(paste0("pep input file not provided or does not exist. --pep ", opt$pep))
}
if (is.null(opt$cds) || !file.exists(opt$cds)) {
  stop(paste0("cds input file not provided or does not exist. --cds ", opt$cds))
}
if (is.null(opt$output)) {
  stop("Output prefix not provided. --output")
}


### CONSTANTS
# Parameters
length <- as.integer(opt$length)
ratio <- as.numeric(opt$ratio)
# Read Files
input_pep <- paste0(opt$pep)
input_cds <- paste0(opt$cds)
# Write Files
output_pep <- paste0(opt$output, input_pep)
output_cds <- paste0(opt$output, input_cds)



### MAIN PROGRAM
if (!is.null(opt$cds)) {
	CDSs <- read.fasta(file=input_cds, as.string=TRUE)
}
Peps <- read.fasta(file=input_pep, as.string=TRUE)
N <- length(Peps)
filtro <- !logical(N)
for (i in 1:N) {
	x <- attr(Peps[[i]], "Annot")
	if (grepl("X", Peps[[i]], fixed=T) || grepl("x", Peps[[i]], fixed=T)) {
		Peps[[i]] <- gsub("x", "X", Peps[[i]])
		while (grepl("XX", Peps[[i]], fixed=T)) {
			Peps[[i]][1] <- gsub("XX", "X", Peps[[i]][1])
		}
		y <- gregexpr("X", Peps[[i]])
		y <- c(0, y[[1]][seq(1:length(y[[1]]))], nchar(Peps[[i]]) + 1)
		d <- diff(y)
		pos <- which(d == max(d))[1]
		pos1 <- y[pos] + 1
		pos2 <- y[pos + 1] - 1
		Peps[[i]] <- substr(Peps[[i]], pos1, pos2)
		j <- which(names(CDSs) == names(Peps)[i])
		CDSs[[j]] <- substr(CDSs[[j]], 3*pos1 - 2, 3*pos2)
		if (sub(".*?type:(.*?) len.*", "\\1", x) == "complete") {
			if (pos1 == 1) {
				x <- gsub("complete", "3prime_partial", x)	
			} else if (grepl("*", Peps[[i]], fixed=T)) {
				x <- gsub("complete", "5prime_partial", x)
			} else {
				x <- gsub("complete", "internal", x)
			}
			attr(Peps[[i]], "Annot") <- paste0(x, " 'note: given part of the sequence was unknown only the longest known fragment was kept'")
			attr(CDSs[[j]], "Annot") <- paste0(x, " 'note: given part of the sequence was unknown only the longest known fragment was kept'")
		} else if (sub(".*?type:(.*?) len.*", "\\1", x) == "3prime_partial") {
			if (pos1 != 1) {
				x <- gsub("3prime_partial", "internal", x)
				attr(Peps[[i]], "Annot") <- paste0(x, " 'note: given part of the sequence was unknown only the longest known fragment was kept'")
				attr(CDSs[[j]], "Annot") <- paste0(x, " 'note: given part of the sequence was unknown only the longest known fragment was kept'")
			}
		} else if (sub(".*?type:(.*?) len.*", "\\1", x) == "5prime_partial") {
			if (!grepl("*", Peps[[i]], fixed=T)) {
				x <- gsub("5prime_partial", "internal", x)
				attr(Peps[[i]], "Annot") <- paste0(x, " 'note: given part of the sequence was unknown only the longest known fragment was kept'")
				attr(CDSs[[j]], "Annot") <- paste0(x, " 'note: given part of the sequence was unknown only the longest known fragment was kept'")
			}
		}
	}
	if (sub(".*?type:(.*?) len.*", "\\1", x) != "complete" && nchar(Peps[i]) < length) {
		filtro[i] <- F
	}
	print(paste0("First loop progress ", round(i*100/N, 4), "%"))
}
Peps <- Peps[filtro]

nombres <- names(Peps)
N <- length(Peps)
for (i in 1:N) {
	nombres[i] <- sub("\\..*", "", names(Peps)[i])
	print(paste0("Second loop progress ", round(i*100/N, 4), "%"))
}
Peps2 <- Peps
names(Peps2) <- nombres
filtro <- logical(N)

transcripts <- unique(nombres)
N <- length(transcripts)
for (i in 1:N) {
	pos <- which(names(Peps2) == transcripts[i])
	n <-length(pos)
	if (n == 1) {
		filtro[pos] <- T
	} else {
		matches <- logical(n)
		for (j in 1:n) {
			matches[j] <- grepl("|", attr(Peps[[pos[j]]], "Annot"), fixed=T)
		}
		if (sum(matches) > 0) {
			maximo <- max(nchar(Peps[pos[matches]]))
			pos2 <- which(nchar(Peps[pos[matches]]) == maximo)[1]
			filtro[pos[matches]][pos2] <- T
			if (sum(!matches) > 0) {
				maximo <- max(nchar(Peps[pos[!matches]]))
				pos2 <- which(nchar(Peps[pos[!matches]]) == maximo)[1]
				filtro[pos[!matches]][pos2] <- T
			}
		} else {
			maximo <- max(nchar(Peps[pos]))
			pos2 <- which(nchar(Peps[pos]) == maximo)[1]
			filtro[pos][pos2] <- T
		}
	}
	print(paste0("Third loop progress ", round(i*100/N, 4), "%"))
}

Peps <- Peps[filtro]
Peps2 <- Peps2[filtro]
transcripts <- unique(names(Peps2))
N <- length(Peps2)
filtro <- logical(N)

N <- length(transcripts)
for (i in 1:N) {
	pos <- which(names(Peps2) == transcripts[i])
	n <-length(pos)
	if (n == 1) {
		filtro[pos] <- T
	} else {
		maximo <- max(nchar(Peps[pos]))
		pos_max <- which(nchar(Peps[pos]) == maximo)[1]
		if (pos_max == 1) {
			pos_min <- 2
		} else {
			pos_min <- 1
		}
		if (nchar(Peps[pos][pos_max]) >= ratio*nchar(Peps[pos][pos_min])) {
			filtro[pos][pos_max] <- T
		} else if(grepl("|", attr(Peps[[pos[1]]], "Annot"), fixed=T)){
			filtro[pos][1] <- T
		} else {
			filtro[pos][2] <- T
		}
	}
	print(paste0("Fourth loop progress ", round(i*100/N, 4), "%"))
}
Peps <- Peps[filtro]

N <- length(Peps)
for (i in 1:N) {
	names(Peps)[i] <- substr(attr(Peps[[i]], "Annot"), 2, nchar(attr(Peps[[i]], "Annot")))
	print(paste0("Fifth loop progress ", round(i*100/N, 4), "%"))
}
write.fasta(Peps, names=names(Peps), as.string=TRUE, file.out=output_pep)

N <- length(CDSs)
for (i in 1:N) {
	names(CDSs)[i] <- substr(attr(CDSs[[i]], "Annot"), 2, nchar(attr(CDSs[[i]], "Annot")))
	print(paste0("Sixth loop progress ", round(i*100/N, 4), "%"))
}

if (!is.null(opt$cds)) {
	pos <- which(names(CDSs) %in% names(Peps))
	CDSs <- CDSs[pos]
	write.fasta(CDSs, names=names(CDSs), as.string=TRUE, file.out=output_cds)
} else {
	warning("CDS file was not provided: CDSs could not be filtered")
}
