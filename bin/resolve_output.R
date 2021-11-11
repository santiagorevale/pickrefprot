#!/usr/bin/env Rscript

### Author: Juan Hurtado

### This program generates output files of PickRepProt
### pipeline.
###
### Proteins present as subjects in a blast output table are
### considered to be Acps (seminal proteins produced by the
### male accessory glands). Based on this, Acp condition is
### added to the table of gene/transcript/protein-ids
### correspondences. And fasta files of Acps transcripts,
### cds, and proteins are generated in the output folder. If
### no blast table is provided the Acp condition is not
### evaluated and no such fasta files are generated.
### The final table of gene/transcript/protein-ids
### correspondences and fasta files of the predicted,
### selected proteins are placed in an output folder. 
### A (species) prefix can be added to each sequence name.
### 
### Run this script from Terminal with, e.g.:
###
###   $ Rscript --vanilla resolve_output.R \
###       -a filtered_blastx.outftm6 \
###       -b selected_transcripts.fasta \
###       -c selected_CDSs.fasta \
###       -d selected_translations.fasta \
###       -e gn_tr_pp.annotated.txt 
###       -f Acp -p Dmoj
###


### LIBRARIES
library("optparse")
library("stringr")
library("seqinr")


### ARGUMENTS
option_list <- list(
  make_option(c("--blastfile", "-a"), default="NULL",
              help = "Blast output file name"),
  make_option(c("--transcripts", "-b"), default="NULL",
              help = "Transcripts fasta file name"),
  make_option(c("--cds", "-c"), default="NULL",
              help = "CDSs fasta file name"),
  make_option(c("--proteins", "-d"), default="NULL",
              help = "Translations fasta file name"),
  make_option(c("--correspondences", "-e"), default="NULL",
              help = "gene/transcript/translationIDs correspondence table file name"),
  make_option(c("--output", "-o"), default="output_files",
              help = "Output folder"),
  make_option(c("--filter", "-f"), default="blast_hit",
              help = "Match criteria: e.g. Acp, Transcription_factor, etc"),
  make_option(c("--prefix", "-p"), default="NULL",
              help = "Name of the target transcripts fasta file")
)
opt <- parse_args(OptionParser(option_list=option_list))



# Validating arguments
if (!is.null(opt$blastfile) && file.exists(opt$blastfile) && file.info(opt$blastfile)$size > 0) {
	warning("Blast results being used for filtering and annotating sequences.")
} else {
	warning("No blast results being used for filtering or annotating sequences, because no file or an empty were provided.")
	opt$blastfile <- NULL
}
if (is.null(opt$transcripts) || !file.exists(opt$transcripts)) {
  stop(paste0("Second input file not provided or does not exist. --transcripts ", opt$transcripts))
}
if (is.null(opt$cds) || !file.exists(opt$cds)) {
  stop(paste0("Third input file not provided or does not exist. --cds ", opt$cds))
}
if (is.null(opt$proteins) || !file.exists(opt$proteins)) {
  stop(paste0("Fourth input file not provided or does not exist. --proteins ", opt$proteins))
}
if (is.null(opt$correspondences) || !file.exists(opt$correspondences)) {
  stop(paste0("Fifth input file not provided or does not exist. --correspondences ", opt$correspondences))
}
if (is.null(opt$output) || dir.exists(opt$output)) {
  stop(paste0("Output folder not provided or does already exist. --output ", opt$output))
}
if (is.null(opt$prefix)) {
  warning("No prefix provided. Are you sure you don't want to use prefixes for sequences IDs/names")
}


### CONSTANTS
# Parameters
filter  <- opt$filter
prefix  <- opt$prefix
queries <- opt$queries
output_folder <- opt$output
# Read Files
input1_file <- opt$blastfile
input2_file <- opt$transcripts
input3_file <- opt$cds
input4_file <- opt$proteins
input5_file <- opt$correspondences
# Write Files
output1A_file <- paste0(output_folder, "/transcripts.fasta")
output1B_file <- paste0(output_folder, "/acps_transcripts.fasta")
output2A_file <- paste0(output_folder, "/cds.fasta")
output2B_file <- paste0(output_folder, "/acps_cds.fasta")
output3A_file <- paste0(output_folder, "/proteins.fasta")
output3B_file <- paste0(output_folder, "/acps_proteins.fasta")
output4_file  <- paste0(output_folder, "/gn_tr_pp.txt")


### MAIN PROGRAM
dir.create(output_folder)

Corr <- read.table(input5_file, header=TRUE)
Corr$gnID <- paste0(prefix, ":", Corr$gnID)
Corr$trID <- paste0(prefix, ":", Corr$trID)
Corr$ppID <- paste0(prefix, ":", Corr$ppID)

if (!is.null(input1_file)) {
	data <- read.table(input1_file, header=FALSE)
	names(data) <- c("qseqid", "qlen", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs") # Nombramos las columnas
	subjects <- paste0(prefix, ":", unique(data$sseqid))
}


N <- length(Corr[,1])
Corr$filtro <- logical(N)
if (!is.null(input1_file)) {
	for(i in 1:N) {
		if(Corr$trID[i] %in% subjects) {
			Corr$filtro[i] <- T
		}
		print(paste0("Loop progress ", round(i*100/N, 4), "%"))
	}
} else {
	Corr$filtro <- rep(NA, N)
}
N <- length(names(Corr))
names(Corr)[N] <- filter
write.table(Corr, file=output4_file, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

Trs <- read.fasta(file=input2_file, as.string=TRUE)
names(Trs) <- paste0(prefix, ":", names(Trs))
write.fasta (Trs, names=names(Trs), as.string=TRUE, file.out=output1A_file)
if (!is.null(input1_file)) {
	Trs <- Trs[names(Trs) %in% subjects]
	write.fasta (Trs, names=names(Trs), as.string=TRUE, file.out=output1B_file)
}
Cds <- read.fasta(file=input3_file, as.string=TRUE)
names(Cds) <- paste0(prefix, ":", names(Cds))
write.fasta (Cds, names=names(Cds), as.string=TRUE, file.out=output2A_file)
if (!is.null(input1_file)) {
	Cds <- Cds[names(Cds) %in% subjects]
	write.fasta (Cds, names=names(Cds), as.string=TRUE, file.out=output2B_file)
}
Pps <- read.fasta(file=input4_file, as.string=TRUE)
names(Pps) <- paste0(prefix, ":", names(Pps))
write.fasta (Pps, names=names(Pps), as.string=TRUE, file.out=output3A_file)
if (!is.null(input1_file)) {
	Pps <- Pps[names(Pps) %in% subjects]
	write.fasta (Pps, names=names(Pps), as.string=TRUE, file.out=output3B_file)
}
