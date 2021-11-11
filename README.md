# PickRepProt

PickRepProt is a pipeline to predict a representative protein sequence per gene from transcriptomes for Drosophila species.

## Dependencies

This pipeline makes use of some third-party software, packages, and databases. The versions provided below are the ones that were used to run the analysis for the publication.

### Software

- nextflow v21.04.3
- blast v2.9.0
- hmmer v3.3.0
- signalp v5.0b
- transdecoder v5.5.0
- r-base v3.6.1

### R packages

- optparse v1.6.4
- stringr v1.4.0
- seqinr v3.6-1


### Databases

- Pfam v33.1


## Prerequisites

- Because SignalP can not be freely distributed, get your own download of version 5.0b from [here](https://services.healthtech.dtu.dk/service.php?SignalP-5.0).

- If you DON'T have all the above mentioned software and packages available in your environment, you could run the pipeline using CONDA. A minimal installer can be downloaded from [here](https://docs.conda.io/en/latest/miniconda.html).


## Installation

1. Get [NextFlow](https://www.nextflow.io/) to run the pipeline.

```bash
curl -s https://get.nextflow.io | bash
```

2. Clone the repo:

```bash
git clone https://github.com/santiagorevale/pickrepprot
```

3. Download and unpack PFam database:

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
```

4. If you don't have `signalp` installed, then unpack it, place the folder on the `assets` directory and create a link to its executable file:

```bash
# you should have signalp-5.0b pre-downloaded in this location
tar xzf signalp-5.0b.Linux.tar.gz
mv signalp-5.0b pickrepprot/assets/
ln -s ../assets/signalp-5.0b/bin/signalp pickrepprot/bin/signalp
```


## Usage

Here's an example command on how to run the pipeline:

```bash
nextflow run \
    pickrepprot/main.nf \
    -profile conda \
    --ref_pfam /path/to/Pfam-A.hmm \
    --ref_translations /path/to/drosophila_translations.fasta \
    --transcripts /path/to/target_transcripts.fasta \
    --genome Dmoj
```


## Output

By default, pipeline results will be stored in the folder `results` in the same path where the pipeline was executed.

A temporary folder named `work` is created as part of the execution. If the pipeline finishes successfully, you could safely remove this folder.

The final output will be located in the `output_files` folder and it will be comprised by the following files:

```bash
# Gene/transcript/protein correspondence table
gn_tr_pp.txt

# Selected transcript, cds, and protein per gene
transcripts.fasta
cds.fasta
proteins.fasta

# Subgroup sequences that are ACP candidates
acps_transcripts.fasta
acps_cds.fasta
acps_proteins.fasta
```

## Troubleshooting

If the pipeline ends with an error, once fixed, you would be able to continue the pipeline from where it finished by adding the `-resume` parameter to the previous command and re-running it. Remember NOT to remove the `work` folder if you want to continue the pipeline from where it ended.

```bash
nextflow run \
    pickrepprot/main.nf \
    -profile conda \
    --ref_pfam /path/to/Pfam-A.hmm \
    --ref_translations /path/to/drosophila_translations.fasta \
    --transcripts /path/to/target_transcripts.fasta \
    --genome Dmoj
    -resume
```


## License

[MIT](LICENSE)


## Authors

- [@santiagorevale](https://www.github.com/santiagorevale)
- [@hurtadojuan](https://github.com/hurtadojuan)


## Citation

If you use PickRepProt for your analysis, you can cite the publication as follows:

> **Research gaps and new insights in the evolution of Drosophila seminal fluid proteins.**
>
> Hurtado, Juan; Almeida, Francisca; Belliard, Silvina; Revale, Santiago; Hasson, Esteban.
>
> _Insect Molecular Biology_ DATE_TBD. doi: [TBD]().


## References

### [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

### Pipeline tools

* [NCBI-Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    > Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. J. Mol. Biol. 1990. 215:403-410.

* [HMMER](http://hmmer.org/)
    > HMMER 3.3 (2019) Copyright (C) (2019) Howard Hughes Medical Institute. Freely distributed under the BSD open source license. http://hmmer.org/
    > Eddy SR. Accelerated profile HMM searches. PLoS Comp. Biol. 2011. 7:e1002195.

* [SignalP](https://services.healthtech.dtu.dk/service.php?SignalP-5.0)
    > Almagro Armenteros JJ, Tsirigos KD, Sønderby CK, Petersen TN, Winther O, Brunak S, et al. SignalP 5.0 improves signal peptide predictions using deep neural networks. Nat. Biotechnol. 2019. 37(4):420–423.

* [Transdecoder](http://transdecoder.github.io)
    > Brian H, Papanicolaou A. Transdecoder (Find Coding Regions Within Transcripts). (n.d.) http://transdecoder.github.io

* [R Project](https://www.r-project.org/)
    > R Core Team (2019) R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

### R packages

* [optparse](https://CRAN.R-project.org/package=optparse)
    > Davis TL. (2019). optparse: Command Line Option Parser. R package version 1.6.4.

* [stringr](https://CRAN.R-project.org/package=stringr)
    > Wickham H. (2019). stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.4.0.

* [seqinr](https://cran.r-project.org/package=seqinr)
    > Charif D, Lobry JR. (2007) SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis. In: Structural approaches to sequence evolution: Molecules, networks, populations (Eds. Bastolla U., Porto, M., Roman, H.E, and Vendruscolo, M.) pp 207-232. Springer Verlag, New York.

### Software packaging/containerisation tools

* [Anaconda](https://anaconda.com)
    > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

* [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)
    > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

### Databases

* [Pfam]
    > Bateman A, Birney E, Durbin R, Eddy SR, Finn RD, Sonnhammer ELL. Pfam 3.1: 1313 multiple alignments and proﬁleHMMs match the majority of proteins. Nucleic Acids Res. 1999. 27:260–262.
    > Mistry J, Chuguransky S, Williams L, Qureshi M, Salazar GA, Sonnhammer ELL, Tosatto SCE, Paladin L, Raj S, Richardson LJ, Finn RD, Bateman A. Pfam: The protein families database in 2021. Nucleic Acids Res. 2021. 49(D1):D412–D419. https://doi.org/10.1093/nar/gkaa913
