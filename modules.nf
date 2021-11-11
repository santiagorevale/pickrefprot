params.options = [:]
options        = params.options ?: ""

process HMMER_HMMPRESS {
    tag "${hmm_file}"
    publishDir path: { params.save_reference ? "${params.outdir}/references/hmmer" : params.outdir },
               saveAs: { params.save_reference ? it : null }, mode: "copy"

    input:
    path hmm_file

    output:
    path "${hmm_file}.*", emit: db

    script:
    """
    hmmpress ${options.args} ${hmm_file}
    """
}


process BLAST_MAKEBLASTDB {
    tag "${fasta_file}"
    publishDir path: { params.save_reference ? "${params.outdir}/references/blast" : params.outdir },
               saveAs: { params.save_reference ? it : null }, mode: "copy"

    input:
    path fasta_file

    output:
    path "${fasta_file}.*", emit: db

    script:
    """
    makeblastdb ${options.args} -in ${fasta_file} -parse_seqids -out ${fasta_file}
    """
}


process TRANSDECODER_LONGORFS {
    tag "${transcripts}"

    input:
    path transcripts

    output:
    path "transdecoder.longorfs/longest_orfs.pep", emit: longest_orfs
    path "transdecoder.longorfs"                 , emit: folder

    script:
    """
    TransDecoder.LongOrfs ${options.args} -t ${transcripts} --output_dir transdecoder.longorfs
    """
}


process BLAST_BLASTP {
    tag "${index_base}, ${fasta_file}"

    input:
    path fasta_file
    path blast_db

    output:
    path "*.blastp.hits", emit: hits

    script:
    index_base  = blast_db[0].toString() - ~/\.\w+$/
    output_file = fasta_file.name - ~/\.\d+\.\w+$/ + ".blastp.hits"
    """
    blastp \
        ${options.args} \
        -query ${fasta_file} -db ${index_base} -out ${output_file} \
        -outfmt "6 qseqid qlen sseqid slen pident length mismatch \
        gapopen qstart qend sstart send evalue bitscore qcovs"
    """
}


process BLAST_BLASTX {
    tag "${index_base}, ${fasta_file}"

    input:
    path fasta_file
    path blast_db

    output:
    path "*.blastx.hits", emit: hits

    script:
    index_base  = blast_db[0].toString() - ~/\.\w+$/
    output_file = fasta_file.name - ~/\.\d+\.\w+$/ + ".blastp.hits"
    """
    blastx \
        ${options.args} \
        -query ${fasta_file} -db ${index_base} -out ${output_file} \
        -outfmt "6 qseqid qlen sseqid slen pident length mismatch \
        gapopen qstart qend sstart send evalue bitscore qcovs"
    """
}


process SCRIPTS_FILTERBLAST {
    tag "${blast_hits}"
    publishDir "${params.outdir}/blast", mode: "copy"

    input:
    path blast_hits

    output:
    path "${blast_hits}.filtered", emit: hits_filtered

    script:
    """
    filter_blast.R ${options.args} -i ${blast_hits} -o ${blast_hits}.filtered
    """
}


process HMMER_HMMSCAN {
    tag "${index_base}, ${fasta_file}"

    input:
    path fasta_file
    path hmm_db

    output:
    path "*.pfam.hits", emit: hits

    script:
    index_base  = hmm_db[0].toString() - ~/\.\w+$/
    output_file = fasta_file.name - ~/\.\d+\.\w+$/ + ".pfam.hits"
    """
    hmmscan \
        ${options.args} \
        --cpu ${task.cpus} \
        -o /dev/null --domtblout ${output_file} ${index_base} ${fasta_file}
    sed -i '/^#/d' ${output_file}
    """
}


process TRANSDECODER_PREDICT {
    tag "${transcripts}"
    publishDir "${params.outdir}/transdecoder", mode: "copy"

    input:
    path  transcripts
    path  blast_hits
    path  pfam_hits
    path  transdecoder_longorfs

    output:
    path  "${transcripts}.transdecoder.{bed,gff3}"
    tuple path("${transcripts}.transdecoder.cds"), path("${transcripts}.transdecoder.pep"), emit: file_pair

    script:
    """
    TransDecoder.Predict \
        ${options.args} \
        -t ${transcripts} \
        --retain_pfam_hits ${pfam_hits} \
        --retain_blastp_hits ${blast_hits} \
        --output_dir ${transdecoder_longorfs}
    """
}


process SCRIPTS_FILTERTRANSDECODER {
    tag "transcripts.fasta, cds.fasta, translations.fasta"
    publishDir "${params.outdir}/transdecoder.filtered", mode: "copy"

    input:
    path  "transcripts.fasta"
    tuple path("cds.fasta"), path("translations.fasta")

    output:
    path "selected.translations.fasta", emit: translations
    path "selected.CDSs.fasta"        , emit: cds
    path "selected.transcripts.fasta" , emit: transcripts
    path "gn_tr_pp.annotated.txt"     , emit: table

    script:
    """
    filter_transdecoder.R \
        ${options.args} \
        -o filtered. \
        -c cds.fasta \
        -p translations.fasta
    
    unique_transcripts.R \
        -o selected. \
        -c transcripts.fasta \
        -b filtered.cds.fasta \
        -a filtered.translations.fasta

    annotate_correspondences.R \
        -a filtered.translations.fasta \
        -b selected.gn_tr_pp.txt \
        -o gn_tr_pp.annotated.txt
    """
}

process SIGNALP_SIGNALP {
    tag "$translations"
    publishDir "${params.outdir}/signalp", mode: "copy"

    input:
    path translations

    output:
    path "summary.signalp5", emit: summary

    script:
    """
    signalp ${options.args} -fasta ${translations} -stdout > summary.signalp5
    """
}


process SCRIPTS_FILTERSIGNALP {
    tag "$translations"
    publishDir "${params.outdir}/signalp", mode: "copy"

    input:
    path translations
    path gntrpp
    path summary

    output:
    path "peptides.signalp.fasta", emit: fasta

    script:
    """
    filter_signalp.R -f ${translations} -s ${summary} -c ${gntrpp} -o peptides.signalp.fasta
    """
}


process RESOLVE_OUTPUT {
    publishDir "${params.outdir}", mode: "copy"

    input:
    val  genome_prefix
    path transcripts
    path cds
    path translations
    path gntrpp
    path acps_blast_hits

    output:
    path 'output_files'

    script:
    def option_acps_blast = acps_blast_hits ? "-a ${acps_blast_hits}" : ""
    """
    resolve_output.R \
        ${options.args} \
        ${option_acps_blast} \
        -p ${genome_prefix} \
        -b ${transcripts} \
        -c ${cds} \
        -d ${translations} \
        -e ${gntrpp}
    """
}
