#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ========================================================================================
//     PRINT HELP
// ========================================================================================
if (params.help) {
    log.info"""
        Usage:

        The typical command for running the pipeline is as follows:

        nextflow run /path/to/pickrepprot \
            --ref_pfam <path_to_file> \
            --ref_translations <path_to_file> \
            --transcripts <path_to_file> \
            --genome <prefix>

        
        Mandatory arguments:
        --ref_pfam                    Path to Pfam-A.hmm file.
        --ref_translations            Path to FastA file with the reference translations.
        --transcripts                 Path to FastA file with the transcripts to analyse.
        --genome                      Genome prefix.

        Optional arguments:
        --ref_acps                    Path to FastA file with the reference ACPs.
        --fasta_chunks                Number of sequences to split a FastA file into to split the processing.
        --outdir                      Output directory for the pipeline (default: ./results)
    """.stripIndent()
    exit 0
}

// ========================================================================================
//     VALIDATE INPUTS
// ========================================================================================
if (!params.genome) {
    exit 1, "Genome prefix not provided: ${params.genome}."
}

if (!params.fasta_chunks.toString().isInteger()) {
    exit 1, "The value provided to --fasta_chunks is not a number: ${params.fasta_chunks}."
}

Channel
    .fromPath(params.transcripts, checkIfExists: true)
    .ifEmpty { exit 1, "Transcripts file not found: ${params.transcripts}" }
    .set { ch_transcripts }

Channel
    .fromPath(params.ref_pfam, checkIfExists: true)
    .ifEmpty { exit 1, "Reference PFAM file not found: ${params.ref_pfam}" }
    .set { ch_ref_pfam }

Channel
    .fromPath(params.ref_translations, checkIfExists: true)
    .ifEmpty { exit 1, "Reference translations file not found: ${params.ref_translations}" }
    .set { ch_ref_translations }

ch_acps = Channel.empty()
if (params.ref_acps) {
    Channel
        .fromPath(params.ref_acps, checkIfExists: true)
        .ifEmpty { exit 1, "Reference ACPs file not found: ${params.ref_acps}" }
        .splitFasta( by: params.fasta_chunks, file: true )
        .set { ch_acps }
}


// ========================================================================================
//     IMPORT LOCAL MODULES
// ========================================================================================

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_SIGNALP  } from './modules' addParams(options: modules['blast_makeblastdb_signalp'])
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_TRANS    } from './modules' addParams(options: modules['blast_makeblastdb_trans'])
include { BLAST_BLASTP as BLAST_BLASTP_ORFS               } from './modules' addParams(options: modules['blast_blastp_orfs'])
include { BLAST_BLASTX as BLAST_BLASTX_ACPS               } from './modules' addParams(options: modules['blast_blastx_acps'])
include { HMMER_HMMPRESS                                  } from './modules' addParams(options: modules['hmmer_hmmpress'])
include { HMMER_HMMSCAN                                   } from './modules' addParams(options: modules['hmmer_hmmscan'])
include { RESOLVE_OUTPUT                                  } from './modules' addParams(options: modules['resolve_output'])
include { SCRIPTS_FILTERBLAST as SCRIPTS_FILTERBLAST_ACPS } from './modules' addParams(options: modules['scripts_filterblast_acps'])
include { SCRIPTS_FILTERBLAST as SCRIPTS_FILTERBLAST_ORFS } from './modules' addParams(options: modules['scripts_filterblast_orfs'])
include { SCRIPTS_FILTERSIGNALP                           } from './modules' addParams(options: modules['scripts_filtersignalp'])
include { SCRIPTS_FILTERTRANSDECODER                      } from './modules' addParams(options: modules['scripts_filtertransdecoder'])
include { SIGNALP_SIGNALP                                 } from './modules' addParams(options: modules['signalp_signalp'])
include { TRANSDECODER_LONGORFS                           } from './modules' addParams(options: modules['transdecoder_longorfs'])
include { TRANSDECODER_PREDICT                            } from './modules' addParams(options: modules['transdecoder_predict'])
//include { GET_SOFTWARE_VERSIONS                      } from './modules'


// ========================================================================================
//     RUN MAIN WORKFLOW
// ========================================================================================
workflow {

    // Prepare Pfam DB
    HMMER_HMMPRESS(
        ch_ref_pfam
    )
    ch_pfam_db = HMMER_HMMPRESS.out.db.collect()

    // Prepare Translations Blast DB
    BLAST_MAKEBLASTDB_TRANS(
        ch_ref_translations
    )
    ch_blast_db_trans = BLAST_MAKEBLASTDB_TRANS.out.db.collect()

    // Get the longest ORFs per transcript
    TRANSDECODER_LONGORFS(
        ch_transcripts
    )
    ch_transdecoder_longorfs = TRANSDECODER_LONGORFS.out.folder
    ch_orfs                  = TRANSDECODER_LONGORFS.out.longest_orfs

    // Split the FastA file with the predicted longest ORFs
    ch_orfs
        .splitFasta( by: params.fasta_chunks, file: true )
        .set { ch_orfs_chunks }

    // Scan the predicted ORFs against Pfam
    HMMER_HMMSCAN(
        ch_orfs_chunks, 
        ch_pfam_db
    )
    ch_orfs_hmmer_hits = HMMER_HMMSCAN.out.hits.collectFile(storeDir: "${params.outdir}/hmmer")

    // Blast the predicted ORFs against the reference translated proteins
    BLAST_BLASTP_ORFS(
        ch_orfs_chunks, 
        ch_blast_db_trans
    )
    ch_orfs_blast_hits = BLAST_BLASTP_ORFS.out.hits.collectFile(storeDir: "${params.outdir}/blast")

    // Filter the Blast hits on the ORFs
    SCRIPTS_FILTERBLAST_ORFS(
        ch_orfs_blast_hits
    )
    ch_orfs_blast_hits_filtered = SCRIPTS_FILTERBLAST_ORFS.out.hits_filtered

    // Run the prediction with transdecoder plus the new evidences
    TRANSDECODER_PREDICT(
        ch_transcripts, 
        ch_orfs_blast_hits_filtered, 
        ch_orfs_hmmer_hits, 
        ch_transdecoder_longorfs
    )
    ch_transdecoder_predict = TRANSDECODER_PREDICT.out.file_pair


    SCRIPTS_FILTERTRANSDECODER(
        ch_transcripts,
        ch_transdecoder_predict
    )
    ch_filtertransdecoder_transcripts  = SCRIPTS_FILTERTRANSDECODER.out.transcripts
    ch_filtertransdecoder_cds          = SCRIPTS_FILTERTRANSDECODER.out.cds
    ch_filtertransdecoder_translations = SCRIPTS_FILTERTRANSDECODER.out.translations
    ch_filtertransdecoder_table        = SCRIPTS_FILTERTRANSDECODER.out.table

    // Only if a list of ACP proteins was provided
    ch_acps_blast_hits_filtered = Channel.empty()
    if (params.ref_acps) {

        SIGNALP_SIGNALP(
            ch_filtertransdecoder_translations
        )
        ch_signalp = SIGNALP_SIGNALP.out.summary

        SCRIPTS_FILTERSIGNALP(
            ch_filtertransdecoder_translations,
            ch_filtertransdecoder_table,
            ch_signalp
        )
        ch_signalp_filtered = SCRIPTS_FILTERSIGNALP.out.fasta

        // Prepare Blast DB for SingalP results
        BLAST_MAKEBLASTDB_SIGNALP(
            ch_signalp_filtered
        )
        ch_blast_db_signalp = BLAST_MAKEBLASTDB_SIGNALP.out.db.collect()

        // Blast the list of provided ACPs against the proteins with a hit against SignalP
        BLAST_BLASTX_ACPS(
            ch_acps, 
            ch_blast_db_signalp
        )
        ch_acps_blast_hits = BLAST_BLASTX_ACPS.out.hits.collectFile(storeDir: "${params.outdir}/blast")

        // Filter the Blast hits on the ORFs
        SCRIPTS_FILTERBLAST_ACPS(
            ch_acps_blast_hits
        )
        ch_acps_blast_hits_filtered = SCRIPTS_FILTERBLAST_ACPS.out.hits_filtered

    }

    // Prepare and save final files
    RESOLVE_OUTPUT(
        params.genome,
        ch_filtertransdecoder_transcripts,
        ch_filtertransdecoder_cds,
        ch_filtertransdecoder_translations,
        ch_filtertransdecoder_table,
        ch_acps_blast_hits_filtered.ifEmpty([])
    )

}





/*
 * Completion e-mail notification
 */
workflow.onComplete {
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    if (workflow.success) {
        log.info "${c_purple}[pickrepprot]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        log.info "${c_purple}[pickrepprot]${c_red} Pipeline completed with errors${c_reset}"
    }
}
