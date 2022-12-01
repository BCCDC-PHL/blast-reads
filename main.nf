#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp }               from './modules/blast_reads.nf'
include { seqtk_seq }           from './modules/blast_reads.nf'
include { blastn as blastn_r1 } from './modules/blast_reads.nf'
include { blastn as blastn_r2 } from './modules/blast_reads.nf'
include { choose_best_hsp_per_query as choose_best_hsp_per_query_r1 } from './modules/blast_reads.nf'
include { choose_best_hsp_per_query as choose_best_hsp_per_query_r2 } from './modules/blast_reads.nf'
include { csvtk_freq as csvtk_freq_r1 } from './modules/blast_reads.nf'
include { csvtk_freq as csvtk_freq_r2 } from './modules/blast_reads.nf'
include { combine_counts }              from './modules/blast_reads.nf'


workflow {

  if (params.samplesheet_input != 'NO_FILE') {
    ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
  } else {
    ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
  }

  ch_blast_db_dir = Channel.fromPath("${params.blast_db_dir}", type: 'dir')
  ch_blast_db_name = Channel.of("${params.blast_db_name}")

  main:

    fastp(ch_fastq)

    seqtk_seq(fastp.out.reads)

    ch_blastn_r1 = blastn_r1(seqtk_seq.out.map{ it -> [it[0], it[1]] }.combine(ch_blast_db_dir).combine(ch_blast_db_name).combine(Channel.of("R1")))

    ch_blastn_r2 = blastn_r2(seqtk_seq.out.map{ it -> [it[0], it[2]] }.combine(ch_blast_db_dir).combine(ch_blast_db_name).combine(Channel.of("R2")))

    ch_blastn_best_hsp_r1 = choose_best_hsp_per_query_r1(ch_blastn_r1)

    ch_blastn_best_hsp_r2 = choose_best_hsp_per_query_r2(ch_blastn_r2)

    ch_blast_counts_r1 = csvtk_freq_r1(ch_blastn_best_hsp_r1)

    ch_blast_counts_r2 = csvtk_freq_r2(ch_blastn_best_hsp_r2)

    combine_counts(ch_blast_counts_r1.map{ it -> [it[0], it[1]] }.join(ch_blast_counts_r2.map{ it -> [it[0], it[1]] }))
    
}
