#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp }       from './modules/blast_reads.nf'
include { seqtk_seq }   from './modules/blast_reads.nf'
include { blastn }      from './modules/blast_reads.nf'
include { csvtk_freq }  from './modules/blast_reads.nf'


workflow {

  if (params.samplesheet_input != 'NO_FILE') {
    ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
  } else {
    ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
  }

  ch_blast_db_dir = Channel.fromPath("${params.blast_db_dir}", type: 'dir')
  ch_blast_db_name = Channel.of("${params.blast_db_name}")

  main:

    // fastp(ch_fastq)

    seqtk_seq(ch_fastq)

    blastn(seqtk_seq.out.map{ it -> [it[0], it[1]] }.combine(ch_blast_db_dir).combine(ch_blast_db_name))

    csvtk_freq(blastn.out)
}
