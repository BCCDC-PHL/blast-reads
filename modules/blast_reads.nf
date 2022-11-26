process fastp {

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_fastp.{json,csv}"

  input:
  tuple val(sample_id), path(reads_1), path(reads_2)

  output:
  tuple val(sample_id), path("${sample_id}_fastp.{json,csv}"), emit: fastp_reports
  tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads

  script:
  """
  fastp \
    --cut_tail \
    -i ${reads_1} \
    -I ${reads_2} \
    -o ${sample_id}_trimmed_R1.fastq.gz \
    -O ${sample_id}_trimmed_R2.fastq.gz

  mv fastp.json ${sample_id}_fastp.json

  fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}_fastp.csv
  """
}


process seqtk_seq {

  tag { sample_id }

  input:
  tuple val(sample_id), path(reads_1), path(reads_2)

  output:
  tuple val(sample_id), path("${sample_id}_R1.fasta"), path("${sample_id}_R2.fasta")

  script:
  """
  seqtk seq -a ${reads_1} > ${sample_id}_R1.fasta
  seqtk seq -a ${reads_2} > ${sample_id}_R2.fasta
  """
}


process blastn {

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_blast.csv"

  input:
  tuple val(sample_id), path(reads), path(blast_db_dir), val(blast_db_name)

  output:
  tuple val(sample_id), path("${sample_id}_blast.csv")

  script:
  """
  export BLASTDB="${blast_db_dir}"
  blastn \
    -num_threads ${task.cpus} \
    -query ${reads} \
    -db ${blast_db_name} \
    -outfmt "10 qseqid sseqid sacc length pident gapopen gaps evalue bitscore staxids sscinames scomnames" \
    -max_target_seqs ${params.max_target_seqs} \
    -out ${sample_id}_blast_noheader.csv

  echo "qseqid,sseqid,sacc,length,pident,gapopen,gaps,evalue,bitscore,staxids,sscinames,scomnames" > ${sample_id}_blast.csv
  cat ${sample_id}_blast_noheader.csv >> ${sample_id}_blast.csv
  """
}


process csvtk_freq {

  tag { sample_id }

  executor 'local'

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_blast_freqs.csv"

  input:
  tuple val(sample_id), path(blast_results)

  output:
  tuple val(sample_id), path("${sample_id}_blast_freqs.csv")

  script:
  """
  csvtk freq -f 'scomnames' ${blast_results} | csvtk sort -k 2:Nr > ${sample_id}_blast_freqs.csv
  """
}
