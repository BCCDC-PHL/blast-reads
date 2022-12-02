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

  tag { sample_id + " / " + read_type }

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${read_type}_blast.tsv"

  input:
  tuple val(sample_id), path(reads), path(blast_db_dir), val(blast_db_name), val(read_type)

  output:
  tuple val(sample_id), path("${sample_id}_${read_type}_blast.tsv"), val(read_type)

  script:
  max_hsps = params.max_hsps == "NO_VALUE" ? "" : "-max_hsps ${params.max_hsps}"
  """
  export BLASTDB="${blast_db_dir}"
  blastn \
    -num_threads ${task.cpus} \
    -query ${reads} \
    -db ${blast_db_name} \
    -outfmt "6 qseqid sseqid sacc sstrand qstart qend sstart send length pident gapopen gaps evalue bitscore staxids sscinames scomnames" \
    -max_target_seqs ${params.max_target_seqs} \
    ${max_hsps} \
    -evalue ${params.evalue} \
    -out ${sample_id}_${read_type}_blast_noheader.tsv

  echo "qseqid,sseqid,sacc,sstrand,qstart,qend,sstart,send,length,pident,gapopen,gaps,evalue,bitscore,staxids,sscinames,scomnames" | tr ',' \$'\\t' > ${sample_id}_${read_type}_blast.tsv
  cat ${sample_id}_${read_type}_blast_noheader.tsv >> ${sample_id}_${read_type}_blast.tsv
  """
}

process choose_best_hsp_per_query {

  tag { sample_id + ' / ' + read_type}

  executor 'local'

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${read_type}_blast_best_hsp_per_query.tsv"

  input:
  tuple val(sample_id), path(blast_report), val(read_type)

  output:
  tuple val(sample_id), path("${sample_id}_${read_type}_blast_best_hsp_per_query.tsv"), val(read_type)

  script:
  """
  choose_best_hsp_per_query.py ${blast_report} > ${sample_id}_${read_type}_blast_best_hsp_per_query.tsv
  """
}


process csvtk_freq {

  tag { sample_id + " / " + read_type }

  executor 'local'

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_${read_type}_blast_counts.tsv"

  input:
  tuple val(sample_id), path(blast_results), val(read_type)

  output:
  tuple val(sample_id), path("${sample_id}_${read_type}_blast_counts.tsv")

  script:
  """
  csvtk freq -t -f 'scomnames' ${blast_results} > ${sample_id}_${read_type}_blast_counts_unsorted.tsv
  echo 'sample_id,name,count' | tr ',' \$'\\t' > ${sample_id}_${read_type}_blast_counts.tsv
  sort -t \$'\\t' -k2,2nr <(tail -qn+2 ${sample_id}_${read_type}_blast_counts_unsorted.tsv) \
    | awk -F \$'\\t' 'BEGIN {OFS=FS}; {print "${sample_id}", \$0}' >> ${sample_id}_${read_type}_blast_counts.tsv
  """
}


process combine_counts {

  tag { sample_id }

  executor 'local'

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_combined_blast_counts.tsv"

  input:
  tuple val(sample_id), path(blast_counts_r1), path(blast_counts_r2)

  output:
  tuple val(sample_id), path("${sample_id}_combined_blast_counts.tsv")

  script:
  """
  echo 'sample_id,name,count' | tr ',' \$'\\t' > ${sample_id}_combined_blast_counts.tsv
  combine_counts_by_name.py ${blast_counts_r1} ${blast_counts_r2} \
    | tail -qn+2 \
    | awk -F \$'\\t' 'BEGIN {OFS=FS}; {print "${sample_id}", \$0}' >> ${sample_id}_combined_blast_counts.tsv
  """
}

