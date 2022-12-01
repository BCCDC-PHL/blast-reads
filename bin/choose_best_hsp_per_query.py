#!/usr/bin/env python

import argparse
import csv
import json
import sys


def parse_blast_report(blast_report_path):
    header = []
    with open(blast_report_path, 'r') as f:
        header = f.readline().strip().split()

    blast_report_by_qseqid = {}
    with open(blast_report_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            qseqid = row['qseqid']
            if qseqid in blast_report_by_qseqid:
                blast_report_by_qseqid[qseqid].append(row)
            else:
                blast_report_by_qseqid[qseqid] = [row]

    return header, blast_report_by_qseqid


def choose_best_hsp(blast_records):
    blast_records_sorted = sorted(blast_records, key=lambda x: float(x['evalue']))
    best_hsp = blast_records_sorted[0]

    return best_hsp


def main(args):
    output_fieldnames, blast_report = parse_blast_report(args.blast_report)

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='excel-tab')
    writer.writeheader()
    
    for qseqid, blast_records in blast_report.items():
        writer.writerow(choose_best_hsp(blast_records))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('blast_report')
    args = parser.parse_args()
    main(args)
