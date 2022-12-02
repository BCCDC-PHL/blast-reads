#!/usr/bin/env python3

import argparse
import csv
import json

def parse_counts_by_name(counts_by_name_path, name_field, count_field):
    """
    """
    counts_by_name = {}
    with open(counts_by_name_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            name = row[name_field]
            count = int(row[count_field])
            counts_by_name[name] = count

    return counts_by_name


def combine_counts_by_name(counts_by_name_1, counts_by_name_2):
    """
    """
    combined_counts_by_name = {}
    for name, count in counts_by_name_1.items():
        if name in counts_by_name_2:
            combined_counts_by_name[name] = count + counts_by_name_2[name]
        else:
            combined_counts_by_name[name] = count

    for name, count in counts_by_name_2.items():
        if name not in combined_counts_by_name:
            combined_counts_by_name[name] = count

    return combined_counts_by_name
    

def main(args):
    counts_by_name_1 = parse_counts_by_name(args.counts_by_name_1, args.name_field, args.count_field)

    counts_by_name_2 = parse_counts_by_name(args.counts_by_name_2, args.name_field, args.count_field)
    counts_by_name_combined = combine_counts_by_name(counts_by_name_1, counts_by_name_2)

    counts_by_name_unsorted = [(name, counts) for name, counts in counts_by_name_combined.items()]
    counts_by_name_sorted = [(name, counts) for name, counts in sorted(counts_by_name_unsorted, key=lambda x: x[1], reverse=True)]

    print('\t'.join([args.name_field, args.count_field]))
    for record in counts_by_name_sorted:
        print('\t'.join([record[0], str(record[1])]))
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()     
    parser.add_argument('--name-field', default="name")
    parser.add_argument('--count-field', default="count")
    parser.add_argument('counts_by_name_1')
    parser.add_argument('counts_by_name_2')
    args = parser.parse_args()
    main(args)
