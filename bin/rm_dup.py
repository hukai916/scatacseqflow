#!/usr/bin/env python

"""
Given bam file, remove duplicate reads according to cell barcode and read coordinates.

Usage:
python rm_dup.py bamfile regular_expression_barcode

Note, pls quote regular_expression_barcode:
python rm_dup.py bamfile "[^:]*"

Codes adapted from https://github.com/tzhu-bio/UMI-ATAC-seq/blob/master/rm_umi_dup.py
"""

import os
import sys
import pysam
import re

bamfile = sys.argv[1]
re_barcode = sys.argv[2]

basename = os.path.basename(bamfile)
output_name = "rm_dup_" + basename

input_file = pysam.AlignmentFile(bamfile, "rb")
output_file = pysam.AlignmentFile(output_name, "wb", template=input_file)

unique_reads = {}

print("Processing ...")
for read in input_file:
    if len(unique_reads) % 100000 == 0:
        if not len(unique_reads) == 0:
            tem = int(len(unique_reads) / 100000) * 100000
            print("\t" + str(tem) + " reads processed ...")

    cell = re.search(re_barcode, read.query_name).group()

    read_id = "_".join([cell, read.reference_name, str(read.reference_start), str(read.template_length)])
    if not read_id in unique_reads:
        unique_reads[read_id] = 1
        output_file.write(read)
    else:
        unique_reads[read_id] += 1
print("Done!")

print("Result from rm_dup: total unique reads:\t\t" + str(len(unique_reads)))
print("Result from rm_dup: total duplicate reads:\t" + str(sum(unique_reads.values()) - len(unique_reads)))
