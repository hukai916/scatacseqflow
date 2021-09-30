#!/usr/bin/env python

"""
Given full whitelist and index fastq, determine the valid barcode pool and its frequency.

Usage:
python get_barcode_pool.py whitelist.txt.gz index.fastq.gz outfile.txt

"""

import sys
import gzip
import pysam
import numpy as np

whitelist_file = sys.argv[1]
index_fastq    = sys.argv[2]
outfile_name   = sys.argv[3]

# read in whitelist and make a dict:
whitelist_dict = {}
if whitelist_file.endswith(".gz"):
    with gzip.open(whitelist_file, "rt") as f:
        for line in f:
            if not line.strip() in whitelist_dict:
                whitelist_dict[line.strip()] = 0
else:
    with open(whitelist_file) as f:
        for line in f:
            if not line.strip() in whitelist_dict:
                whitelist_dict[line.strip()] = 0

# count barcode from index fastq file:
not_in_whitelist = 0 # count how many barcodes are not in whitelist
with pysam.FastxFile(index_fastq) as f:
    for read in f:
        if read.sequence in whitelist_dict:
            whitelist_dict[read.sequence] += 1
        else:
            not_in_whitelist += 1

total_valid_barcode = sum(whitelist_dict.values())
total_unique_valid_barcode = np.count_nonzero(list(whitelist_dict.values()))

# store output:
if outfile_name.endswith(".gz"):
    out_file = gzip.open(outfile_name, "wt")
else:
    out_file = open(outfile_name, "wt")

for barcode in whitelist_dict:
    if not whitelist_dict[barcode] == 0:
        out_file.write("\t".join([barcode, str(whitelist_dict[barcode] / total_valid_barcode)]) + "\n")
out_file.close()

# output some statistics:
print("total valid bacodes: ", total_valid_barcode)
print("total unique valid barcodes: ", total_unique_valid_barcode)
print("not_in_whitelist / total valid barcodes: ", not_in_whitelist/total_valid_barcode)
print("If above value > 0.25, you may have used the wrong whitelist.")
print("Done!")
