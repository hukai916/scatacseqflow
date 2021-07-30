#!/usr/bin/env python

import sys
from collections import OrderedDict

filename = sys.argv[1]

dict = OrderedDict()

for line in open(filename):
    gene_id = line.split("\t")[8].split(";")[0]
    if not gene_id in dict:
        dict[gene_id] = [line]
    else:
        dict[gene_id].append(line)

for key in dict:
    for x in dict[key]:
        print(x, end="")
