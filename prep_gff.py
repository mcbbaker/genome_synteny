#!/usr/bin/env python3

# purpose: prep file that contains gff3 info (grepped for 'gene') to put in the format to run mcscanx
# usage: python3 prep_gff.py gff-file sp-name (ie 'lc', NOT with a number)

import sys
import re

gff = open(sys.argv[1],'r')
sp = sys.argv[2]

for line in gff:
    line = line.strip('\n')
    line = line.split('\t')
    chrom = re.search("^.+\.Chr(\d)$",line[0])
    # only keep Chrom entries
    if chrom:
        name = re.search("^ID=(.+);Name=",line[8])
        print(sp + chrom.group(1), name.group(1), line[3], line[4], sep='\t')
    # don't print out entries from unitigs
    else:
        pass

gff.close()
