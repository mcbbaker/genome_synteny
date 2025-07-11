#!/usr/bin/env python3

# purpose: remove blast hits that contain a gene id (given in a file of gene ids to remove)
# usage: python3 remove_blacklist_hits.py blast-file gene-ids

import sys

blast = open(sys.argv[1],'r')
remove_genes = open(sys.argv[2],'r')

remove_list = []
for line in remove_genes:
    line = line.strip('\n')
    remove_list.append(line)

for line in blast:
    line = line.strip('\n')
    line = line.split('\t')
    # query or target are in list of genes to remove, don't print line (remove hit from blast file)
    if (line[0] in remove_list) or (line[1] in remove_list):
        pass
    # neither query or target are in list of genes to remove, print line normally
    else:
        print(*line,sep='\t')

blast.close()
remove_genes.close()
