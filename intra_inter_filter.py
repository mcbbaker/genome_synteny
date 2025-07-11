# usage: python3 intra_inter_filter.py collinearity-file same-chrom-score diff-chrom-score
# purpose: removes alignments from an MCScanX-derived collinearity file that have a score less than the score provided as input to the script for each type
#          (blocks between same chromosome and blocks between different chromosomes). blocks between different chromosomes tend to need to
#          be filtered at a higher score. 

import sys
import re

# open collinearity file
col = open(sys.argv[1],'r')
same_score = float(sys.argv[2])
diff_score = float(sys.argv[3])

align_reached = False
print_cur = False

for line in col:
    line = line.strip('\n')
    align = re.search("## Alignment (\d+): score=(.+) e_value=.+ N=\d+ (..)(\d)&(..)(\d) ", line)

    if align_reached == False:
        # still not an alignment line, print line
        if align == None:
            print(line)

        # hit first alignment line
        if align:
            align_reached = True
            
            # check if same chromosome block or different chromosome block
            if int(align.group(4)) == int(align.group(6)):
                score = same_score
            else:
                score = diff_score

            # greater than score and not same species, print alignment line and change printing bool
            if float(align.group(2)) >= score and align.group(3) != align.group(5):
                print_cur = True
                print(line)
            # less than score, keep printing bool at false
            else:
                print_cur = False

    else:
        # if a new alignment line
        if align:
            # check if same chromosome block or different chromosome block
            if int(align.group(4)) == int(align.group(6)):
                score = same_score
            else:
                score = diff_score

            if float(align.group(2)) >= score and align.group(3) != align.group(5):
                print_cur = True
                print(line)
            else:
                print_cur = False
        # if intermediate line
        else:
            if print_cur == True:
                print(line)
            else:
                pass

col.close()
