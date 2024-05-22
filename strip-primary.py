#!/usr/bin/env python

# purpose: remove all blast hits that contain a query or subject with a non-primary transcript (only hits with both the query and subject
# being the primary transcript .1 are retained). also strip the transcript from query and subject to conform with mcscanx criteria.  
# usage: python3 strip-primary.py blast-hits

import sys
import re

output = open(sys.argv[1],'r')

for line in output:
    line = line.strip('\n')
    line = line.split('\t')
    # start of string followed by at least 1 character, ending with a period followed by at least one digit and then end of the string
    query = re.search("^(.+)\.(\d+)$",line[0])
    subject = re.search("^(.+)\.(\d+)$",line[1])
    # strip transcript
    line[0] = query.group(1)
    line[1] = subject.group(1)
    if (int(query.group(2)) != 1 or int(subject.group(2)) != 1):
        pass
    else:
        print(*line,sep='\t')


output.close()
