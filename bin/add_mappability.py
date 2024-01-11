#!/usr/bin/env python

# this simply parses the maf and mappability files and 
# adds a mappability column to the maf file
import sys

args = sys.argv

location = {}
with open(args[1], 'r') as f:
    lines = f.readlines()
    lines = [line.strip() for line in lines]
    lines = [line.split('\t') for line in lines]
    for line in lines:
        location[line[0]] = line[4]

out_file = args[2].split('.')[0] + '.adj.maf'
with open(out_file, 'w') as f:
    with open(args[2], 'r') as g:
        for number,line in enumerate(g, start=1):
            line = line.strip()
            if number == 1:
                f.write(line + '\n')
            if number == 2:
                f.write(line + '\t' + 'mappability' + '\n')
            if number > 2:
                line_list = line.split('\t')
                maf_loc = line_list[4] + '_' + line_list[5]
                if maf_loc in location:
                    f.write(line + '\t' + location[maf_loc] + '\n')
                else:
                    f.write(line + '\t' + 'NA' + '\n')
