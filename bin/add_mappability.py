#!/usr/bin/env python

# this simply parses the maf and mappability files and 
# adds a mappability column to the maf file
import sys

args = sys.argv
print(args)

location = {}
with open(args[1], 'r') as f:
    lines = f.readlines()
    lines = [line.strip() for line in lines]
    lines = [line.split('\t') for line in lines]
    for line in lines:
        my_key = line[0].split("_")
        my_key = [my_key[0], int(my_key[1])]
        # the key is the location
        my_key = tuple(my_key)
        # the value is the mappability
        location[my_key] = line[4]

out_file = args[2].split('.')[0] + '.adj.maf'

found_location = []
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
                # this ugly if else series allows matching near the change
                maf_loc = [line_list[4],line_list[5]]
                if tuple([maf_loc[0],maf_loc[1]]) in location:
                    f.write(line + '\t' + location[tuple([maf_loc[0],maf_loc[1]])] + '\n')
                elif tuple([maf_loc[0],int(maf_loc[1]) + 1]) in location:
                    f.write(line + '\t' + location[tuple([maf_loc[0],int(maf_loc[1]) + 1])] + '\n')
                elif tuple([maf_loc[0],int(maf_loc[1]) - 1]) in location:
                    f.write(line + '\t' + location[tuple([maf_loc[0],int(maf_loc[1]) - 1])] + '\n')
                else:
                    f.write(line + '\t' + 'NA' + '\n')