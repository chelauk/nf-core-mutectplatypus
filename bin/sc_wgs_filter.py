#!/usr/bin/env python

############################################################
#  Minimal filter for sc_wgs 
#
# Filters the VCF file produced by Platypus
#
# The filtering criteria are:
#   1) The variant is not a known GL variant.
#
############################################################

import sys
import argparse

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Filter Platypus VCF file")
    parser.add_argument("input_file", help="Input VCF file produced by Platypus")
    parser.add_argument("normal_name", help="Name of the normal sample")
    
    args = parser.parse_args()
    
    # Open the input VCF file and the output filtered VCF file
    input_file = args.input_file
    normal_name = args.normal_name
    
    if not input_file.endswith(".vcf"):
        print("Error: Input file should be a VCF file")
        sys.exit(1)
    
    output_file = input_file.replace(".vcf", "_filtered.vcf")
    
    with open(output_file, 'w') as platypus_pass:
        with open(input_file, 'r') as platcalls:
            # Copy the header of the Platypus file
            line = ''
            while line[:6] != "#CHROM":
                line = platcalls.readline()
                platypus_pass.write(line)
            
            # Identify the normal sample location in the headers
            header = line.strip().split('\t')[9:]
            if normal_name not in header:
                print(f"Error: Normal sample '{normal_name}' not found in the VCF header")
                sys.exit(1)
                
            normIx = header.index(normal_name)
            
            # Process each record
            line = platcalls.readline()
            while line:
                record = line.strip().split('\t')
                
                # Split the individual sample information for each sample
                for k in range(len(record) - len(header), len(record)):
                    record[k] = record[k].split(':')
                    if ',' in record[k][4]:
                        record[k][4] = max(record[k][4].split(','))
                        record[k][5] = max(record[k][5].split(','))
                
                allSamples = record[-len(header):]
                normSample = allSamples[normIx]
                
                # Filter: not germline
                if not record[0].startswith('GL'):
                    platypus_pass.write(line)
                
if __name__ == "__main__":
    main()