#!/usr/bin/env python

############################################################
# filter_platypus.py Daniel Nichol 01/05/2018 George Cresswell 12/09/2018
#
# Filters the vcf file produced by platypus
#
# The filtering criteria are:
#   1) The FILTER tag for the call is in filterList (below)
#   2) The variant is not a known GL variant, or aligned to the decoy genome.
#   3) A genotype is called for all samples (i.e. no sample assigned './.')
#   4) The genotype phred scores (GQ) are >=10 for all samples
#   5) The number of reads covering the variant site is >=10 in all samples.
#   6) The tumour has a 'tef' times enrichment in VAF from in the tumour with the
#      largest VAF for the mutation compared to the normal
#   7) At least one tumour samples has >=3 reads for the variant.
#
############################################################
import sys
import pandas as pd

#####################################################################
# Filter list - the list of FILTER entries that are acceptable.
# Amend this list to include/exclude specific FILTER values.
#####################################################################

# We are going to remove list filtering for now, because it seems we are only filtering out the "strandbias" and "MQ" filter
#  ...so instead we will track how many come through with these flags

# filterList =   ['PASS', 'alleleBias', 'Q20', 'Q20;alleleBias','QD', 'Q20;QD', 'QD;alleleBias',
#                 'Q20;QD;alleleBias', 'SC', 'SC;Q20', 'SC;alleleBias', 'SC;Q20;alleleBias', 'SC;QD',
#                 'SC;Q20;QD', 'SC;QD;alleleBias', 'SC;Q20;QD;alleleBias', 'HapScore', 'Q20;HapScore',
#                 'HapScore;alleleBias', 'Q20;HapScore;alleleBias', 'QD;HapScore', 'Q20;HapScore;QD',
#                 'QD;HapScore;alleleBias', 'Q20;HapScore;QD;alleleBias', 'SC;HapScore', 'SC;Q20;HapScore',
#                 'SC;HapScore;alleleBias', 'SC;Q20;HapScore;alleleBias', 'SC;HapScore;QD', 'SC;Q20;HapScore;QD',
#                 'SC;HapScore;QD;alleleBias', 'SC;Q20;HapScore;QD;alleleBias',
#                 # These are experimental pass flag and IS NOT comprehesive and needs to made so
#                 'HapScore;badReads;QD;alleleBias',
#                 'Q20;HapScore;badReads;QD;alleleBias',
#                 'Q20;HapScore;badReads;QD',
#                 'SC;Q20;badReads;alleleBias',
#                 'HapScore;badReads;alleleBias',
#                 'SC;badReads;QD',
#                 'SC;badReads',
#                 'SC;Q20;badReads',
#                 'badReads;QD',
#                 'SC;Q20;badReads;QD;alleleBias',
#                 'Q20;badReads;alleleBias',
#                 'SC;badReads;QD;alleleBias',
#                 'SC;badReads;alleleBias',
#                 'SC;Q20;badReads;QD',
#                 'badReads',
#                 'Q20;badReads',
#                 'Q20;badReads;QD;alleleBias',
#                 'badReads;QD;alleleBias',
#                 'badReads;alleleBias',
#                 'Q20;badReads;QD']


#####################################################################
# Parse the arguments
#####################################################################
if len(sys.argv) != 4:
    print("python filter_platypus.py <input> <normal_name> <tumour_enrichment>")
    exit()
elif len(sys.argv) == 4:
    platypus_file = sys.argv[1]
    normal_name = sys.argv[2]
    tef = float(sys.argv[3])

#####################################################################
# Filter the platypus file
#####################################################################
with open(platypus_file[0:-4]+"_filtered.vcf", 'w') as platypus_pass:
    with open(platypus_file, 'r') as platcalls:

        # Copy the header of the platypus file.
        line = ''
        while line[0:6] != "#CHROM":
            line = platcalls.readline()
            platypus_pass.write(line)

        # Identify the normal sample location in the headers
        header = line[:-1].split('\t')[9:]
        samples = len(header)
        normIx = header.index(normal_name)

        line = platcalls.readline()
        while line:
            record = line[:-1].split('\t')
            # Split the (colon (:) separated) individual sample information for each sample.
            for k in range(len(record[:-samples]), len(record)):
                record[k] = record[k].split(':')
                # If multiple alts exist, we take NR, NV to be the maximum amongst them
                if ',' in record[k][4]:
                    record[k][4] = max(record[k][4].split(','))
                    record[k][5] = max(record[k][5].split(','))

            allSamples = record[-samples:]
            tumSamples = [allSamples[i]
                          for i in range(len(allSamples)) if i != normIx]
            normSample = allSamples[normIx]

            # Filters: FILTER passed and not decoy/germline
            # if (record[6] in filterList) and (record[0]!='hs37d5') and (record[0][0:2]!='GL'):

            # Filter here by searching for the string of FAIL flags, e.g. MQ and strandBias
            MQ_fail = pd.Series(record[6]).str.contains('MQ').tolist()
            SB_fail = pd.Series(record[6]).str.contains('strandBias').tolist()

            if MQ_fail[0]:
                print("Mapping Quality failure..."+record[0]+":"+record[1])

            if SB_fail[0]:
                print("Strand Bias failure..."+record[0]+":"+record[1])

            if not MQ_fail[0] and not SB_fail[0] and (record[0] != 'hs37d5') and (record[0][0:2] != 'GL'):
                # Change to at least one sample having a genotype quality >=60 (i.e. 1 in a million change of being wrong)
                GQsPass = any([float(el[3]) >= 60 for el in allSamples])
                NRsPass = all([float(el[4]) >= 10 for el in allSamples])

                # GQs (genotype phred scores) all >= 10, and
                # NRs (number of reads at site) all >= 10
                if GQsPass and NRsPass:
                    allGTs = [el[0] for el in allSamples]
                    tumGTs = [el[0] for el in tumSamples]
                    normGT = normSample[0]

                    # If the genotype for all samples exist:
                    if ("./." not in allGTs) and (normGT == "0/0") and (tumGTs.count("0/0") != (samples-1)):

                        # What is the maximum vaf of the variant in the tumour samples
                        tvafs = [float(el[5]) / float(el[4])
                                 for el in tumSamples]
                        mxtvaf = max(tvafs)

                        # What is the normal vaf?
                        nvaf = float(normSample[5]) / float(normSample[4])

                        # Enrichment compared to normal
                        # >=3 variant reads for at least one tumour sample
                        if float(normSample[5]) == 0:
                            normNVpass = True
                        else:
                            normNVpass = (mxtvaf / nvaf) >= tef
                            # Some print lines for trouble shooting
                            # print("TEF filter found passable mutation")
                            # print(normNVpass)
                            # print(mxtvaf)
                            # print(nvaf)
                            # print(line)
                        tumNVpass = any(
                            [float(el[5]) >= 3 for el in tumSamples])

                        if normNVpass and tumNVpass:
                            platypus_pass.write(line)

            line = platcalls.readline()
