#!/usr/bin/env python

import sys
import getopt
import string
from optparse import OptionParser
from collections import Counter
import random

def filter_vcf(data, maxHetRate, minMAF, maxMAF):
    minAF = float(minMAF)
    maxAF = float(maxMAF)
    maxHetRate = float(maxHetRate)

    # load the snp data
    for line in data:
        if line[0] == "#":
            continue

        # split the line by tabs
        v = line.rstrip().split('\t')

        # number of heterozygotes at locus.
        numHet = 0

        # the number of samples that are not ./. at the locus
        numInformative = 0;
        for i in range(9,len(v)):
            gt = v[i].split(":")[0]
            if gt == "./.":
                v[i] = -1
            else:
                v[i] = gt.replace(v[3],'0')
                v[i] = gt.replace(v[4],'1')
                v[i] = str(int(v[i][0]) + int(v[i][2]))
                numInformative += 1

            if v[i] == '1':
                numHet += 1

        # store the genotypes for each sample
        gen = v[9:]
        maxHetCount = maxHetRate * numInformative

        # now calculate the allele frequencies for the total set, and for each subpopulation
        # check with original VCF file to make sure that we get the same numbers
        # first get the total AF of snp1
        tCount = float(0)
        for i in range(0,len(gen)):
            if gen[i] != -1:
                tCount = tCount + float(gen[i])
        tAF = tCount / (2*numInformative)

        # conditional filtering
        # max hets
        if numHet > maxHetCount:
            continue

        # allele freqs
        if (tAF < minAF or tAF > maxAF):
            continue
 
        # print output
        print '\t'.join(map(str, [v[0],
                                  int(v[1]) - 1,
                                  int(v[1]) - 1 + len(v[3]),
                                  v[2]] +
                            v[9:]))
        
#---------------------------------------------------------------------------
# argument parsing

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg		

def main():	
    usage = """%prog -i <file> 

filterVcf.py
Author: Ira Hall and Colby Chiang
Description: filter a VCF by variant population frequencies 
    """
    parser = OptionParser(usage)

    parser.add_option("-i", "--file", dest="inFile", 
        help="VCF file [stdin]",
        metavar="FILE")
    parser.add_option('-m', '--maxHetRate', dest='maxHetRate', default=1, help='maximum heterozygosity rate for a variant')
    parser.add_option('-p', '--minMAF', dest='minMAF', default=0, help='minimum minor allele frequency')
    parser.add_option('-q', '--maxMAF', dest='maxMAF', default=1, help='maximum minor allele frequency')
    # parser.add_option('-a', '--alleleFreqDiff', dest='alleleFreqDiff', default=None, help='maximum amount of allele frequency between subpopulations abs(overallAlleleFreq - subpopFreq) / overallAlleleFreq [None]')

    (opts, args) = parser.parse_args()

    # bail if no input and no stdin
    if opts.inFile is None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            data = sys.stdin
    else:
        data = open(opts.inFile)

    # process the fil
    filter_vcf(data, opts.maxHetRate, opts.minMAF, opts.maxMAF)

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
