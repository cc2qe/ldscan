#!/usr/bin/env python

import sys
import getopt
import string
from optparse import OptionParser
from collections import Counter
import random

def filter_vcf(data, sampleFile, maxHetRate, minMAF, maxMAF, AFD):
    minAF = float(minMAF)
    maxAF = float(maxMAF)
    AFD = float(AFD)
    maxHetRate = float(maxHetRate)

    # map from VCF column to sample name
    sample_to_col = dict()
    
# load the sample data
    file = open(sampleFile,'r')
    sample = []
    for line in file:
        sample.append(line.strip().split('\t'))

# get vectors which are the indices of the samples for each of ASN, AMR, AFR & EUR superpopulations;
    amr = set()
    afr = set() 
    asn = set()
    eur = set()

    for i in xrange(len(sample)):
        if sample[i][2] == 'AMR':
            amr.add(sample[i][0]);
        elif sample[i][2] == 'AFR':
            afr.add(sample[i][0]);
        elif sample[i][2] == 'ASN':
            asn.add(sample[i][0]);
        elif sample[i][2] == 'EUR':
            eur.add(sample[i][0])

    # load the snp data
    snp = []
    for line in data:
        if line[0] == "#":
            if line[:6] == "#CHROM":
                v = line.rstrip().split('\t')
                for i in xrange(9,len(v)):
                    sample_to_col[v[i]] = i
            continue

        # split the line by tabs
        v = line.rstrip().split('\t')

        # number of heterozygotes at locus.
        numHet = 0

        # the number of samples that are not ./. at the locus
        numInformative = 0;
        for i in range(9,len(v)):
            gt = v[i].split(":")[0]
            print gt
            if gt == "./.":
                v[i] = -1
            else:
                v[i] = gt.replace(v[3],'0')
                v[i] = gt.replace(v[4],'1')
                v[i] = str(int(v[i][0]) + int(v[i][2]))
                numInformative += 1

            if v[i] == '1':
                numHet += 1

        gen = v[9:]
        print gen
        maxHetCount = maxHetRate * numInformative

        # now calculate the allele frequencies for the total set, and for each subpopulation
        # check with original VCF file to make sure that we get the same numbers
        # first get the total AF of snp1
        tCount = float(0)
        for i in range(0,len(gen)):
            if gen[i] != -1:
                tCount = tCount + float(gen[i])
        tAF = tCount / (2*numInformative)

        # now calculate the total variance:
        # how to calculate variance: 
        # 1) first calculate the mean M
        # 2) for each value Vi, Di = (Vi-M)^2
        # 3) variance is the mean of the D values
                # first the mean:
                # tMean = tCount / numSamples
                # now the difference:
                #t = 0
                #for i in range(0,len(gen)):
                #    t = t + ((float(gen[i]) - tMean)**2)
                #Vt = t / numInformative;

        # now get af for subpopulations
        eurCount = float(0)
        # number of informative eur at that locus
        eurInf = 0
        for sample in eur:
            val = sample_to_col[sample]
            if gen[val] != -1:
                eurCount = eurCount + float(gen[val])
                eurInf += 1
        if eurInf != 0:
            eurAF = eurCount / (2*eurInf)
        else: eurAF = -1

        afrCount = float(0)
        afrInf = 0
        for sample in afr:
            val = sample_to_col[sample]
            if gen[val] != -1:
                afrCount = afrCount + float(gen[val])
                afrInf += 1
        if afrInf != 0:
            afrAF = afrCount / (2*afrInf)
        else: afrAF = -1

        asnCount = float(0)
        asnInf = 0
        for sample in asn:
            val = sample_to_col[sample]
            if gen[val] != -1:
                asnCount = asnCount + float(gen[val])
                asnInf += 1
        if asnInf != 0:
            asnAF = asnCount / (2*asnInf)
        else: asnAF = -1

        amrCount = float(0)
        amrInf = 0
        for sample in amr:
            val = sample_to_col[sample]
            if gen[val] != -1:
                amrCount = amrCount + float(gen[val])
                amrInf += 1
        if amrInf != 0:
            amrAF = amrCount / (2*amrInf)
        else: amrAF = -1

        # now calculate variance among subpopulations
        #if tAF >= minAF and tAF <= maxAF and (abs(tAF-eurAF)/tAF) <= AFD and  (abs(tAF-afrAF)/tAF) <= AFD and (abs(tAF-asnAF)/tAF) <= AFD and (abs(tAF-amrAF)/tAF) <= AFD:
        #    print "tAF=" + str(tAF) + "\t" + "eurAF=" + str(eurAF) + "\t" + "afrAF=" + str(afrAF) + "\t" + "asnAF=" + str(asnAF) + "\t" + "amrAF=" + str(amrAF)
        # + "\t" + "totalVariance=" + str(Vt)

        if numHet > maxHetCount:
            continue

        if (AFD != None):
            if (tAF < minAF or
                tAF > maxAF or
                (abs(tAF-eurAF)/tAF) > AFD or
                (abs(tAF-afrAF)/tAF) > AFD or
                (abs(tAF-asnAF)/tAF) > AFD or
                (abs(tAF-amrAF)/tAF) > AFD):
                continue

        # print output
        print '\t'.join(map(str, [v[0],
                                  int(v[2])-1,
                                  v[2],
                                  v[5]] +
                            v[8:]))

        # print '\t'.join(map(str, [v[0],
        #                           v[2],
        #                           v[5],
        #                           v[6],
        #                           round(tAF,2),
        #                           round(amrAF,2),
        #                           round(afrAF,2),
        #                           round(asnAF,2),
        #                           round(eurAF,2)] +
        #                     v[8:]))

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
        help="A VCF file [stdin]",
        metavar="FILE")
        
    parser.add_option("-s", "--sampleFile", dest="sampleFile", 
        help="A tab delimited file with sample, population, and superpopulation",
        metavar="FILE")        
    parser.add_option('-m', '--maxHetRate', dest='maxHetRate', default=1, help='maximum heterozygosity rate for a variant (default: 1)')
    parser.add_option('-p', '--minMAF', dest='minMAF', default=0, help='minimum minor allele frequency (default: 0')
    parser.add_option('-q', '--maxMAF', dest='maxMAF', default=1, help='maximum minor allele frequency (default: 1)')
    parser.add_option('-a', '--alleleFreqDiff', dest='alleleFreqDiff', default=None, help='maximum amount of allele frequency between subpopulations abs(overallAlleleFreq - subpopFreq) / overallAlleleFreq [None]')

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
    filter_vcf(data, opts.sampleFile, opts.maxHetRate, opts.minMAF, opts.maxMAF, opts.alleleFreqDiff)

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
