#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-02-16 14:59 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
coordToInt.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: convert BEDPE coords to integers based on concatenated chromosomes")
    parser.add_argument('-g', '--genome',
                        metavar='FILE', type=argparse.FileType('r'),
                        required=True,
                        help='tab delimited genome file of chrom<tab>chromSize')
    parser.add_argument('input', nargs='?',
                        type=argparse.FileType('r'),
                        default=None,
                        help='BEDPE file to read. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin

    # send back the user input
    return args

# primary function
def coord_to_int(in_file, genome):
    # running sum of chrom lengths
    runsum = 0
    # start number of each chrom
    chrom_starts = dict()
    
    for line in genome:
        v = line.rstrip().split("\t")
        runsum += int(v[1])
        chrom_starts[v[0]] = runsum

    for line in in_file:
        v = line.rstrip().split("\t")
        chrom_a = v[0]
        start_a = int(v[1])
        end_a = int(v[2])
        chrom_b = v[3]
        start_b = int(v[4])
        end_b = int(v[5])
        id = v[6]
        misc = v[7:]

        start_a = chrom_starts[chrom_a] + start_a
        end_a = chrom_starts[chrom_a] + end_a

        start_b = chrom_starts[chrom_b] + start_b
        end_b = chrom_starts[chrom_b] + end_b

        print "\t".join(map(str,
                            [(start_a + end_a) / 2,
                             (start_b + end_b) / 2,
                             id]
                            + misc))

    # close the input file
    in_file.close()
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    coord_to_int(args.input, args.genome)

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
