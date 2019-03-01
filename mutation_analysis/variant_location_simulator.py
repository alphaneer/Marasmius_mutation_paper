#!/usr/bin/env python
# -*- coding: utf-8 -*-

# variant_location_simulator.py
# Version 0.1
# Author: Markus Hiltunen
# E-mail: markus.hiltunen@ebc.uu.se
#
# This script randomizes genomic locations for BAMsurgeon to simulate variants
# at.
#
# LICENSING

import argparse
import random
import pysam

parser = argparse.ArgumentParser(description="Takes a reference genome and simulates SNV positions in a bed output.")
parser.add_argument("-r", "--reference", help="Fasta genome file. Required.", type = str)
parser.add_argument("-d", "--variant_density", help="Simulate an SNV on average every nth bp [100000]", default = 100000, type = int)
parser.add_argument("-o", "--out", help="Prefix for output.", type = str)
parser.add_argument("-i","--indels", help="Add to produce indels instead of SNPs.", action='store_true')
parser.add_argument("-f", "--indelrate", help="Insertion/deletion ratio [0.33]", default = 0.33, type = float)
args = parser.parse_args()

def simindels(nuclstr, positions):
    '''
    Takes a nucleotide string and list of positions to simulate indels at
    these positions
    '''
    indels = [] # To return data from

    # Randomize if insertion or deletion
    delchance = 100-int(args.indelrate*100)
    inchance = int(args.indelrate*100)

    for pos in positions:
        type = random.choices(['in', 'del'], [inchance, delchance])[0]
        VAF = random.uniform(0.3,0.7)

        # Randomize length of indel, max 10 bp
        length = random.randrange(1,10)

        if type == "del":
            indels.append("{}\t{}\t{}\t{}".format(pos, pos+length, VAF, "DEL"))

        # Simulate tandem duplications
        else:
            basestoinsert = nuclstr[pos-length:pos]

            if len(basestoinsert) > 0:
                indels.append("{}\t{}\t{}\t{}\t{}".format(pos, pos+1, VAF, "INS", basestoinsert))
    return indels

def simloc(nuclstr):
    '''
    Takes a nucleotide string and simulates variant locations on average every
    <args.variant_density> bp
    '''
    locations = []
    n = 0
    while n < len(nuclstr):
        loc = random.randint(n, n+args.variant_density)
        locations.append(loc)
        n += args.variant_density
    return locations

def getOut():
    if args.out:
        outprefix = args.out
    else:
        if "/" in args.reference:
            outprefix = args.reference.split("/")[-1].split(".fasta")[0]
        else:
            outprefix = args.reference.split(".fasta")[0]
    return outprefix


def main():
    print("\nReading {}...".format(args.reference))
    ref = pysam.FastaFile(args.reference)
    print("Simulating variants from {} on average every {} kb...\n".format(args.reference, str(args.variant_density)))

    outlist = []
    for k in ref.references:
        seq = ref.fetch(k)
        variants = simloc(seq) # Returns a list of variant positions

        if args.indels == True:
            for indel in simindels(seq, variants):
                outlist.append("{}\t{}".format(k,indel))

        else:
            for var in variants:
                VAF = random.gauss(0.4374843111024255, 0.13189590475114896) # Numbers from empirical data (newasm calls)
                outlist.append("{}\t{}\t{}\t{}".format(k,var,var,VAF))
    ref.close()

    outprefix = getOut()

    with open(outprefix+".bed", "w") as out:
        out.write("\n".join(outlist))


if __name__ == "__main__":
    main()
