#!/usr/bin/env python
# -*- coding: utf-8 -*-

# genomic_mutation_distribution_chi2.py
# Version 0.1
# Author: Markus Hiltunen
# E-mail: markus.hiltunen@ebc.uu.se
#
# This script was used in the Marasmius mutation paper to investigate
# the genomic distribution of mutations. It divides the genome into bins
# of size -b, sums of the number of mutations in each bin, and
# compares all bins to each other with a chi2 test.
# Also capable of producing a simple plot (not used in the paper).
#
# LICENSING


import argparse
from scipy.stats import chisquare
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
matplotlib.pyplot.switch_backend('agg')

parser = argparse.ArgumentParser("For testing the genomic distribution of \
                                mutations with a chi2 test")
parser.add_argument("-i","--input", \
                    help="Input tab delimited file of mutation locations, e.g. \
                    vcf format", \
                    type = str)
parser.add_argument("-g","--genome", \
                    help="Input tab delimited file of scaffolds and their \
                    lengths. Fasta index .fai works", \
                    type = str)
parser.add_argument("-b","--binsize", \
                    help="Size of bins to split the genome into and use \
                    for chi2 test [100000]", \
                    default = 100000, \
                    type = int)
parser.add_argument("-p","--plot_mutations", \
                    help="Add to plot mutation bins in histograms.", \
                    action='store_true')
args = parser.parse_args()

def plotBins(bins):
    '''
    Takes the filled bins and plots a histogram.
    '''
    scaffolds = set([x.split(":")[0] for x in list(bins)])

    # Create a subplot for every scaffold
    for sc in scaffolds:
        dat = bins.loc[:,[sc in s for s in bins.columns]]
        x = np.array(dat.values.tolist()[0])
        x[x==0] = 'nan'

        y = [n for n in range(0,len(x))]
        plt.scatter(y,x, c = 'g')
        xvals = x[np.logical_not(np.isnan(x))]


        plt.title(sc)
        plt.xlabel("Coordinate (kb)")
        plt.ylabel("# mutations/{} kb".format(str(int(args.binsize/1000))))

        xlabs = [str(int(val*args.binsize/1000)) for val in y]
        plt.xticks(y, xlabs, rotation = 90)
        plt.margins(0.2)
        plt.subplots_adjust(bottom=0.15)
        plt.ylim(0,np.amax(xvals)+1)
        plt.xlim(0,np.amax(y)+1)

        #plt.show()
        plt.savefig("test.pdf")
        plt.clf()

def readTabFile(f):
    '''
    Reads the tab delimited file of genome scaffolds and coordinates.
    '''
    with open(f, "r") as tabfile:
        results = {}

        for line in tabfile:
            line = line.strip()
            f1 = line.split("\t")[0]
            f2 = int(line.split("\t")[1])

            results[f1] = f2

    return results

def readMutationFile(f):
    '''
    Reads the tab delimited file of variants provided.
    '''
    with open(f, "r") as tabfile:
        results = []

        for line in tabfile:
            # If using vcf format, get rid of header lines
            if not line.startswith("#"):
                line = line.strip()
                f1 = line.split("\t")[0]
                f2 = int(line.split("\t")[1])

                results.append((f1,f2))

    return results

def createBins():
    '''
    Creates bins from the genome.
    '''
    genome = readTabFile(args.genome)
    bins = {}

    for k,v in genome.items():
        binstart = 1

        if args.binsize > v:
            continue

        while binstart < v:

            if k not in bins.keys():
                bins[k] = [(binstart, binstart + args.binsize-1)]

            else:
                l = bins[k]
                l.append((binstart, binstart + args.binsize-1))
                bins[k] = l

            binstart += args.binsize

    return bins

def binMutations(bins):
    '''
    Puts the mutations from the given file into each bin
    '''
    mutations = readMutationFile(args.input)

    binlist = []
    for k,v in bins.items():
        for val in v:
            binlist.append("{}:{}-{}".format(k,val[0],val[1]))

    df = pd.DataFrame(np.zeros((1,len(binlist))), columns=binlist)

    for m in mutations:

        if m[0] not in bins.keys():
            continue

        for b in bins[m[0]]:

            if m[1] >= b[0] and m[1] < b[1]:
                df["{}:{}-{}".format(m[0], b[0], b[1])] += 1
                break

    return df


def main():
    print("Creating bins")
    bins = createBins()

    print("Putting mutations into bins")
    mutbins = binMutations(bins)

    print("Calculating chi2")
    chi2 = chisquare(mutbins.values.tolist()[0])
    print("\n{}\n".format(chi2))

    if args.plot_mutations:
        print("Plotting histogram")
        plotBins(mutbins)



if __name__ == "__main__":
    main()
