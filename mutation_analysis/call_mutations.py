#!/usr/bin/env python
# -*- coding: utf-8 -*-

# call_mutations.py
# Version 0.1
# Author: Markus Hiltunen
# E-mail: markus.hiltunen@ebc.uu.se
#
# This script was used in the Marasmius mutation paper to filter out
# possible mutations from other variants in a variant table produced by
# GATK variantsToTable. The principle is to look for variants that are present
# in a subset of sampling points in a given fairy ring, e.g. present in
# sample 01N but absent in 01E, 01S, etc.
# Variants that conform to this principle are further filtered, as there will be
# LOTS of false positives otherwise to look through manually in e.g. IGV.
# Filtering is done in two steps: coverage (options -n and -x for minimum and
# maximum respectively) and read frequency (option f).
#
# LICENSING

from operator import truediv # To be able to divide lists
from scipy import stats as ss
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", \
                    help="Input table", \
                    type = str)
parser.add_argument("-o","--output", \
                    help="Output prefix", \
                    type = str)
parser.add_argument("-n","--minimum_coverage", \
                    help="Minimum coverage to keep a variant [20]", \
                    type = int, \
                    default = 20)
parser.add_argument("-x","--maximum_coverage", \
                    help="Maximum coverage to keep a variant [60]", \
                    type = int, \
                    default = 60)
parser.add_argument("-f","--frequency_cutoff", \
                    help="Minimum variant read frequency to keep a variant [0.2]", \
                    type = float, \
                    default = 0.2)
parser.add_argument("-s","--skip_indels", \
                    help="Add to skip indels (useful for pileup data).", \
                    action='store_true')
args = parser.parse_args()

def getOut():
    '''
    To create a prefix for output files.
    '''
    if args.output:
        outprefix = args.output
    else:
        outprefix = "potential_mutations"
    return outprefix

def test_signif(var,cov):
    '''
    Tests for significant skew of reads towards one variant
    Allows one error per 100 000 variants
    True means site passed check (= is not significantly different
    from a 0.5 frequency)
    '''
    alpha = 0.00001
    p = ss.binom_test(var,cov,1.0/2,alternative="two-sided")
    if p < alpha:
        return False
    else:
        return True

def main():
    rings = [ "01", "03", "04", "05", "06", "08" ]

    for ring in rings:
        outfile = getOut() + "_ring" + ring + ".table"
        outlines = []
        print("Starting ring {}".format(ring))

        with open(args.input, "r", encoding="utf-8") as input:
            for line in input:
                line = line.strip()
                fields = line.split("\t")
                fields = list(filter(None, fields))

                # If skipping indels
                if args.skip_indels == True and len(fields[3]) != len(fields[4]):
                    continue

                if fields[0] == "CHROM":
                    n = 0
                    samples = []
                    genotype_positions = []
                    for i in fields:
                        if ring in i:
                            samples.append(fields[n])
                            genotype_positions.append(n)
                        n += 1

                    outlines.append("\t".join(fields[0:6])+"\t"+"\t".join(samples))
                    #print("\t".join(fields[0:6])+"\t"+"\t".join(samples))
                    continue

                qual = float(fields[5])
                pos = int(fields[1])
                tig = fields[0]
                minpos = min(genotype_positions)
                maxpos = max(genotype_positions) + 1
                genotypes = fields[minpos:maxpos]

                # Convert to list of tuples
                genotypes = [ (int(i.split(",")[0]), int(i.split(",")[1])) \
                            if i != "NA" else (0,0) for i in genotypes ]
                variant_alleles, ref_alleles, cov = [], [], []

                for g in genotypes:
                    variant_alleles.append(g[1])
                    ref_alleles.append(g[0])
                    cov.append(sum(g))

                # Divides two lists element if cov > 0
                variant_frequencies = [ai/bi if bi > 0 else ai for ai,bi in zip(variant_alleles,cov)]

                # If all frequencies are the same, continue
                if all(x == variant_frequencies[0] for x in variant_frequencies):
                    continue

                # Look for homozygous sites
                if 1.0 in variant_frequencies or 0.0 in variant_frequencies:
                    counter = 0
                    homoz_checker, heteroz_checker, minor_allele_checker, freq_checker = \
                    False, False, False, False
                    totcov = 0
                    tot_ref, tot_alt = 0, 0

                    for i in variant_frequencies:

                        # Check total ref/alt frequency at this site,
                        # for heterozygous samples with more than 1 read
                        # containing variant allele
                        if i != 1.0 and i != 0.0:
                            totcov = totcov + cov[counter]

                            # Add number of variant reads to total variant
                            # allele counter
                            tot_alt += variant_alleles[counter]
                            tot_ref += ref_alleles[counter]

                            if i > args.frequency_cutoff \
                            and i < 1 - args.frequency_cutoff:
                                freq_checker = True

                        # Check that at least one homozygous site has
                        # sufficient coverage to consider
                        if i == 0 or i == 1:
                            if int(cov[counter]) >= args.minimum_coverage \
                            and int(cov[counter]) <= args.maximum_coverage:
                                #print(int(n_o_reads[counter]))
                                homoz_checker = True

                        # Check that heterozygous site has at least
                        # args.coverage_cutoff coverage
                        elif int(cov[counter]) >= args.minimum_coverage \
                        and int(cov[counter]) <= args.maximum_coverage:
                            heteroz_checker = True

                        counter += 1

                    # If checkers are true, continue to check minor_allele
                    # read frequency
                    if homoz_checker == True \
                    and heteroz_checker == True \
                    and freq_checker == True:

                        # False positives often have several samples with a
                        # low frequency of the variant in the reads.
                        # Real variants, not caused by sporadic read mapping,
                        # should follow a binomial distribution of frequency
                        # in the reads.
                        # This is solved by calculating a p-value for the
                        # read frequency across all samples where the
                        # variant is found. Samples with a single read
                        # containing the variant are ignored as this could be
                        # the result of sequencing error or index hopping
                        minor_allele = tot_alt if tot_alt < tot_ref else tot_ref
                        if minor_allele > 1:
                            if test_signif(minor_allele, totcov) == True:

                                # Convert genotypes back to string
                                genotypes = [ "{},{}".format(i[0],i[1]) for i in genotypes ]
                                outlines.append("\t".join(fields[0:6])+"\t"+"\t".join(genotypes))
                                #print("\t".join(fields[0:6])+"\t"+"\t".join(genotypes))

        with open(outfile, "w") as out:
            out.write("\n".join(outlines))

if __name__ == "__main__":
    main()
