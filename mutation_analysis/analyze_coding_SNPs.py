#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
analyze_coding_SNPs.py
Version 0.1
Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

This script was used in the Marasmius mutation paper to investigate
the presence of SNPs in coding regions, and to calculate the dN/dS ratio.
In short, it reads a vcf and collects single nucleotide variants,
compares these to coding regions in a gff file and predicts amino acid
changes. Also requires the genome in a fasta file to get initial gene
sequences from. The genes in the gff file have to have entries for both gene
and CDS (e.g. maker output).

Copyright (c) 2019, Johannesson lab
Licensed under the MIT license. See LICENSE file.
"""

import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Analyzes SNPs in coding regions \
                                to predict their functions. Outputs dN/dS \
                                statistics.")
parser.add_argument("reference", \
                    help="Reference genome fasta file. Required.", \
                    type = str)
parser.add_argument("annotation", \
                    help="Genome annotation gff file. Required.", \
                    type = str)
parser.add_argument("variants", \
                    help="Variant call file. Required", \
                    type = str)
args = parser.parse_args()

def readGeneticCode():
    '''
    Creates a dict with the genetic code.
    '''
    global genetic_code
    genetic_code = {"TTT":"Phe",
                    "TCT":"Ser",
                    "TAT":"Tyr",
                    "TGT":"Cys",
                    "TTC":"Phe",
                    "TCC":"Ser",
                    "TAC":"Tyr",
                    "TGC":"Cys",
                    "TTA":"Leu",
                    "TCA":"Ser",
                    "TAA":"Ter",
                    "TGA":"Ter",
                    "TTG":"Leu",
                    "TCG":"Ser",
                    "TAG":"Ter",
                    "TGG":"Trp",
                    "CTT":"Leu",
                    "CCT":"Pro",
                    "CAT":"His",
                    "CGT":"Arg",
                    "CTC":"Leu",
                    "CCC":"Pro",
                    "CAC":"His",
                    "CGC":"Arg",
                    "CTA":"Leu",
                    "CCA":"Pro",
                    "CAA":"Gln",
                    "CGA":"Arg",
                    "CTG":"Leu",
                    "CCG":"Pro",
                    "CAG":"Gln",
                    "CGG":"Arg",
                    "ATT":"Ile",
                    "ACT":"Thr",
                    "AAT":"Asn",
                    "AGT":"Ser",
                    "ATC":"Ile",
                    "ACC":"Thr",
                    "AAC":"Asn",
                    "AGC":"Ser",
                    "ATA":"Ile",
                    "ACA":"Thr",
                    "AAA":"Lys",
                    "AGA":"Arg",
                    "ATG":"Met",
                    "ACG":"Thr",
                    "AAG":"Lys",
                    "AGG":"Arg",
                    "GTT":"Val",
                    "GCT":"Ala",
                    "GAT":"Asp",
                    "GGT":"Gly",
                    "GTC":"Val",
                    "GCC":"Ala",
                    "GAC":"Asp",
                    "GGC":"Gly",
                    "GTA":"Val",
                    "GCA":"Ala",
                    "GAA":"Glu",
                    "GGA":"Gly",
                    "GTG":"Val",
                    "GCG":"Ala",
                    "GAG":"Glu",
                    "GGG":"Gly"}


def calcdnds(syn_sites,nonsyn_sites,syn_subst,nonsyn_subst):
    '''
    Calculates the dN/dS ratio from given files.
    '''
    pS, pN = syn_subst/syn_sites, nonsyn_subst/nonsyn_sites
    dS, dN = -3/4*np.log(1-(4*pS)/3), -3/4*np.log(1-(4*pN)/3)
    dnds = dN/dS
    return dnds

def collectVariants():
    '''
    Reads the vcf to collect the single nucleotide variants in a dict.
    '''
    with open(args.variants, "r") as vcf:
        variants = {}
        for line in vcf:
            line = line.strip()
            if not line.startswith("#"):
                fields = line.split("\t")
                if len(fields[3]) == 1 and len(fields[4]) == 1:
                    variants[(fields[0],fields[1])] = (fields[3],fields[4])

        return variants

def readFasta():
    '''
    Reads the input fasta file into a dict.
    '''
    fa = {}
    with open(args.reference, "r") as fasta:
        sequence = None
        for line in fasta:
            line = line.strip()

            if line.startswith(">") and sequence == None:
                header = line[1:]
                sequence = []

            elif line.startswith(">") and sequence != None:
                # If new fasta entry, add old one to dict,
                # pick new header and reset sequence
                fa[header] = "".join(sequence)
                header = line[1:]
                sequence = []

            else:
                sequence.append(line)

        # Last passthrough won't have any new entries,
        # just add the remaining sequence
        fa[header] = "".join(sequence)

    return fa

def reverse_complement(nuclstring):
    '''
    Reverse complements input nucleotide sequence.
    '''
    rev_comped = ""
    for l in reversed(nuclstring):
        if l == "A" or l == "a":
            rev_comped += "T"
        elif l == "T" or l == "t":
            rev_comped += "A"
        elif l == "C" or l == "c":
            rev_comped += "G"
        elif l == "G" or l == "g":
            rev_comped += "C"
        elif l == "N" or l == "n":
            rev_comped += "N"
    return rev_comped

def collectTranscripts(genes,fasta):
    '''
    Collects spliced gene sequences.
    '''
    spliced_genes = {}

    for k,exons in genes.items():
        exon_sequence_list = []
        #if k[3] == "-":
        #    exons = exons[::-1]

        for ex in exons:
            seq = fasta[ k[0] ][ int(ex[0])-1:int(ex[1]) ]
            exon_sequence_list.append(seq.upper())

        spliced_seq = "".join(exon_sequence_list)

        if k[3] == "-":
            spliced_seq = reverse_complement(spliced_seq)

        spliced_genes[k] = spliced_seq

    return spliced_genes

def swapbases(variants, genome):
    '''
    Swaps bases at variant positions in the genome fasta.
    '''
    variant_genome = {}
    for k,v in variants.items():
        tig = k[0]
        coord = int(k[1])
        ref_base = v[0]
        alt_base = v[1]

        if tig not in variant_genome.keys():
            variant_genome[tig] = "".join([ genome[tig][:coord-1], \
                                            alt_base, \
                                            genome[tig][coord:] ])
        else:
            s = variant_genome[tig]
            variant_genome[tig] = "".join([ s[:coord-1], \
                                            alt_base, \
                                            s[coord:] ])

    return variant_genome

def findVariantGenes(variants, genes):
    '''
    Filters the genes dict to only include genes with variants in them.
    Also returns the variants in coding regions.
    '''
    filtered_genes = {}
    coding_variants = []

    for key in variants.keys():
        tig = key[0]
        coord = int(key[1])

        for k,v in genes.items():

            # Check if tig is the same and coord is within the gene boundaries
            if k[0] == tig:
                if int(k[1]) < coord and int(k[2]) > coord:

                    # If so, go through exons to check if variant
                    # is in an exon or intron
                    for exon in v:
                        c1, c2 = int(exon[0]), int(exon[1])

                        if c1 < coord and c2 > coord:
                            filtered_genes[k] = v
                            coding_variants.append(key)

    return filtered_genes, coding_variants

def codonsplit(dnaseq):
    '''
    Splits a DNA sequence string into a list of codons.
    '''
    return [dnaseq[i:i+3] for i in range(0, len(dnaseq), 3)]

def translate(geneseq):
    '''
    Translates a gene sequence into a protein sequence.
    '''
    protseq = []
    for codon in codonsplit(geneseq):
        protseq.append( genetic_code[codon] )
    return "".join(protseq)

def calcPossibleChanges(geneseq):
    '''
    Calculates the number of possible synonymous and nonsynonymous changes.
    in a gene sequence.
    '''
    totsyn = 0
    totnonsyn = 0
    for codon in codonsplit(geneseq):
        nonsyn = 0
        for idx,base in enumerate(codon):
            bases = ["A", "T", "C", "G"]
            bases.remove(base)
            for b in bases:
                mutated_codon = "".join( [codon[:idx], b, codon[idx+1:]] )
                if translate(mutated_codon) != translate(codon):
                    nonsyn += 1/3
        syn = 3 - nonsyn
        totsyn += syn
        totnonsyn += nonsyn

    return totsyn, totnonsyn


def collectCDS():
    '''
    Collects coding region coordinates.
    '''
    with open(args.annotation, "r") as gff:
        genes = {}
        coord1,coord2,strand = 0,0,"+"
        for line in gff:
            line = line.strip()
            if not line.startswith("#"):
                fields = line.split("\t")
                tig = fields[0]
                if fields[2] == "gene":
                    coord1,coord2,strand = fields[3],fields[4],fields[6]
                elif fields[2] == "CDS":
                    if (tig,coord1,coord2,strand) not in genes.keys():
                        genes[(tig,coord1,coord2,strand)] = [ (fields[3], \
                                                            fields[4]) ]

                    else:
                        genes[(tig,coord1,coord2,strand)].append((fields[3], \
                                                                fields[4]))

        return genes

def calcObserved(cds,variant_cds):
    '''
    Calculates the number of synonymous and nonsynonymous changes
    in two gene sequences.
    '''
    syn = 0
    nonsym = 0

    for k,v in cds.items():
        ref_seq = v
        alt_seq = variant_cds[k]
        variant_codons = codonsplit(alt_seq)

        for idx, codon in enumerate(codonsplit(ref_seq)):

            if codon != variant_codons[idx]:
                if genetic_code[codon] != genetic_code[ variant_codons[idx] ]:
                    nonsym += 1
                else:
                    syn += 1
    return syn, nonsym

def calcNonCDSGenomeSize(genome, n_coding_bases):
    '''
    Calculates the noncoding size of the genome.
    '''
    genomesize = 0
    for value in genome.values():
        genomesize += len(value)
    return genomesize - n_coding_bases

def main():
    variants = collectVariants()
    genome = readFasta()
    genes = collectCDS()
    readGeneticCode()
    filtered_genes, coding_variants = findVariantGenes(variants,genes)
    n_coding = len(coding_variants)
    n_noncoding = len(variants) - n_coding

    variant_genome = swapbases(variants, genome)
    cds = collectTranscripts(filtered_genes,genome)
    variant_cds = collectTranscripts(filtered_genes,variant_genome)

    obs_syn, obs_nonsyn = calcObserved(cds,variant_cds)

    synonymous = 0
    nonsynonymous = 0

    for gene in cds.values():
        syn, nonsyn = calcPossibleChanges(gene)
        synonymous += syn
        nonsynonymous += nonsyn

    print( "Number of variants in coding regions: {}".format(str(n_coding)))
    print( "Number of variants in non-coding regions: {}\n".format(str(n_noncoding)))

    print( "Synonymous sites: " + str(synonymous))
    print( "Nonsynonymous sites: " + str(nonsynonymous))
    print("Synonymous substitutions: " + str(obs_syn))
    print("Nonsynonymous substitutions: " + str(obs_nonsyn))
    print("dN/dS: " + str(calcdnds(synonymous,nonsynonymous,obs_syn,obs_nonsyn)))

if __name__ == "__main__":
    main()
