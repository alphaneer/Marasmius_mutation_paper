#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
hisat2-bed2gff.py
Version 0.1
Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

This script converts the bed file produced by hisat2 with the option
--novel-splicesite-outfile to a gff format compatible with e.g. GeneMark-ET.

Copyright (c) 2019, Johannesson lab
Licensed under the MIT license. See LICENSE file.
"""


import argparse
parser = argparse.ArgumentParser(description="Converter of bed file produced by \
                                            hisat2 with the option \
                                            --novel-splicesite-outfile to a gff \
                                            file compatible with e.g. GeneMark-ET.")
parser.add_argument("bed", help="Input bed file. Required.", type = str)
parser.add_argument("-o","--output", help="Output prefix", default = "output", \
                    type = str)

args = parser.parse_args()

def main():
    with open(args.bed, "r") as bed:
        gff = []
        for line in bed:
            line = line.strip()
            fields = line.split("\t")
            contig = fields[0]

            # +2 because of a bug in hisat2 that shifts the intron starting
            # splice site
            start = int(fields[1]) + 2
            stop = fields[2]
            orient = fields[3]

            gff_line = "\t".join([contig,"hisat2", "intron", str(start),stop, \
                                ".", orient,".","."])
            gff.append(gff_line)

        gff = "\n".join(gff)

    with open(args.output + ".gff", "w") as out:
        out.write(gff)

if __name__ == "__main__":
    main()
